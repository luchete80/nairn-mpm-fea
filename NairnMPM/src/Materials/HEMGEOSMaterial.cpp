/********************************************************************************
	HEMGEOSMaterial
	nairn-mpm-fea

	Created by John Nairn, Feb 18, 2013.
	Copyright (c) 2013 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "Materials/HEMGEOSMaterial.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Exceptions/MPMWarnings.hpp"

// class statics
int HEMGEOSMaterial::warnExcessiveX = -1;

#pragma mark HEMGEOSMaterial::Constructors and Destructors

// Constructor
HEMGEOSMaterial::HEMGEOSMaterial(char *matName,int matID) : HEIsotropic(matName,matID)
{
	gamma0=1.64;		// dimensionless
	C0=4004000.;		// mm/sec
	S1=1.35;			// dimsionless
	S2=0.;				// dimsionless
	S3=0.;				// dimsionless
	Kbulk = 1.;			// not used
	Kmax = 0;			// maxium relative increase allows in K (<1 not limit)
}

#pragma mark HEMGEOSMaterial::Initialization

// Read material properties
char *HEMGEOSMaterial::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
	//	G found in IsotropicMat
	//	rho, Cv found in material base class
	
	// here are the rest
    if(strcmp(xName,"gamma0")==0)
    {	input=DOUBLE_NUM;
        return((char *)&gamma0);
    }
    
    else if(strcmp(xName,"C0")==0)
    {	input=DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&C0,gScaling,1.e3);
    }
    
    else if(strcmp(xName,"S1")==0)
    {	input=DOUBLE_NUM;
        return((char *)&S1);
    }
	
    else if(strcmp(xName,"S2")==0)
    {	input=DOUBLE_NUM;
        return((char *)&S2);
    }
	
    else if(strcmp(xName,"S3")==0)
    {	input=DOUBLE_NUM;
        return((char *)&S3);
    }
	
    else if(strcmp(xName,"Kmax")==0)
    {	input=DOUBLE_NUM;
        return((char *)&Kmax);
    }
	
	return(HEIsotropic::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *HEMGEOSMaterial::VerifyAndLoadProperties(int np)
{
	// call plastic law first
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
    
	// check properties (need here because IsotropicMat is skipped
	if(G1<0) return "The shear modulus, G1, is missing";
	
	// G in specific units using initial rho (F/L^2 L^3/mass)
	G1sp = G1/rho;
	
    // Use in place of C0^2. Units are L^2/sec^2 = F/L^2 L^3/mass
	// Equal to reduced bulk modulus
    Ksp = C0squared = C0*C0;
	
    // Initial bulk modulus
	Kbulk = rho*C0squared;
	
	// expansion coefficients - affect on pressure is handled by eos, but
	// needed for shear parts (which in large deformation have Jres component)
	double effAlpha = (heatCapacity*gamma0)/C0squared;
	CTE1 = effAlpha/3.;
	
	// this material not coupled to moisture expansion or poroelasticity
	betaI = 0.;
	CME1 = 0.;
#ifdef POROELASTICITY
	alphaPE = 0.;
#endif
	
    // for Cp-Cv (units nJ/(g-K^2))
    Ka2sp = Ksp*CTE1*CTE1;
	
	// warning
	if(Kmax>1.)
	{	if(warnExcessiveX<0)
			warnExcessiveX = warnings.CreateWarning("Compressive strain has exceeded MGEOS law range in HEMGEOSMaterial",-1,0);
		Xmax = GetMGEOSXmax(gamma0,S1,S2,S3,Kmax);
	}
	else
		Xmax = 10.;		// i.e. never
	
	// skip Hyperelstic methods
	return MaterialBase::VerifyAndLoadProperties(np);
}

// print mechanical properties to the results
void HEMGEOSMaterial::PrintMechanicalProperties(void) const
{
	// core properties
	PrintProperty("C0",C0*UnitsController::Scaling(1.e-3),UnitsController::Label(ALTVELOCITY_UNITS));
	PrintProperty("gam0",gamma0,"");
	PrintProperty("K",Kbulk*UnitsController::Scaling(1.e-6),"");
    PrintProperty("G1",G1*UnitsController::Scaling(1.e-6),"");
	cout << endl;
    
	PrintProperty("S1",S1,"");
	PrintProperty("S2",S2,"");
	PrintProperty("S3",S3,"");
	cout << endl;
	
	PrintProperty("E",9.*Kbulk*G1/(3.*Kbulk+G1)*UnitsController::Scaling(1.e-6),"");
    PrintProperty("nu",(3.*Kbulk-2.*G1)/(6.*Kbulk+2.*G1),"");
	switch(UofJOption)
	{   case J_MINUS_1_SQUARED:
			PrintProperty("U(J)",UofJOption,"[ = (K/2)(J-1)^2 ]");
			break;
			
		case LN_J_SQUARED:
			PrintProperty("U(J)",UofJOption,"[ = (K/2)(ln J)^2 ]");
			break;
			
		case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
		default:
			PrintProperty("U(J)",UofJOption,"[ = (K/2)((1/2)(J^2-1) - ln J) ]");
			break;
	}
	cout << endl;
    
	// effective volumetric CTE (in ppm/K) alpha = rho0 gamma0 Cv / K
	PrintProperty("a",1.e6*CTE1,"");
	PrintProperty("T0",thermal.reference,"K");
	cout <<  endl;
	
	// Kmax
	if(Kmax>1.)
	{	PrintProperty("Kmax",Kmax," K0");
		PrintProperty("Xmax",Xmax,"");
	}
	else
		cout << "Kmax= no limit";
	cout <<  endl;
    
    plasticLaw->PrintYieldProperties();
	
	// skip super class, but call it's superclass
	HyperElastic::PrintMechanicalProperties();
}

// Print transport properties
void HEMGEOSMaterial::PrintTransportProperties(void) const
{
	// Conductivity constants
	if(ConductionTask::active)
	{	MaterialBase::PrintTransportProperties();
	}
	else if(!ConductionTask::adiabatic)
	{	PrintProperty("Cv",heatCapacity*UnitsController::Scaling(1.e-6),UnitsController::Label(HEATCAPACITY_UNITS));
		cout << endl;
	}
}

// if analysis not allowed, throw an exception
// throws CommonException()
void HEMGEOSMaterial::ValidateForUse(int np) const
{
	if(thermal.reference<=0)
	{	throw CommonException("MGEOSMaterial material requires the simulation to set the stress free temperature in degrees K",
							  "HEMGEOSMaterial::ValidateForUse");
	}
    
    if(np==PLANE_STRESS_MPM)
    {	throw CommonException("MGEOSMaterial material has not yet been updated to do plane stress calculations",
                              "HEMGEOSMaterial::ValidateForUse");
    }
	
	// call super class
	HEIsotropic::ValidateForUse(np);
}


#pragma mark HEMGEOSMaterial::Custom Methods

// Get plastic properties or NULL on memory error
void *HEMGEOSMaterial::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	HEPlasticProperties *p = (HEPlasticProperties *)matBuffer;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer,offset);
	// Gred and Kred found in UpdatePressure() - do not use before that
	return p;
}

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy due to dilation
// 4. Call plasticLaw to see if it wants to change the shear modulus
// J = V(T,c)/V0(T,c), Jeff = V(T,c)/V0(Tref,cref), Jres = V0(T,c)/V0(Tref,cref) = J/Jeff
// Jtot = V(T,c)/V0(Trec,cref), Jres = V0(T,c)/V0(Tref,cref) for free expansion, J = V(T,c)/V0(T,c)
// Jn+1 = (detdF/detdFres) Jn, Jresn+1 = detdFres Jresn, Jtot = detdF Jtotn
// detdFres = (1+dres)^3 (approximately)
// Here Tref and cref are starting conditions and T and c are current temperature and moisture
void HEMGEOSMaterial::UpdatePressure(MPMBase *mptr,double J,double detdF,int np,double Jeff,
									 double delTime,HEPlasticProperties *p,ResidualStrains *res,double detdFres,int offset,
									 double &dTq0,double &AVEnergy) const
{
    // J is total volume change - may need to reference to free-swelling volume if that works
	// Note that swelling looks like a problem because the sums of strains needs to be adjusted
	//		to stress-free state
	
	// previous pressure
	double P,P0 = mptr->GetPressure();
	double delV = 1. - 1./detdF;
    
    // M-G EOS
    // Want specific pressure or pressure over current density (using J = rho0/rho)
    if(J<1.)
	{	// new compression J(k+1) = 1-x(k+1)
		double x = 1.-J;
		
		if(x>Xmax)
		{	// law not valid if gets too high
			if(warnings.Issue(warnExcessiveX,-1)==GAVE_WARNING)
			{
#pragma omp critical (output)
				{	cout << "# Excessive x = " << x << " causing Kred to increase more than " << Kmax << " fold" << endl;
					mptr -> Describe();
				}
			}
			
			// reset
			x = Xmax;
		}
		
		// get reduced bulk modulus (used in AV and maybe hardening)
		p->Kred = C0squared*GetMGEOSKRatio(x,gamma0,S1,S2,S3);
		
		// compression law
		// denominator = 1 - S1*x - S2*x^2 - S3*x^3
		double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
		
		// current effective and reduced (by rho0) bulk modulus
		double Klawred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
        
		// Pressure from bulk modulus and an energy term
		double e = mptr->GetInternalEnergy();
		P = J*(Klawred*x + gamma0*e);
		
		// particle isoentropic temperature increment
		//double Kratio = p->Kred/Ksp;
		//dTq0 += -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
		dTq0 += -J*gamma0*mptr->pPreviousTemperature*delV;
    }
    else
    {	// In tension hyperelastic law (times J to get Kirchoff pressure)
		P = -J*GetVolumetricTerms(Jeff,C0squared);
		//P = -J*C0squared*(Jeff-1.);
		
		// elastic particle isentropic temperature increment
		double Kratio;				// = rho_0 K/(rho K_0)
		switch(UofJOption)
		{   case J_MINUS_1_SQUARED:
				Kratio = Jeff;
				break;
				
			case LN_J_SQUARED:
				Kratio = (1-log(Jeff))/(Jeff*Jeff);
				break;
				
			case HALF_J_SQUARED_MINUS_1_MINUS_LN_J:
			default:
				Kratio = 0.5*(Jeff + 1./Jeff);
				break;
		}
		
		// get reduced bulk modulus (used in AV and maybe hardening)
		p->Kred = C0squared*Kratio;
		
		// particle isentropic temperature increment
		dTq0 = -J*Kratio*gamma0*mptr->pPreviousTemperature*delV;
    }
	
    // artifical viscosity
	// delV is total incremental volumetric strain = total Delta(V)/V
    if(delV<0. && artificialViscosity)
	{	double QAVred = GetArtificialViscosity(delV/delTime,sqrt(p->Kred*J),mptr);
        AVEnergy += fabs(QAVred*delV);
		P += QAVred;
    }
    
    // set final pressure
    mptr->SetPressure(P);
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic terms, deviatoric later
    double avgP = 0.5*(P0+P);
	double delVres = 1. - 1./detdFres;
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-avgP*delVres);

	// SCGL and SL shear modulus and save Gratio = Jeff G/G0 for later calculations
    // Note: Jeff in Gred and Gratio is so that where they are used, they give
    //          specific Cauchy stress
    p->Gred = G1sp * plasticLaw->GetShearRatio(mptr,avgP,Jeff,p->hardProps,offset);
}

#pragma mark HEMGEOSMaterial::Accessors

// return unique, short name for this material
const char *HEMGEOSMaterial::MaterialType(void) const { return "Hyperelastic MGEOS Material"; }

// Calculate current wave speed in mm/sec. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
// but K/rho = Kred*J and G/rho = Gred*J (in mm^2/sec^2)
double HEMGEOSMaterial::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{
    // compressive volumetric strain x = 1-J
    double J = mptr->GetRelativeVolume();
    double x = 1. - J;
    
    // get K/rho0, but this ignores slope of energy term
    double KcurrRed,pressure;
    if(x>0.)
	{	// get bulk modulus
		if(x<=Xmax)
		{	KcurrRed = C0squared*GetMGEOSKRatio(x,gamma0,S1,S2,S3);
		}
		else
		{	KcurrRed = C0squared*Kmax;
			x = Xmax;
		}
		
		// compression law for pressure
		// denominator = 1 - S1*x - S2*x^2 - S3*x^3
		double denom = 1./(1. - x*(S1 + x*(S2 + x*S3)));
        
		// current effective and reduced (by rho0) bulk modulus
		double Klawred = C0squared*(1.-0.5*gamma0*x)*denom*denom;
		
		// get pressure
		double e = mptr->GetInternalEnergy();
		pressure = J*(Klawred*x + gamma0*e);
    }
    else
    {   // In tension hyperelastic law P = - K0(J-1)
		// does not account for Jeff and residual stresses
        KcurrRed = C0squared;
		pressure = -J*KcurrRed*(J-1.);
    }
    KcurrRed *= J;          // convert to K/rho
    
    // get G/rho at current pressure
    double GcurrRed = J*(G1sp * plasticLaw->GetShearRatio(mptr,pressure,J,NULL,offset));
    
    // return current save speed
    return sqrt((KcurrRed + 4.*GcurrRed/3.));
}

// not supported because moisture swelling may change MGEOS law (doe not work in Diffusion either)
bool HEMGEOSMaterial::SupportsDiffusion(void) const { return false; }

