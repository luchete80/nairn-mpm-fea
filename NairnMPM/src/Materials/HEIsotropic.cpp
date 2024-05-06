/********************************************************************************
    HEIsotropic.cpp
    nairn-mpm-fea
    
    Created by John Nairn, Sept 27, 2011.
    Copyright (c) 2011 John A. Nairn, All rights reserved.

	Dependencies
		HyperElastic.hpp (MaterialBase.hpp)
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "HEIsotropic.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Materials/LinearHardening.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "System/UnitsController.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"

#pragma mark HEIsotropic::Constructors and Destructors

// Constructor
HEIsotropic::HEIsotropic(char *matName,int matID) : HyperElastic(matName,matID)
{
   	G1 = -1.;			// required
    
	// JAN: hard-code linear hardening - future can set to other laws
	plasticLaw = new LinearHardening(this);
}


#pragma mark HEIsotropic::Initialization

// Read material properties
char *HEIsotropic::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    input=DOUBLE_NUM;
    
    if(strcmp(xName,"G1")==0 || strcmp(xName,"G")==0)
        return UnitsController::ScaledPtr((char *)&G1,gScaling,1.e6);
    
    // look for different plastic law
    if(strcmp(xName,"Hardening")==0)
    {	input = HARDENING_LAW_SELECTION;
        return (char *)this;
    }
    
    // JAN: Move yielding properties to the hardening law (yield, Ep, Khard)
	// check plastic law
    char *ptr = plasticLaw->InputHardeningProperty(xName,input,gScaling);
    if(ptr != NULL) return ptr;
    
    return(HyperElastic::InputMaterialProperty(xName,input,gScaling));
}

// Allows any hardening law
bool HEIsotropic::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID)
{   delete plasticLaw;
    plasticLaw = pLaw;
    return true;
}

// return plastic law ID (or 0 if none)
HardeningLawBase *HEIsotropic::GetPlasticLaw(void) const { return plasticLaw; }

// verify settings and maybe some initial calculations
const char *HEIsotropic::VerifyAndLoadProperties(int np)
{
	// call plastic law first
    const char *ptr = plasticLaw->VerifyAndLoadProperties(np);
    if(ptr != NULL) return ptr;
    
    // user must enter G1 and Kbulk
    if(G1<0. || Kbulk < 0. )
		return "HEIsotropic Material needs non-negative G1 and K";
    
	// G in specific units using initial rho (F/L^2 L^3/mass)
 	G1sp = G1/rho;
	
    // heating gamma0 (dimensionless)
    double alphaV = 3.e-6*aI;
    gammaI = Kbulk*alphaV/(rho*heatCapacity);
	
	// must call super class
	return HyperElastic::VerifyAndLoadProperties(np);
}

// plane stress not allowed in viscoelasticity
// throws CommonException()
void HEIsotropic::ValidateForUse(int np) const
{	if(np==PLANE_STRESS_MPM)
    {	throw CommonException("HEIsotropic materials cannot be used in plane stress MPM yet",
                                "HEIsotropic::ValidateForUse");
    }
	
	//call super class (why can't call super class?)
	HyperElastic::ValidateForUse(np);
}

// print mechanical properties to the results
void HEIsotropic::PrintMechanicalProperties(void) const
{	
    PrintProperty("G1",G1*UnitsController::Scaling(1.e-6),"");
    PrintProperty("K",Kbulk*UnitsController::Scaling(1.e-6),"");
    cout << endl;
	double calcE = 9.*Kbulk*G1/(3.*Kbulk+G1);
    PrintProperty("E",calcE*UnitsController::Scaling(1.e-6),"");
    PrintProperty("nu",(3.*Kbulk-2.*G1)/(6.*Kbulk+2.*G1),"");
    cout << endl;
    
    PrintProperty("a",aI,"");
	PrintProperty("gam0",gammaI,"");
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
    
	// JAN: print hardening law properties
	plasticLaw->PrintYieldProperties();
	
	// call superclass here if it is not Material base
	HyperElastic::PrintMechanicalProperties();
}

#pragma mark HEIsotropic:History Data Methods

// return number of bytes needed for history data
int HEIsotropic::SizeOfHistoryData(void) const
{	return (plasticLaw->HistoryDoublesNeeded()+2)*sizeof(double); }

// First ones for hardening law. Particle J appended at the end
char *HEIsotropic::InitHistoryData(char *pchr,MPMBase *mptr)
{	J_History = plasticLaw->HistoryDoublesNeeded();
	double *p = CreateAndZeroDoubles(pchr,J_History+2);
    plasticLaw->InitPlasticHistoryData(p);
	p[J_History]=1.;					// J
	p[J_History+1]=1.;					// Jres
	return (char *)p;
}

// reset history data
void HEIsotropic::ResetHistoryData(char *pchr,MPMBase *mptr)
{	double *p = (double *)pchr;
	plasticLaw->InitPlasticHistoryData(p);
	p[J_History]=1.;					// J
	p[J_History+1]=1.;					// Jres
}

// Number of history variables - plastic plus 2
int HEIsotropic::NumberOfHistoryDoubles(void) const { return J_History+2; }

#pragma mark HEIsotropic:Methods

// buffer size for mechanical properties
int HEIsotropic::SizeOfMechanicalProperties(int &altBufferSize) const
{   altBufferSize = plasticLaw->SizeOfHardeningProps();
    return sizeof(HEPlasticProperties);
}

// Get elastic and plastic properties, return null on error
void *HEIsotropic::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	HEPlasticProperties *p = (HEPlasticProperties *)matBuffer;
 	p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer,offset);
	double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),1.,p->hardProps,offset);
	
	// Find current p->Gred and p->Kred (reduced moduli)
	p->Gred = G1sp*Gratio;
	p->Kred = Ksp;
    
    return p;
}

/* Take increments in strain and calculate new Particle: strains, rotation strain,
        stresses, strain energy,
    dvij are (gradient rates X time increment) to give deformation gradient change
    For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
    This material tracks pressure and stores deviatoric stress only in particle stress
        tensor
*/
void HEIsotropic::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,
									 ResidualStrains *res,int historyOffset) const
{
	HEPlasticProperties *p = (HEPlasticProperties *)properties;

    // store initial stress
    Tensor *sp = mptr->GetStressTensor();
    Tensor st0 = *sp;
    
    // Compute Elastic Predictor
    // ============================================
    
	// Update total deformation gradient, and calculate trial B
    Tensor Btrial;
	double detdF = IncrementDeformation(mptr,du,&Btrial,np);
    
	// Deformation gradients and Cauchy tensor differ in plane stress and plane strain
    // This code handles plane strain, axisymmetric, and 3D - Plane stress is blocked
	double J = detdF * mptr->GetHistoryDble(J_History,historyOffset);
	mptr->SetHistoryDble(J_History,J,historyOffset);						// Stocking J
    
    // J is determinant of F (or sqrt root of determinant of B), Jeff is normalized to residual stretch
    double Jres = mptr->GetHistoryDble(J_History+1,historyOffset);
    double dJres = GetIncrementalResJ(mptr,res,Jres);
    Jres *= dJres;
	mptr->SetHistoryDble(J_History+1,Jres,historyOffset);
    double Jeff = J/Jres;
	
	// Get hydrostatic stress component in subroutine
	double dTq0 = 0.,dispEnergy = 0.;
    UpdatePressure(mptr,J,detdF,np,Jeff,delTime,p,res,dJres,historyOffset,dTq0,dispEnergy);
    
    // Others constants
    double J23 = pow(J, 2./3.)/Jres;
    
    // find Trial (Cauchy stress)/rho0
	// (Trial_s/rho0 = Trial_s*rho/(rho*rho0) = (Trial_tau*rho/rho0^2) = (1/J)*(Trial_tau/rho0)
    Tensor stk = GetTrialDevStressTensor(&Btrial,J*J23,np,p->Gred);
    
    // Checking for plastic loading
    // ============================================
    
    // Get magnitude of the deviatoric stress tensor
    // ||s|| = sqrt(s.s)
    
	// Set alpint for particle
	HardeningAlpha alpha;
	plasticLaw->UpdateTrialAlpha(mptr,np,&alpha,historyOffset);
	
	// Trial stress state
    double magnitude_strial = GetMagnitudeS(&stk,np);
    double gyld = plasticLaw->GetYield(mptr,np,delTime,&alpha,p->hardProps);
    double ftrial = magnitude_strial-SQRT_TWOTHIRDS*gyld;
    //cout << "  #magnitude_strial =   "<< magnitude_strial<< "  GetYield =   "<< gyld<< "  ftrial =   "<< ftrial<< endl;
    //cout << "  #yldred =   "<< yldred << "  Epred =   "<< Epred << "  gyld =   "<< gyld <<"  alpint =   "<< alpint<< "  ftrial =   "<< ftrial<< endl;
    
    // these will be needed for elastic or plastic
    Tensor *pB = mptr->GetAltStrainTensor();
    
    //============================
    //  TEST
    //============================
    if(ftrial<=0.)
	{	// if elastic
        //============================
		
		// save on particle
        *pB = Btrial;
        
        // Get specifique stress i.e. (Cauchy Stress)/rho = J*(Cauchy Stress)/rho0 = (Kirchoff Stress)/rho0
		// The deviatoric stress was calculated as (Cauchy Stress)/rho, so need to scale by J to get correct stress
        sp->xx = J*stk.xx;
		sp->yy = J*stk.yy;
        sp->xy = J*stk.xy;
        sp->zz = J*stk.zz;
       
        // work energy per unit mass (U/(rho0 V0)) and we are using
        // W(F) as the energy density per reference volume V0 (U/V0) and not current volume V
		double workEnergy = 0.5*((st0.xx+sp->xx)*du(0,0)
								  + (st0.yy+sp->yy)*du(1,1)
								  + (st0.zz+sp->zz)*du(2,2)
								  + (st0.xy+sp->xy)*(du(1,0)+du(0,1)));
        if(np==THREED_MPM)
		{	sp->xz = J*stk.xz;
			sp->yz = J*stk.yz;
			workEnergy += 0.5*((st0.yz+sp->yz)*(du(2,1)+du(1,2))
                                       + (st0.xz+sp->xz)*(du(2,0)+du(0,2)));
        }
        mptr->AddWorkEnergy(workEnergy);

		// residual energy or sigma.deres - it is zero here for isotropic material
		// because deviatoric stress is traceless and deres has zero shear terms
		// residual energy due to pressure was added in the pressure update
		
		// heat energy
		IncrementHeatEnergy(mptr,dTq0,dispEnergy);
        
		// give material chance to update history variables that change in elastic updates
		plasticLaw->ElasticUpdateFinished(mptr,np,delTime,historyOffset);
		
        return;
    }
    
    // Plastic behavior - Return Mapping algorithm 
    //=====================================================
    // if plastic
    
	// JAN: Use hardening law method (which can now use other laws too)
    double Ie1bar = (Btrial.xx+Btrial.yy+Btrial.zz)/(3.*J23);
    double MUbar = Jres*p->Gred*Ie1bar;
    
    // Find  lambda for this plastic state
	double dlambda = plasticLaw->SolveForLambdaBracketed(mptr,np,magnitude_strial,&stk,
														 MUbar,1.,1.,delTime,&alpha,p->hardProps,historyOffset);
    
    // update deviatoric stress (need to scale by J to get to Kirchoff stress/rho
	Tensor nk = GetNormalTensor(&stk,magnitude_strial,np);
    //cout << "nk.xx  =    " << nk.xx << "nk.xy  =    " << nk.xy << endl;
	double twoMuLam = 2.*MUbar*dlambda;
    sp->xx = J*(stk.xx - twoMuLam*nk.xx);
    sp->yy = J*(stk.yy - twoMuLam*nk.yy);
    sp->zz = J*(stk.zz - twoMuLam*nk.zz);
    sp->xy = J*(stk.xy - twoMuLam*nk.xy);
    if(np == THREED_MPM)
    {   sp->xz = J*(stk.xz - twoMuLam*nk.xz);
        sp->yz = J*(stk.yz - twoMuLam*nk.yz);
    }
    
    // save on particle
	double twoThirdsLamI1bar = 2.*dlambda*Ie1bar;
    pB->xx = Btrial.xx - twoThirdsLamI1bar*nk.xx;
    pB->yy = Btrial.yy - twoThirdsLamI1bar*nk.yy;
    pB->zz = Btrial.zz - twoThirdsLamI1bar*nk.zz;
	pB->xy = Btrial.xy - twoThirdsLamI1bar*nk.xy;
    if(np == THREED_MPM)
    {   pB->xz = Btrial.xz - twoThirdsLamI1bar*nk.xz;
        pB->yz = Btrial.yz - twoThirdsLamI1bar*nk.yz;
    }
    /* Old method collecting B from stresses
	pB->xx = (sp->xx/p->Gred+Ie1bar)*J23;
	pB->yy = (sp->yy/p->Gred+Ie1bar)*J23;
	pB->zz = (sp->zz/p->Gred+Ie1bar)*J23;
	pB->xy = sp->xy*J23/p->Gred;
    if(np == THREED_MPM)
    {   pB->xz = sp->xz*J23/p->Gred;
        pB->yz = sp->yz*J23/p->Gred;
    }
    */
    
    // strain energy per unit mass (U/(rho0 V0)) and we are using
    double workEnergy = 0.5*((st0.xx+sp->xx)*du(0,0)
                               + (st0.yy+sp->yy)*du(1,1)
                               + (st0.zz+sp->zz)*du(2,2)
                               + (st0.xy+sp->xy)*(du(1,0)+du(0,1)));
    if(np==THREED_MPM)
    {   workEnergy += 0.5*((st0.yz+sp->yz)*(du(2,1)+du(1,2))
                                   + (st0.xz+sp->xz)*(du(2,0)+du(0,2)));
    }
    
    // total work
    mptr->AddWorkEnergy(workEnergy);

    // residual energy or sigma.deres - it is zero here for isotropic material
    // because deviatoric stress is traceless and deres has zero shear terms
    // residual energy due to pressure was added in the pressure update
    
    // Plastic work increment per unit mass (dw/(rho0 V0)) (nJ/g)
    dispEnergy += dlambda*(sp->xx*nk.xx + sp->yy*nk.yy + sp->zz*nk.zz + 2.*sp->xy*nk.xy);
    if(np==THREED_MPM)  dispEnergy += 2.*dlambda*(sp->xz*nk.xz + sp->yz*nk.yz);
	
    // Subtract q.dalpha to get final disispated energy per unit mass (dPhi/(rho0 V0)) (nJ/g)
    dispEnergy -= dlambda*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
    
	// The cumulative dissipated energy is tracked in plastic energy
    mptr->AddPlastEnergy(dispEnergy);
    
    // heat energy
    IncrementHeatEnergy(mptr,dTq0,dispEnergy);
    
	// update internal variables in the plastic law
	plasticLaw->UpdatePlasticInternal(mptr,np,&alpha,historyOffset);
}

#pragma mark HEIsotropic::Custom Methods

// To allow better subclassing it is better to separate out calculations
//  of dilation energy. This version updates pressure (i.e. dilational
//  contribution to normal stress) and adds incremental energy to strain energy
// J = V(T,c)/V0(T,c), Jeff = V(T,c)/V0(Tref,cref), Jres = V0(T,c)/V0(Tref,cref) = J/Jeff
// Jn+1 = (detdF/detdFres) Jn, Jresn+1 = detdFres Jresn, Jtot = detdF Jtotn
// detdFres = (1+dres)^3 (approximately)
// Here Tref and cref are starting conditions and T and c are current temperature and moisture
void HEIsotropic::UpdatePressure(MPMBase *mptr,double J,double detdF,int np,double Jeff,
								 double delTime,HEPlasticProperties *p,ResidualStrains *res,double detdFres,int offset,
								 double &dTq0,double &AVEnergy) const
{
	double Kterm = J*GetVolumetricTerms(Jeff,p->Kred);       // times J to get Kirchoff stress
    double P0 = mptr->GetPressure();
    
    // artifical viscosity
	// delV is total incremental volumetric strain = total Delta(V)/V, here it is (Vn+1-Vn)/Vn+1
	double delV = 1. - 1./detdF;
    double QAVred = 0.;
    if(delV<0. && artificialViscosity)
	{	QAVred = GetArtificialViscosity(delV/delTime,sqrt(p->Kred*J),mptr);
        AVEnergy += fabs(QAVred*delV);
    }
    double Pfinal = -Kterm + QAVred;
    mptr->SetPressure(Pfinal);
    
    // work energy is dU = -P dV + s.de(total)
	// Here do hydrostatic term
    // Internal energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
	// Also get residual energy from -P (Vresn+1-Vresn)/Vresn+1 = -P delVres
    double avgP = 0.5*(P0+Pfinal);
	double delVres = 1. - 1./detdFres;
    mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-avgP*delVres);
	
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
	dTq0 = -J*Kratio*gammaI*mptr->pPreviousTemperature*delV;
}

// get trial deviatoric stress tensor based on trial B
// Note that input J23 is actually J^(2/3)/Jres
Tensor HEIsotropic::GetTrialDevStressTensor(Tensor *B,double J23,int np,double Gred) const
{
    // Trial Kirchhoff Stress Tensor / rho0 = J Cauchy/rho0 = Cauchy/rho
    Tensor strial;
	double J23normal = 3.*J23;
    strial.xx = Gred*(2.*B->xx-B->yy-B->zz)/J23normal;
    strial.yy = Gred*(2.*B->yy-B->xx-B->zz)/J23normal;
    //strial.zz = Gred*(2.*B->zz-B->xx-B->yy)/J23normal;
	strial.zz = -(strial.xx+strial.yy);			// slightly faster and making use of traceless
    strial.xy = Gred*B->xy/J23;
    if(np==THREED_MPM)
    {   strial.xz = Gred*B->xz/J23;
        strial.yz = Gred*B->yz/J23;
    }
	else
	{	strial.xz = 0.;
		strial.yz = 0.;
	}
    
    return strial;
}

// Get magnitude of s = sqrt(s.s) when s is a deviatoric stress
double HEIsotropic::GetMagnitudeS(Tensor *st,int np) const
{
	double s,t;
	
	switch(np)
    {   case THREED_MPM:
            s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
            t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
            break;
            
		default:
			s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
			t = st->xy*st->xy;
			break;
	}
	return sqrt(s+t+t);
}

// Implementation of the normal to the yield surface
Tensor HEIsotropic::GetNormalTensor(Tensor *strial,double magnitude_strial,int np) const
{
    // Trial specific Kirchhoff Stress Tensor
    Tensor n;
    n.xx = strial->xx/magnitude_strial;
    n.yy = strial->yy/magnitude_strial;
    n.zz = strial->zz/magnitude_strial;
    n.xy = strial->xy/magnitude_strial;
    if(np==THREED_MPM)
    {   n.xz = strial->xz/magnitude_strial;
        n.yz = strial->yz/magnitude_strial;
    }
	else
	{	n.xz = 0.;
		n.yz = 0.;
	}
    //cout << "strial.yy  =    " << strial->yy << "      strial.zz  =    " << strial->zz << "magnitude_strial  =    " << magnitude_strial << endl;
    //cout << "n.xx  =    " << n.xx << "n.xy  =    " << n.xy << endl;
    return n;
}

#pragma mark HEIsotropic::Accessors

// convert J to K using isotropic method
Vector HEIsotropic::ConvertJToK(Vector d,Vector C,Vector J0,int np)
{	double nuLS = (3.*Kbulk-2.*G1)/(6.*Kbulk+2.*G1);
	return IsotropicJToK(d,C,J0,np,nuLS,G1);
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor HEIsotropic::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{	return GetStressPandDev(sp,pressure,mptr);
}

// store a new total stress on a particle's stress and pressure variables
void HEIsotropic::SetStress(Tensor *spnew,MPMBase *mptr) const
{	SetStressPandDev(spnew,mptr);
}

// Increment thickness (zz) stress through deviatoric stress and pressure
void HEIsotropic::IncrementThicknessStress(double dszz,MPMBase *mptr) const
{	IncrementThicknessStressPandDev(dszz,mptr);
}

// return unique, short name for this material
const char *HEIsotropic::MaterialType(void) const { return "Hyperelastic Isotropic"; }

// calculate wave speed in mm/sec (props in mass/(L sec^2) and rho in mass/L^3)
double HEIsotropic::WaveSpeed(bool threeD,MPMBase *mptr) const
{	return sqrt((Kbulk+4.*G1/3.)/rho);
}

// Calculate current wave speed. Uses sqrt((K+4G/3)/rho) which is dilational wave speed
// Adjusts K, but not sure how to change G
double HEIsotropic::CurrentWaveSpeed(bool threeD,MPMBase *mptr,int offset) const
{   double Jeff = mptr->GetHistoryDble(J_History,offset)/mptr->GetHistoryDble(J_History+1,offset);
    double Kratio;
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
    return sqrt((Kratio*Kbulk+4.*G1/3.)/rho);
}

// if a subclass material supports artificial viscosity, override this and return TRUE
bool HEIsotropic::SupportsArtificialViscosity(void) const { return true; }

// store elastic B in alt strain while F has elastic + plastic deformation
int HEIsotropic::AltStrainContains(void) const { return LEFT_CAUCHY_ELASTIC_B_STRAIN; }
