/****************************************************************************************************************************
    MaterialBase.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Wed Jan 23 2002.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
	
	Material Classes
		+ indicates actual material (others cannot be materials)
		(MPM) or (FEA) indicate only that code, otherwise in both
	
	MaterialBase (base class)
****************************************************************************************************************************/

#include "stdafx.h"
#include "Materials/MaterialBase.hpp"
#include "System/UnitsController.hpp"

// Constants for output fields for material properties
#define PROP_LABEL_LENGTH 4
#define PROP_COLUMN 20

// Materials globals
MaterialBase **theMaterials=NULL;		// list of materials
int nmat=0;								// number of materials

#pragma mark MaterialBase::Constructors and Destructors

// Constructors
MaterialBase::MaterialBase()
{	// unknow material and nothing else gets set
	name=new char[20];
	strcpy(name,"UNKNOWN MATERIAL");
	materialID = -1;
}

// Constructors
// throws std::bad_alloc
MaterialBase::MaterialBase(char *matName,int matID)
{
	name=new char[strlen(matName)+1];
    strcpy(name,matName);
	materialID = matID;
	concSaturation=1.;
	betaI=0.;
	red=-1.;
#ifdef MPM_CODE
    rho=1.;
    KIc=KIIc=JIc=JIIc=-1.;				// traction laws assumes -1
	delIc=delIIc=-1.;					// traction laws assumes -1
	nmix=1.;							// mixed mode exponent (if used)
	maxLength=-1.;
	initTime=-1.;
	initSpeed=1.;
	diffusionCon=0.;
	// For diffusion. Materials that support poroelasticy change it in VerifyAndLoadProps()
	diffusionCT=1.;
    // KIexp and KIIexp are the parameters of the empirical fracture criterion:
    // (KI/KIc)^KIexp+(KII/KIIc)^KIIexp=1.
    // Set the default values to 2 for elliptical fracture locus.
    KIexp=KIIexp=2.;
    criterion[0]=criterion[1]=UNSPECIFIED;
	constantDirection = false;
	constantTip = false;
	growDir.x=0.;
	growDir.y=0.;
	matPropagateDirection[0]=matPropagateDirection[1]=UNSPECIFIED;
	tractionMat[0]=tractionMat[1]=0;
	field=-1;							// material velocity field
	shareMatField=0;					// share field with another material
	activeField=-1;
	kCond=0.;
	heatCapacity=UnitsController::Scaling(1.e6);		// keep one because needed by ideal gas
	lastFriction=NULL;
    artificialViscosity = false;
    avA1 = 0.2;
    avA2 = 2.0;
	matPdamping=0;
	matUsePDamping=false;
    allowsCracks= true;

#ifdef POROELASTICITY
	Darcy=0.;
	alphaPE=0.;
	Ku=UnitsController::Scaling(1.e6);			// Legacy MPa and default 1 MPa
#endif
#endif
}

// Destructor (and it is virtual)
MaterialBase::~MaterialBase()
{	delete [] name;
}

#pragma mark MaterialBase::Initialization

// print to output window
void MaterialBase::PrintMaterial(int num) const
{
	// material name
    cout << "Material " << num << ": " << name << endl;
    cout << "     " << MaterialType() << " with:" << endl;
	
	// call class to list mechanical properties
	PrintMechanicalProperties();
	
	// properties in common base class (differs for MPM and FEA)
	PrintCommonProperties();
    
#ifdef MPM_CODE
	// Transport Properties
	PrintTransportProperties();
#endif
}

// print mechanical properties to output window
void MaterialBase::PrintMechanicalProperties(void) const {}

// print property with units (option) and align to columns for material properties
// In field of width PROP_COLUMN
void MaterialBase::PrintProperty(const char *propName,double value,const char *unitsStr)
{
	char prop[200];
	
	// name
	strcpy(prop,propName);
	int length=(int)strlen(prop);
	while(length<PROP_LABEL_LENGTH)
	{	strcat(prop," ");
		length++;
	}
	
	// value
	char valueStr[50];
	sprintf(valueStr,"= %g",value);
	strcat(prop,valueStr);
	
	// units
	if(strlen(unitsStr)>0)
	{	strcat(prop," ");
		strcat(prop,unitsStr);
	}
	
	// pad to column
	length=(int)strlen(prop);
	int padLength = (length<PROP_COLUMN) ? PROP_COLUMN : 2*PROP_COLUMN ;
	while(length<padLength)
	{	strcat(prop," ");
		length++;
	}
	
	cout << prop;
}

// print text aligned with material columns, left or right justified
// In field of width PROP_COLUMN, will use two if needed
void MaterialBase::PrintProperty(const char *text,bool rightJustify)
{
	char prop[200];
	
	// up to two columns wide
	int length=(int)strlen(text);
	int padLength = (length<PROP_COLUMN) ? PROP_COLUMN : 2*PROP_COLUMN ;
	strcpy(prop,"");
	
	// right justify
	if(rightJustify)
	{	int pad=padLength-length-1;
		while(pad>0)
		{	strcat(prop," ");
			pad--;
		}
		strcat(prop,text);
		strcat(prop," ");
	}
	
	// left justify
	else
	{	strcat(prop,text);
		while(length<padLength)
		{	strcat(prop," ");
			length++;
		}
	}
	
	cout << prop;
}

    
// return material type and ID
const char *MaterialBase::MaterialType(void) const { return "Unknown Material Type"; }
const int MaterialBase::MaterialID(void) const { return materialID; }

