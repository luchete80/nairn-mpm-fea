/********************************************************************************
    MpsController.cpp
    NairnFEA
    
    Created by John Nairn on 6/27/05.
    Copyright (c) 2005 John A. Nairn, All rights reserved.
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Read_MPM/MpsController.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Global_Quantities/ThermalRamp.hpp"
#include "NairnMPM_Class/ResetElementsTask.hpp"

#include "System/MPMPrefix.hpp"

MpsController *mpCtrl=NULL;

/********************************************************************************
	MpsController: methods
********************************************************************************/

// add new material point
void MpsController::AddMaterialPoint(MPMBase *newMpt,double conc,double temp)
{
	AddObject(newMpt);
	newMpt->SetConcentration(conc,false);
	newMpt->SetTemperature(temp,thermal.reference);
}

// set position (if "pt") or velocity (if "vel")
int MpsController::SetPtOrVel(char *xName,Vector *value)
{
	MPMBase *newMpt=(MPMBase *)lastObject;
	if(newMpt==NULL) return FALSE;		// error if no particle to set
	if(strcmp(xName,"pt")==0)
	{	newMpt->SetPosition(value);
		newMpt->SetOrigin(value);
		if(!ResetElementsTask::ResetElement(newMpt)) return FALSE;
		newMpt->SetElementCrossings(0);
		
	}
	else
		newMpt->SetVelocity(value);
	return TRUE;
}

// set particle mass
void MpsController::SetPtMass(double pmass)
{
	MPMBase *newMpt=(MPMBase *)lastObject;
	if(newMpt==NULL) return;	// ignore if no particle to set
	newMpt->mp=pmass;
}

// assemble into array used in the code
int MpsController::SetMPArray(void)
{
	// make 0-based array of material points
	if(numObjects==0) return false;
	mpm = new (std::nothrow) MPMBase *[numObjects];
	if(mpm==NULL) return false;
	
	// fill the array
	MPMBase *obj=(MPMBase *)firstObject;
	nmpms=0;
	while(obj!=NULL)
	{	mpm[nmpms]=obj;
		nmpms++;
		obj=(MPMBase *)obj->GetNextObject();
	}
	return true;
}
