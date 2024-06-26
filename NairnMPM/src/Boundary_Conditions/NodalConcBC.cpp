/********************************************************************************
    NodalConcBC.cpp
    nairn-mpm-fea
    
    Created by John Nairn on Thu Apr 1 2004.
    Copyright (c) 2004 John A. Nairn, All rights reserved.
********************************************************************************/
#if defined ( _MSC_VER) || defined (__APPLE__) 
#include "stdafx.h"
#endif
#include "Boundary_Conditions/NodalConcBC.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "System/UnitsController.hpp"


// Nodal concentration BC globals
NodalConcBC *firstConcBC=NULL;
NodalConcBC *lastConcBC=NULL;
NodalConcBC *firstRigidConcBC=NULL;
NodalConcBC *reuseRigidConcBC=NULL;
NodalConcBC *firstDiffBC[NUM_DUFFUSION_OPTIONS];

#pragma mark NodalConcBC: Constructors and Destructors

NodalConcBC::NodalConcBC(int num,int setStyle,double concentration,double argTime,int phase)
                    : NodalValueBC(num,setStyle,concentration,argTime)
{
    phaseStyle = phase;
    if(phaseStyle<FRACTURE_PHASE_FIELD)
    {   phaseStyle = MOISTURE_DIFFUSION;
#ifdef POROELASTICITY
        if(DiffusionTask::HasPoroelasticity())
            phaseStyle = POROELASTICITY_DIFFUSION;
#endif
    }
    
    // set fixed direction
    switch(phaseStyle)
    {   case MOISTURE_DIFFUSION:
        case POROELASTICITY_DIFFUSION:
            nd[num]->SetFixedDirection(CONC_DIRECTION);
        case FRACTURE_PHASE_FIELD:
            nd[num]->SetFixedDirection(FRACTURE_PHASE_DIRECTION);
        case BATTERY_PHASE_FIELD:
            nd[num]->SetFixedDirection(BATTERY_PHASE_DIRECTION);
        case CONDUCTION_PHASE_FIELD:
            nd[num]->SetFixedDirection(POISSON_DIFF_DIRECTION);
        default:
            nd[num]->SetFixedDirection(CONC_DIRECTION);
    }
}

// print it (if can be used)
BoundaryCondition *NodalConcBC::PrintBC(ostream &os)
{
    double scale = 1.;
#ifdef POROELASTICITY
    if(phaseStyle==POROELASTICITY_DIFFUSION) scale = UnitsController::Scaling(1.e-6);
#endif
    char nline[200];
	size_t nlsize=200;
    snprintf(nline,nlsize,"%7d %2d %15.7e %15.7e %3d",nodeNum,style,scale*GetBCValueOut(),GetBCFirstTimeOut(),phaseStyle);
    os << nline;
    PrintFunction(os);
	return (BoundaryCondition *)GetNextObject();
}

// when getting total diffusion flux, add qreaction to input double
// if matchID==0 include it, otherwise include only if matchID equals bcID
NodalConcBC *NodalConcBC::AddConcReaction(double *totalReaction,int matchID)
{	if(bcID==matchID || matchID==0)
		*totalReaction += qreaction;
	return (NodalConcBC *)GetNextObject();
}

#pragma mark NodalConcBC: ACCESSORS

// set value (and scale legacy MPa to Pa)
void NodalConcBC::SetBCValue(double bcvalue)
{	double rescale = DiffusionTask::RescalePotential(phaseStyle);
	BoundaryCondition::SetBCValue(rescale*bcvalue);
}

// get set direction (changes in constructor above too)
int NodalConcBC::GetSetDirection(void) const
{   // depends on style
    switch(phaseStyle)
    {   case MOISTURE_DIFFUSION:
        case POROELASTICITY_DIFFUSION:
            return CONC_DIRECTION;
        case FRACTURE_PHASE_FIELD:
            return FRACTURE_PHASE_DIRECTION;
        case BATTERY_PHASE_FIELD:
            return BATTERY_PHASE_DIRECTION;
        case CONDUCTION_PHASE_FIELD:
            return POISSON_DIFF_DIRECTION;
        default:
            return CONC_DIRECTION;
    }
}

/**********************************************************
	Sum all reaction diffusion forces for all grid BCs
	If ID is not zero, only includes those with matching ID
	Result is in transport content/sec
*/
double NodalConcBC::TotalConcReaction(int matchID,DiffusionTask *qdiff)
{
	double reactionTotal = 0.;
	NodalConcBC *nextBC = (NodalConcBC *)qdiff->GetFirstBCPtr();
	
	// Add each value
	while(nextBC!=NULL)
		nextBC=nextBC->AddConcReaction(&reactionTotal,matchID);
	
	return reactionTotal;
}
