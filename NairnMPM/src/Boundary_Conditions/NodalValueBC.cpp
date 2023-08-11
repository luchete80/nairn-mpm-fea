/********************************************************************************
    NodalValueBC.cpp
    nairn-mpm-fea
 
    Parent class to setting scalar value on grid such as for
        temperature and concentration.

    Created by John Nairn on Sep 20, 2017.
    Copyright (c) 2017 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Boundary_Conditions/NodalValueBC.hpp"
#include "Nodes/NodalPoint.hpp"

#pragma mark NodalValueBC: Constructors and Destructors

NodalValueBC::NodalValueBC(int num,int setStyle,double concentration,double argTime)
                : BoundaryCondition(setStyle,concentration,argTime)
{
    nodeNum=num;
}

#pragma mark NodalValueBC: Methods

// save nodal concentration and zero it
NodalValueBC *NodalValueBC::CopyNodalValue(NodalPoint *nd,TransportField *gTrans)
{
    //cout << fmobj->mstep << ": copy node " << nd->num << " value " << gTrans->gTValue << endl;
    valueNoBC = gTrans->gTValue;
    return (NodalValueBC *)GetNextObject();
}

// restore nodal concentration and get initial force to cancel no-BC result
NodalValueBC *NodalValueBC::PasteNodalValue(NodalPoint *nd,TransportField *gTrans)
{
    //cout << fmobj->mstep << ": paste node " << nd->num << " value " << valueNoBC << endl;
    gTrans->gTValue = valueNoBC;
    return (NodalValueBC *)GetNextObject();
}

// initialize reaction flow at constant temperature boundary conditions
void NodalValueBC::InitQReaction(void) { qreaction = 0.; }

// add flow required to bring global nodal temperature to the BC temperature
void NodalValueBC::SuperposeQReaction(double qflow) { qreaction += qflow; }

