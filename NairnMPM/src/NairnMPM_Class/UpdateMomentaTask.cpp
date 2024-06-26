/********************************************************************************
	UpdateMomentaTask.cpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.
 
	The tasks are:
	-------------
	* Updates momenta on the nodes using
		pk(i+1) = pk(i) + ftot(i)*dt
	  for each velocity field on each node.
	* If transport activated, find transport rate by dividing
	  transport flow by transport mass
	* Once get new momenta, check for material contact. Crack contact
	  is checked in a separate step outside the main loop. The material
	  contact checks all nodes. The crack contact	looks only at nodes known
	  to have cracks
	  Note: If either contact changes momenta, change force too to keep consistent with
	  momentum change (because not in post-update tasks)
********************************************************************************/

#if defined ( _MSC_VER) || defined (__APPLE__)
#include "stdafx.h"
#endif
#include "NairnMPM_Class/UpdateMomentaTask.hpp"
#include "NairnMPM_Class/NairnMPM.hpp"
#include "Custom_Tasks/TransportTask.hpp"
#include "Nodes/NodalPoint.hpp"
#include "Cracks/CrackNode.hpp"
#include "Exceptions/CommonException.hpp"
#include "Nodes/MaterialContactNode.hpp"
#include "Boundary_Conditions/NodalVelBC.hpp"
#include "Global_Quantities/BodyForce.hpp"

#pragma mark CONSTRUCTORS

UpdateMomentaTask::UpdateMomentaTask(const char *name) : MPMTask(name)
{
}

#pragma mark REQUIRED METHODS

// Update grid momenta and transport rates
// throws CommonException()
bool UpdateMomentaTask::Execute(int taskOption)
{	
#pragma omp parallel for
	for(int i=1;i<=*nda;i++)
	{   NodalPoint *ndptr = nd[nda[i]];
        
		// update nodal momenta
		ndptr->UpdateMomentum(timestep);
        
        // get grid transport rates
		TransportTask::UpdateTransportOnGrid(ndptr);
	}

	// contact and BCs
	ContactAndMomentaBCs(UPDATE_MOMENTUM_CALL);
	
	// impose transport BCs on the grid
	TransportTask::TransportGridBCs(mtime,timestep,UPDATE_MOMENTUM_CALL);
    
    return true;
}

// Do contact calculations and impose momenta conditions
// passType == MASS_MOMENTUM_CALL, UPDATE_MOMENTUM_CALL, UPDATE_STRAINS_LAST_CALL
// For MASS_MOMENTUM_CALL and UPDATE_STRAINS_LAST_CALL, always impose momenta conditions, but:
//		MASS_MOMENTUM_CALL only does symetry BCs if using USL- or USL+ because no strain being found
//		UPDATE_STRAINS_LAST_CALL only occurs for USL+ and USAVG+
//		These only change p that is needed for subsequent strain update
// For UPDATE_MOMENTUM_CALL, update gets correct momenta, so only done if contact
//		might have change them (and change p and f to keep consistent)
void UpdateMomentaTask::ContactAndMomentaBCs(int passType)
{
	// material contact
	bool hasContact = MaterialContactNode::ContactOnKnownNodes(timestep,passType);
	
	// adjust momenta and forces for crack contact on known nodes
	hasContact = CrackNode::ContactOnKnownNodes(timestep,passType) ? true : hasContact;
    
    // In theory, if the above did not change momenta (i.e., hasCracks=false) and this
    // if the UPDATE_MOMENTUM_CALL, the velocity should be correct (the update should
    // have changed initial velocities into final velocites. In other words, the
    // code can exit. A method equivalent to follow call was add in revision 3158
    // of OSPARTICULAS. It was noticed in revision 3618 that exiting now changed
    // the results. I am not sure which method is best, but it is likly there is
    // no problem making sure velocities are exact as set. The following line
    // was therefore commented out in revision 3619
    //if(UPDATE_MOMENTUM_CALL && !hasContact) return;

	// For FLIP, reimpose velocity BCs. For FMPM impose them unless using option
    // that only does velocity BCs in the particle update
	if(bodyFrc.GetXPICOrder()<=1 || bodyFrc.GridBCOption()!=GRIDBC_VELOCITY_ONLY)
    {	NodalVelBC::GridVelocityConditions(passType);
    }
}
