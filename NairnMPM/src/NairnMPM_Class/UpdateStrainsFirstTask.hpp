/********************************************************************************
	UpdateStrainsFirstTask.hpp
	nairn-mpm-fea

	Created by John Nairn on July 22, 2010
	Copyright (c) 2010 John A. Nairn, All rights reserved.

	Dependencies
		MPMTask, CommonTask
********************************************************************************/

#ifndef _UPDATESTRAINSFIRSTTASK_

#define _UPDATESTRAINSFIRSTTASK_

#include "NairnMPM_Class/MPMTask.hpp"

class UpdateStrainsFirstTask : public MPMTask
{
	public:
        static void **matBuffer;
        static void **altBuffer;
	
		// constructor
		UpdateStrainsFirstTask(const char *);
	
		// required methods
		virtual bool Execute(int);
	
#ifdef RESTART_OPTION
        virtual bool BlockRestartOption(void) const;
#endif

        // class methods
        static void FullStrainUpdate(double,int,int,bool);
        static void CreatePropertyBuffers(int);
	
	protected:
};

extern UpdateStrainsFirstTask *USFTask;

#endif
