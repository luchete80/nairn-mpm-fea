/********************************************************************************
    NodalValueBC.hpp
    nairn-mpm-fea
 
    Created by John Nairn on Sep 20, 2017.
    Copyright (c) 2017 John A. Nairn, All rights reserved.
	
    Dependencies
        BoundaryCondition.hpp
 ********************************************************************************/

#ifndef _NODALVALUEBC_

#define _NODALVALUEBC_

#include "Boundary_Conditions/BoundaryCondition.hpp"
#include "System/DataTypes.hpp"

class NodalPoint;

class NodalValueBC: public BoundaryCondition
{
    public:
        double valueNoBC;
        double nodalActivity;
    
        // constructors and destructors
        NodalValueBC(int,int,double,double);
    
        // methods
        virtual NodalValueBC *CopyNodalValue(NodalPoint *,TransportField *);
        virtual NodalValueBC *CopyNodalValue(NodalPoint *,TransportField *,bool);
		virtual NodalValueBC *PasteNodalValue(NodalPoint *,TransportField *);
	
		// reaction flow calculation
		void InitQReaction(void);
        void InitQReaction(double scale);
		void SuperposeQReaction(double);
	
	protected:
		double qreaction;
	
};

#endif

