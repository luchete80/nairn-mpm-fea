/********************************************************************************
	LinearSoftening.hpp
	nairn-mpm-fea

	Created by John Nairn, Jan 21 2017.
	Copyright (c) 2017 John A. Nairn, All rights reserved.

	Dependencies
		SofteningLaw.hpp
********************************************************************************/

#ifndef _LinearSoftening_
#define _LinearSoftening_

#include "Materials/SofteningLaw.hpp"

class LinearSoftening : public SofteningLaw
{
	public:
		// methods
		virtual double GetFFxn(double,double) const;
		virtual double GetFpFxn(double,double) const;
		virtual double GetDDelta(double,double,double,double) const;
        virtual double GetDDeltaElastic(double,double,double,double,double) const;
		virtual double GetGToDelta(double,double) const;
		virtual double GetGoverGc(double,double) const;
		virtual double GetMaxSlope(double) const;
		virtual double GetDeltaFromDamage(double,double,double,double);
		virtual double GetRdFxn(double,double,double) const;
		virtual double GetPhiFxn(double,double) const;
	
		// accessors
		virtual const char *GetSofteningLawName(void) const;
		virtual double GetDeltaMax(double) const;
		virtual bool IsLinear(void) const;
		virtual double GetEtaStability(void) const;
};

#endif
