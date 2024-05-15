/********************************************************************************
    JohnsonCook.hpp
    nairn-mpm-fea
    
    Created by John Nairn, August 12, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		HardeningLawBase.hpp, MaterialBase.hpp
********************************************************************************/

#ifndef _GMTHARDENING_

#define _GMTHARDENING_

#define GMT_ID 10

#include "Materials/HardeningLawBase.hpp"

// plastic law properties
typedef struct {
	double hmlgTemp;
	double TjcTerm;
} GMTProperties;

class GMT : public HardeningLawBase
{
    public:
        // contructors
        GMT();
        GMT(MaterialBase *);
        
        // initialize
        virtual char *InputHardeningProperty(char *,int &,double &);
        virtual void PrintYieldProperties(void) const;
		virtual const char *VerifyAndLoadProperties(int);
	
		// copy of properties
        virtual int SizeOfHardeningProps(void) const;
		virtual void *GetCopyOfHardeningProps(MPMBase *,int,void *,int);
		virtual void DeleteCopyOfHardeningProps(void *,int) const;
    
        // hardening law core methods
        virtual double GetYield(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetKPrime(MPMBase *,int,double,HardeningAlpha *,void *) const;
        virtual double GetK2Prime(MPMBase *,double,double,HardeningAlpha *,void *) const;
        virtual double GetYieldIncrement(MPMBase *,int,double,HardeningAlpha *,void *) const;
	
		// return mapping
        virtual double SolveForLambdaBracketed(MPMBase *,int,double,Tensor *,double,double,
											   double,double,HardeningAlpha *a,void *,int) const;
    
        // accessors
        virtual const char *GetHardeningLawName(void) const;
    
    protected:
		double Bjc,Cjc,njc,ep0jc,Tmjc,mjc,Djc,n2jc;
    double n1gmt,n2gmt, m1gmt,m2gmt,C1gmt,C2gmt, I1gmt, I2gmt;
    double er_min, er_max, e_min, e_max, T_min, T_max;  
		double C1gmt_red,edotMin,eminTerm;

};

#endif
