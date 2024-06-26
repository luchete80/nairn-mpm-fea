/********************************************************************************
    ShapeController.hpp
    NairnFEA
    
    Created by John Nairn on 8/8/07.
    Copyright (c) 2007 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _SHAPECONTROLLER_

#define _SHAPECONTROLLER_

#include "System/MPMPrefix.hpp"

class CommonReadHandler;
#ifdef MPM_CODE
	class MPMBase;
#endif

class ShapeController
{
	public:
	
		// contructors
		ShapeController(int);
		ShapeController(int,double,double,double,double);
		virtual ~ShapeController();
    
        // initialize
        virtual void SetProperty(const char *,char *,CommonReadHandler *);
        virtual void SetProperty(const char *,double);
        virtual void SetProperty(char *,CommonReadHandler *);
        virtual void SetParameter(const char *,const char *);
        virtual bool FinishParameter(void);
        virtual bool FinishSetup(void);
        virtual bool HasAllParameters(void);
	
		// initialize, but do not override without converting to virtual
		void SetScaling(double);
	
		// ShapeController methods only (non-virtual)
		bool ShapeContainsPoint(Vector&);
		void resetNodeEnumerator(void);
		void AddCutoutShape(ShapeController *);
	
		// methods
		virtual bool ContainsPoint(Vector &);
		virtual const char *startNodeEnumerator(int,int);
		virtual int nextNode(void);
        void resetElementEnumerator(void);
        int nextElement(void);
	
        // MPM only methods
#ifdef MPM_CODE
		void setNetBC(bool);
		double particleCount(void);
		virtual void resetParticleEnumerator(void);
		virtual int nextParticle(void);

#endif
    
        // accessors
        virtual const char *GetShapeName(void);
		virtual bool IsRealShape(void);
		virtual void DescribeShape(const char *);
        virtual bool Is2DShape(void);
        virtual char *GetContextInfo(void);
        // base class only (non virtual)
        int GetSourceBlock(void);
        bool RequiredBlock(int);
		ShapeController *GetParentShape(void) const;
		void SetParentShape(ShapeController *);
	
	protected:
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double distScaling;
		int sourceBlock,nodeNum,elemNum;
		bool twoDShape;
#ifdef MPM_CODE
		int particleNum,numParticles;
#endif
		vector< ShapeController * > children;
		ShapeController *parentShape;
	
};

extern ShapeController *theShape;

#endif

