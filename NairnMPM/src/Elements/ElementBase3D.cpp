/********************************************************************************
    ElementBase3D.cpp
    nairn-mpm-fea
    
    Created by John Nairn on 7/20/06.
    Copyright (c) 2001 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"
#include "Elements/ElementBase3D.hpp"
#include "Nodes/NodalPoint.hpp"
#ifdef MPM_CODE
	#include "NairnMPM_Class/MeshInfo.hpp"
#endif

#pragma mark ElementBase3D::Constructors and Destructors

// main constructor - pass up the chain
ElementBase3D::ElementBase3D(int eNum,int *eNode) : ElementBase(eNum,eNode)
{
}

#pragma mark ElementBase3D::Methods

//	Find extent of this element - called once at start (and must be called)
void ElementBase3D::FindExtent(void)
{
	// find x-y extent in super class
	ElementBase::FindExtent();
	
    // find z extent of element
    int i,numnds=NumberNodes();
	double zNode;
    zmin=zmax=nd[nodes[0]]->z;
    for(i=1;i<numnds;i++)
    {   zNode=nd[nodes[i]]->z;
        if(zNode>zmax) zmax=zNode;
        if(zNode<zmin) zmin=zNode;
    }
	
    // set grid tolerance (1/10 minimum grid spacing)
    double range=TOLERANCE_RATIO*(zmax-zmin);
    if(range<gridTolerance) gridTolerance=range;
}

// Find center of mass of element (3D), and needed before extent is known
// If extent known, use GetXYZCentroid() instead
void ElementBase3D::FindCentroid(Vector *center) const
{
    cout << "Getting centroid" << endl;
    int i,numnds=NumberNodes();
    double xtot=nd[nodes[0]]->x;
    double ytot=nd[nodes[0]]->y;
    double ztot=nd[nodes[0]]->z;
    for(i=1;i<numnds;i++)
    {   xtot+=nd[nodes[i]]->x;
        ytot+=nd[nodes[i]]->y;
        ztot+=nd[nodes[i]]->z;
    }
    center->x=xtot/(double)numnds;
    center->y=ytot/(double)numnds;
    center->z=ztot/(double)numnds;
}

#pragma mark ElementBase3D::Accessors

// face nodes not meaningful
int ElementBase3D::FaceNodes(void) { return 0; }

// centroid (but possibly may not be in the element if it is distorted
void ElementBase3D::GetXYZCentroid(Vector *center)
{	center->x=(xmin+xmax)/2.;
	center->y=(ymin+ymax)/2.;
	center->z=(zmin+zmax)/2.;
}

// depth - 3D element return z extent
double ElementBase3D::GetDeltaZ(void) const { return zmax-zmin; }
bool ElementBase3D::IntersectsBox(Vector orig,double xlength,double ylength) const
{	if(orig.z<zmin) return false;
	if(orig.z>zmax) return false;
	return ElementBase::IntersectsBox(orig,xlength,ylength);
}
// 3D range (overides 2D in base)
void ElementBase3D::GetRange(int ax,double &amin,double &amax) const
{	if(ax==0)
	{	amin = xmin;
		amax = xmax;
	}
	else if(ax==1)
	{	amin = ymin;
		amax = ymax;
	}
	else
	{	amin = zmin;
		amax = zmax;
	}
}

// check if this GIMP element is on the edge of the grid
// assumes a generated 3D structured grid
bool ElementBase3D::OnTheEdge(void)
{	// now adds extra element for POINT_GIMP too and treats edge as left the grid
	//if(useGimp == POINT_GIMP) return FALSE;
	return mpmgrid.EdgeElement3D(num);
}

// for structured grid, return 0-terminated list of neighbors
void ElementBase3D::GetListOfNeighbors(int *theList) const { mpmgrid.ListOfNeighbors3D(num,theList); }


