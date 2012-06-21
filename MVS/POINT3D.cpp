#include <math.h>
#include <iomanip>
#include "POINT3D.h"

//Defnitions for POINT3D
POINT3D::POINT3D(const DP xx, const DP yy, const DP zz, const float rr, const float gg, const float bb, const float nx, const float ny, const float nz, const float maskk, const double cont, const int rid)
{
	x=xx;
	y=yy;
	z=zz;
	r=rr;
	g=gg;
	b=bb;
	normal_x=nx;
	normal_y=ny;
	normal_z=nz;
//	pde1 = pd1;
//	pde2 = pd2;
	mask = maskk;
	contour = cont;
	RegionId = rid;
}
	
POINT3D::POINT3D(const POINT3D &point) 
{
	if(this != &point)
	{
	    x=point.x;
	    y=point.y;
	    z=point.z;
		r=point.r;
		g=point.g;
		b=point.b;
		normal_x=point.normal_x;
		normal_y=point.normal_y;
		normal_z=point.normal_z;
		pde1 = point.pde1;
		pde2 = point.pde2;
		contour = point.contour;
		BCDArea = point.BCDArea;
	    averagedAbGradient=point.averagedAbGradient;
		mask = point.mask;
		RegionId = point.RegionId;
	}
}

POINT3D& POINT3D::operator=(const POINT3D temp)
{
	x=temp.x;
	y=temp.y;
	z=temp.z;
	r=temp.r;
	g=temp.g;
	b=temp.b;
	normal_x=temp.normal_x;
	normal_y=temp.normal_y;
	normal_z=temp.normal_z;
	pde1 = temp.pde1;
	pde2 = temp.pde2;
	contour = temp.contour;
	BCDArea = temp.BCDArea;
	averagedAbGradient=temp.averagedAbGradient;
	mask = temp.mask;
	RegionId = temp.RegionId;
	return *this;
}

POINT3D POINT3D::operator/(const DP c)
{
	DP d=c;	

	if (d == 0.0) d=TINY;
	return POINT3D(x/d,y/d,z/d,r,g,b,normal_x,normal_y,normal_z,mask,contour,RegionId);
	//operator/is only used for x y z
}	

POINT3D& POINT3D::operator/=(const DP c)
{
	*this=(*this)/c;
	return *this;
}
/*
GRID GRID::operator +(GRID g)
{
	return GRID(x+g.x,y+g.y,z+g.z,dist+g.dist,row+g.row,col+g.col,dep+g.dep,tag+g.tag);
}
*/