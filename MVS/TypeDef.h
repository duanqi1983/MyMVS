#ifndef __TYPEDEF_H
#define __TYPEDEF_H

using namespace std;
#include <vector>
#include "POINT3D.h"
#include <ctime>

typedef vector<vector<vector<VECTOR3D> > > vec3DVec3d;
typedef vector<vector<VECTOR3D> > vec2DVec3d;
typedef vector<VECTOR3D> vec1DVec3d;
typedef vector<POINT3d> vecPoints;
typedef vector<GRID> vecGrids;
typedef vector<vector<vector<GRID> > > vec3DGRIDs;//三维的GRID数组。
typedef vector<vector<GRID> > vec2DGRIDs;//二维的GRID数组
typedef vector<GRID> vec1DGRIDs;//一维的GRID数组
typedef	vector<float> vecFloats;
typedef vector<double> vecDoubles;
typedef vector<vector<int> > vecTriangles;
typedef vector<int> vecInts;
typedef vector<ONEDISK> vecOneDisks2;//not used yet.
typedef vector<vecInts> vecOneDisks;
typedef vector<DISKTRIANGLE> DISK;
typedef vector<DISK> vecONEDISKS;
typedef vector<DISKTRIANGLEANISOTROPIC> DISKANISOTROPIC;
typedef vector<DISKANISOTROPIC> vecONEDISKSANISOTROPIC;
typedef vector<ONERINGVERTEX> ONERING;
typedef vector<ONERING> vecONERINGS;
typedef vector<TRIANGLE> vecTriangles2;
typedef vector<TRIANGLEWITHBC> vecTrianglesWithBC;
typedef vector<VERTEX3D> vecVertex3D;
typedef vector<TRIANGLEPPIBGRADIENT> vecTPPIBG;
typedef vector<COLORRGB> vecCOLORRGB;

#define EPS 1.0e-14
#define INFINITY (100000.0)
#define SURFACE 1
#define DISTANCEFUNCTION 2
#define DISTANCEFUNCTION_DATA 3
#define PI 3.1415926535
#define ALPHA 0.00000001
#define OMEGA 0.5
#define CFL 0.05
#define DISPLAYEPSILON 1.0//原来为0.75.
#define SOLVEBAND 5.0
#define INITIALIZEDBAND 1.25//原来为1.05
#define CONTOURVIEW 12.0
//the following grids is for sphere data:(没有溢出，最小梯度模超过0.5f.)
//#define NX (148) //maybe 128
//#define NY (148) //maybe 128
//#define NZ (148) //maybe 128
//the following grids is for Bretzel2 data:
#define NX (240)
#define NY (144)
#define NZ (80)
//the following grids is for horse data:
//#define NX (158) //m_extnx
//#define NY (135) //m_extny
//#define NZ (81)  //m_extnz
//the following grids is for horse data(high resolution)://
//#define NX (220) //m_extnx
//#define NY (187) //m_extny
//#define NZ (112)  //m_extnz
//the following grids is for venus_high_res data:
//#define NX (115) //m_extnx
//#define NY (91) //m_extny
//#define NZ (216)  //m_extnz
//the following grids is for vase data:
//#define NX (123) //m_extnx
//#define NY (205) //m_extny
//#define NZ (118)  //m_extnz
//the following grids is for vase_new data://
//#define NX (216) //m_extnx
//#define NY (130) //m_extny
//#define NZ (125)  //m_extnz

#define NT 10
#define EXTGRIDNUMBER 8 //每个方向的边界处扩展的网格数，有利于高精度的求解
#define MIN(a,b)  (((a)<=(b))?(a):(b))
#define MAX(a,b)  (((a)>=(b))?(a):(b))
#define MAX3(a,b,c) (MAX(MAX(a,b),c))
#define MIN3(a,b,c) (MIN(MIN(a,b),c))
#define POWER(x)  ((x)*(x))
#define POWER3(x) ((x)*(x)*(x))
#define POWER4(x) ((x)*(x)*(x)*(x))
#define SIGN(a)   (((a)>0)? 1: (((a)<0)? -1:0))

#endif