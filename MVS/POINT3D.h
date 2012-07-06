#ifndef _POINT_H_
#define _POINT_H_
#include <limits>
#include <iostream>
#include <vector>

using namespace std;

typedef double DP;	// Coincide with the definiton in "vecmat.h"
const DP TINY=numeric_limits<DP>::epsilon();

struct VECTOR3D{
	double x;
	double y;
	double z;
};
typedef VECTOR3D VECTOR3D;
typedef VECTOR3D VERTEX3D;
typedef VECTOR3D DUALVERTEX3D;

class POINT3D {
public:
	DP x,y,z;
	double r,g,b;//为该三维点定义颜色
	double normal_x,normal_y,normal_z;
	VECTOR3D pde1,pde2;//2 principle directions...
	double divergence;// added on 28/12/2009.
	double BCDArea;
	double averagedAbGradient;//averaged absolute gradient;
	double contour;
	float mask;//mask for inpainting images.
	bool IsImageMapped;
	bool IsOnBoundary;//是否为边界点？
	int RegionId;
	POINT3D(const DP xx=0., const DP yy=0., const DP zz=0., const float rr=0.,const float gg=0.,const float bb=0.,const float nx=1.0,const float ny=1.0,const float nz=1.0,const float maskk=1.0, const double cont=0.0, const int rid=1);
	POINT3D(const POINT3D &point);

	POINT3D & operator=(const POINT3D q);
	POINT3D operator/(const DP c);
	POINT3D & operator/=(const DP c);

	// add by duan qi for alm_mesh_refinement
	double ref_x, ref_y, ref_z, ref_light_x, ref_light_y, ref_light_z;
	double ps_normal_x, ps_normal_y, ps_normal_z;
<<<<<<< HEAD
	double light_r, light_theta, light_phi; // sphere coordinate, x=r*sin(theta)*cos(phi), y=r*sin(theta)*sin(phi), z = r*cos(phi)
=======
	double light_x, light_y, light_z;
>>>>>>> Sphere-Coordinates-Light
	std::vector<double> intensity_list;
	std::vector<VECTOR3D> Camera_center_list;
	double tar_x, tar_y, tar_z;
	double ps_weight; double intensity;
	double x_divergence, y_divergence, z_divergence;
	double Voronoi_Area;
};

typedef POINT3D VECTOR;
typedef POINT3D POINT3d;//原来是typedef POINT3D POINT;

struct TRIANGLE{
	int ver0;
	int ver1;
	int ver2;//三个顶点序号
	VECTOR3D grad;//the gradient of current quantity in this triangle, it's a constant vector.
	VECTOR3D grad1;
	VECTOR3D grad2;
	double edgeindicator;//the edge indicator of this triangle, used in edge detection of images.
};
typedef TRIANGLE TRIANGLE;
typedef TRIANGLE ORIENTEDTRIANGLE;


struct EDGE{
	int ver0;
	int ver1;
};
typedef EDGE EDGE;

struct ORIENTEDEDGE{
	int begin;
	int end;
};
typedef ORIENTEDEDGE ORIENTEDEDGE;

struct INDEXLIST{
	int index;
	INDEXLIST *pNext;
};
typedef INDEXLIST INDEXLIST;

struct ONEDISK{
	int vertexIndex;
	bool vertexIsOnBoundary;
	INDEXLIST *pTriList;
};
typedef ONEDISK ONEDISK;

struct DISKTRIANGLE{
	int itri;//该disk中此三角片的序号。
	float ver0co;//lb算子在本三角片的ver0顶点处的组合系数。
	float ver1co;//lb算子在本三角片的ver1顶点处的组合系数。
	float ver2co;//lb算子在本三角片的ver2顶点处的组合系数。
};
typedef DISKTRIANGLE DISKTRIANGLE;

struct DISKTRIANGLEANISOTROPIC{
	int iDiskCent;//disk center vertex index.
	int itri;
	float d1ver0co;//各向异性耗散对应d1方向的v0系数。
	float d1ver1co;//各向异性耗散对应d1方向的v1系数。
	float d1ver2co;//各向异性耗散对应d1方向的v2系数。
	float d2ver0co;//各向异性耗散对应d2方向的v0系数。
	float d2ver1co;//各向异性耗散对应d2方向的v1系数。
	float d2ver2co;//各向异性耗散对应d2方向的v2系数。
};
typedef DISKTRIANGLEANISOTROPIC DISKTRIANGLEANISOTROPIC;

struct ONERINGVERTEX{
	int iver;//一环中某个顶点的序号。
	double co;//lb算子在该顶点上的组合系数。
};
typedef ONERINGVERTEX ONERINGVERTEX;

struct GRID{
	//int row,col,dep;
	int image_x,image_y;//距离最近的图象像素格点。(-1,-1)for not image mapped。
	double image_dist;//与映射好的图象曲面片之间的距离。
	//double x,y,z;
	float r,g,b;//
	bool toMapImage;//是否需要映射图象像素。
	bool isImageMapped;//表明此网格点是否已经映射了像素。
	bool isInitialized;//是否初始化过，指曲面数据附近
	int tag;//-1 for in;0 for on;+1 for out
	double dist;//the distance from the grid to the given data
	float ramda;/*inpainting indicator:0 for inpaintng,1 for uninpainting */
	double contour;//用来提取图像边缘的网格属性，实际上是一个定义于曲面网格上的
	//隐函数，其0等值集表示图像的边缘。
};

typedef GRID GRID;


struct DUALEDGE{
	VERTEX3D ver0;
	VERTEX3D ver1;
};
typedef DUALEDGE DUALEDGE;
typedef DUALEDGE EDGEWITHCOOR;
struct ORIENTEDDUALEDGE{
	VERTEX3D begin;
	VERTEX3D end;
};
typedef ORIENTEDDUALEDGE ORIENTEDDUALEDGE;
typedef ORIENTEDDUALEDGE ORIENTEDEDGEWITHCOOR;
struct TRIANGLEWITHBC{
	int ver0;
	int ver1;
	int ver2;
	DUALVERTEX3D bc;
};
typedef TRIANGLEWITHBC TRIANGLEWITHBC;
typedef TRIANGLEWITHBC ORIENTEDTRIANGLEWITHBC;
struct TRIANGLEWITHCOOR{
	VERTEX3D ver0;
	VERTEX3D ver1;
	VERTEX3D ver2;
};
typedef TRIANGLEWITHCOOR TRIANGLEWITHCOOR;
struct TRIANGLEPPIBGRADIENT{//一个三角形内的primal-primal插值基函数的梯度,每个顶点一个。
	VECTOR3D v0;
	VECTOR3D v1;
	VECTOR3D v2;
};
typedef TRIANGLEPPIBGRADIENT TRIANGLEPPIBGRADIENT;

struct FLOATSTRIPLE3{
	float r;
	float g;
	float b;
};
typedef FLOATSTRIPLE3 COLORRGB;

struct REACTIONDIFFUSIONTEXTUREPARAM{
	double Da;//diffusion rate of chemical a;
	double Db;//diffusion rate of chemical b;
	double ar;//action rate between a and b;
	double gr;//growth rate of a;
	double decay;//decay of b;
};
typedef REACTIONDIFFUSIONTEXTUREPARAM RDTEXTUREPARAM;

struct REACTIONDIFFUSIONTEXTUREPARAM_MeinhardtStrips{
	double c;
	double alpha;
	double Dg;
	double tho0;
	double beta;
	double gamma;
	double Ds;
	double tho1;
};
typedef REACTIONDIFFUSIONTEXTUREPARAM_MeinhardtStrips REACTIONDIFFUSIONTEXTUREPARAM_MeinStrips;

struct FRAME{
	//a frame is a orthornormal basis of the space at current position.
	POINT3d p;
	VECTOR3D e1;
	VECTOR3D e2;
	VECTOR3D e3;
};
typedef FRAME FRAME;

struct MATRIX3BY3{
	double a11;double a12;double a13;
	double a21;double a22;double a23;
	double a31;double a32;double a33;
};
typedef MATRIX3BY3 MATRIX3BY3;
typedef MATRIX3BY3 DIFFUSIONTENSOR;

struct FABDIFFUSIONPARAM{
	float k_fore;//相对于most absolute intrinsic gradient的倍数.
	float k_back;//相对于most absolute intrinsic gradient的倍数.
	float bandwidth;//相对于most absolute intrinsic gradient的倍数.
	unsigned int nforee;//指数.
	unsigned int nbacke;//指数.
	float alpharatio;//ratio.
};
typedef FABDIFFUSIONPARAM FABDIFFUSIONPARAM;

struct GRID_WITH_COORD{
	GRID g;
	VECTOR3D p;
};
typedef GRID_WITH_COORD GRID_WITH_COORD;

struct CUBE{
	GRID_WITH_COORD c,e,n,en,u,ue,un,uen;
	VECTOR3D nc,ne,nn,nen,nu,nue,nun,nuen;
};
typedef CUBE CUBE;

struct PLANE{
	//给出两种描述平面的方法：4系数法和点法向量法。
	double A;
	double B;
	double C;
	double D;
	VECTOR3D p;//平面上一点
	VECTOR3D n;//法向量
};
typedef PLANE PLANE;

#endif //for _POINT_H_
