// TriangularMesh.cpp: implementation of the TriangularMesh class.
//
//////////////////////////////////////////////////////////////////////

#include "TriangularMesh.h"
#include <time.h>
#include "speigen.h"

//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define DEBUG_NEW new(THIS_FILE, __LINE__) 
//#define new DEBUG_NEW
//#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TriangularMesh::TriangularMesh()
{
	m_nTime = 1;//50 for reaction-diffusion texturing.
	m_Lambda = 1000;//1.
	m_alpha = 0.00001;//for TV diffusion term, to avoid zero division, 原为0.01.
	m_tStep = 0.00005; // 0.5 for reaction-diffusion texturing.
	                  // 0.00005 for denoising. if the noise is high, timestep can be larger to speed up.
	                  // time steps for inpainting are always large.
	                  // 0.00005 for Anisotropic scale spaces.
	                  // 0.0005 for Aniso SS antialiasing + denoising.
	                  // 0.002 for GCF scale spaces.
	                  // 0.01 for Perona-Malik-I scale spaces.
	                  // 0.03 for BFB scale spaces.
	                  // 0.001 for ITV scale spaces.
	                  // 0.0005 for LHF scale spaces.
	                  // 0.05 for edge detection (0.01 for the example of LeafonHorse).
	SetRDTextureParam(1.0,0.0625,0.025,16,12);
	m_ObjectDrawing = 2;// to show the mesh and the contour.
}

TriangularMesh::~TriangularMesh()
{

}

double TriangularMesh::DotProduct(VECTOR3D v1, VECTOR3D v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

VECTOR3D TriangularMesh::CrossProduct(VECTOR3D v1, VECTOR3D v2)
{
	VECTOR3D res;
	res.x = v1.y*v2.z - v1.z*v2.y;
    res.y = v1.z*v2.x - v1.x*v2.z;
	res.z = v1.x*v2.y - v1.y*v2.x;
	return res;
}

bool TriangularMesh::LoadFromFile(const char* filename)
{
	ifstream ifp;
	ifp.open(filename);
	if(ifp.fail())
	{
		ifp.close();
		return false;
	}
	POINT3d tmp;
	TRIANGLE tmp2;
	COLORRGB crgb;
	int feature_id,ss,tt;//
	double s0,t0,s1,t1,s2,t2;//variables for loading data models from Turk.
	//ifp.open(filename);
	m_vertices.clear();
	m_triangles.clear();
	m_tempRGB.clear();
	m_tempContour.clear();

	m_WeightedCurveLength.clear();//
	m_tSteps.clear();//
	m_WeightedCurveLength_Forobservation.clear();//
	m_tSteps_Forobservation.clear();//
	m_CPUtime.clear();//

	
	int i; int poly;
	double aaa;
	ifp>>m_vnum>>m_trinum;
    for(i=0;i<m_vnum;i++)
	{
	   ifp>>tmp.x>>tmp.y>>tmp.z>>tmp.r>>tmp.g>>tmp.b;//for textured meshes.
//	   ifp>>tmp.x>>tmp.y>>tmp.z>>tmp.r>>tmp.g>>tmp.b>>tmp.contour;//for textured meshes with contour data.
//	   ifp>>tmp.x>>tmp.y>>tmp.z>>tmp.contour;// for meshes with only contour data.
//     ifp>>tmp.x>>tmp.y>>tmp.z>>tmp.r>>tmp.g>>tmp.b>>tmp.mask;//for textured distorted meshes.
//	   ifp>>tmp.x>>tmp.y>>tmp.z;
//	   tmp.mask = 1.0f;
//	   tmp.z = 0;
	   tmp.RegionId = 0;
	   m_vertices.push_back(tmp);
	}

//	Contour2Intensity(0.3,0.9);
/*	for(i=0;i<m_vnum;i++)
		if(m_vertices[i].x<-0.3)
		{
			m_vertices[i].r = 0.8;
			m_vertices[i].g = 0.8;
			m_vertices[i].b = 0.8;
		}
		else if(m_vertices[i].x<0)
		{
			m_vertices[i].r = 0.6;
			m_vertices[i].g = 0.6;
			m_vertices[i].b = 0.6;
		}
		else if(m_vertices[i].x<0.3)
		{
			m_vertices[i].r = 0.2;
			m_vertices[i].g = 0.2;
			m_vertices[i].b = 0.2;
		}
		else
		{
			m_vertices[i].r = 0.6;
			m_vertices[i].g = 0.6;
			m_vertices[i].b = 0.6;
		}*/
/*  if the mesh is not textured, color it with initial color.
    for(i=0;i<m_vnum;i++)
	{
		if(m_vertices[i].r>=0.50)
		{
			m_vertices[i].r = 0.5;
            m_vertices[i].g = m_vertices[i].r;
            m_vertices[i].b = m_vertices[i].r;
		}
		else
		{
			m_vertices[i].r = 0.1;
            m_vertices[i].g = m_vertices[i].r;
            m_vertices[i].b = m_vertices[i].r;
		}
	}
*/	

	for(i=0;i<m_trinum;i++)
	{
//		ifp>>poly>>tmp2.ver0>>tmp2.ver1>>tmp2.ver2;//for venus_head and goodbunny.
		ifp>>tmp2.ver0>>tmp2.ver1>>tmp2.ver2;//for horse.
		//ifp>>poly>>tmp2.ver0>>tmp2.ver1>>tmp2.ver2
		//	>>feature_id>>ss>>s0>>tt>>t0
		//	>>ss>>s1>>tt>>t1>>ss>>s2>>tt>>t2;
		m_triangles.push_back(tmp2);
	}

	ifp.close();
	UnifyData();
/*
	for(i=0;i<m_vnum;i++)
		if((m_vertices[i].z)>0.25)
		{
            m_vertices[i].r = 0.1;
            m_vertices[i].g = m_vertices[i].r;
            m_vertices[i].b = m_vertices[i].r;
		}
		else if((m_vertices[i].z)>0)
		{
            m_vertices[i].r = 0.4;
            m_vertices[i].g = m_vertices[i].r;
            m_vertices[i].b = m_vertices[i].r;
		}
        else if((m_vertices[i].z)>-0.25)
		{
			m_vertices[i].r = 0.7;
            m_vertices[i].g = m_vertices[i].r;
            m_vertices[i].b = m_vertices[i].r;
		}
		else
		{
			m_vertices[i].r = 0.9;
            m_vertices[i].g = m_vertices[i].r;
            m_vertices[i].b = m_vertices[i].r;
		}
*/		
//	AddSaltPepperNoise(0.1);
//	AddUniformNoise(0.1);
	AddGaussianNoise(0.05);
//  for gray image .............
	for(i=0;i<m_vnum;i++)
	{
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}
//	AddGaussianNoise2VertexGeometry(0.002);

/*  the following is for initial contouring. */
//    SetInitialContourFunction();

//	WriteData();

//the following color and mask configuration is for inpainting.
/*
    for(i=0;i<m_trinum;i++)
	{
		if(i%240==0)
		{
			tmp2 = m_triangles[i];
			m_vertices[tmp2.ver0].r = 0.5;
			m_vertices[tmp2.ver0].g = 0.5;
			m_vertices[tmp2.ver0].b = 0.5;
			m_vertices[tmp2.ver0].mask = 0.0f;
			m_vertices[tmp2.ver1].r = 0.5;
			m_vertices[tmp2.ver1].g = 0.5;
			m_vertices[tmp2.ver1].b = 0.5;
			m_vertices[tmp2.ver1].mask = 0.0f;
			m_vertices[tmp2.ver2].r = 0.5;
			m_vertices[tmp2.ver2].g = 0.5;
			m_vertices[tmp2.ver2].b = 0.5;
			m_vertices[tmp2.ver2].mask = 0.0f;
		}
	}
	WriteData();*/
/*	
	int random;
	for(i=0;i<300;i++)
	{
		random = int(double(rand())/RAND_MAX*m_trinum);
		tmp2 = m_triangles[random];
		m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		m_vertices[tmp2.ver0].mask = 0.0f;
		m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
		m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		m_vertices[tmp2.ver1].mask = 0.0f;
		m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		m_vertices[tmp2.ver2].mask = 0.0f;

		if(random+1024<m_trinum && random-1024>0)
		{
		    tmp2 = m_triangles[random+512];//512 is the total number of a row triangles.
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
		    tmp2 = m_triangles[random+1];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random-1];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random-512];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random-1024];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random+1024];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random-513];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random-511];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random+513];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
			tmp2 = m_triangles[random+511];
		    m_vertices[tmp2.ver0].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver0].g = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].b = m_vertices[tmp2.ver0].r;
		    m_vertices[tmp2.ver0].mask = 0.0f;
		    m_vertices[tmp2.ver1].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver1].g = m_vertices[tmp2.ver1].r;
	 	    m_vertices[tmp2.ver1].b = m_vertices[tmp2.ver1].r;
		    m_vertices[tmp2.ver1].mask = 0.0f;
		    m_vertices[tmp2.ver2].r = (double)(rand())/RAND_MAX;
		    m_vertices[tmp2.ver2].g = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].b = m_vertices[tmp2.ver2].r;
		    m_vertices[tmp2.ver2].mask = 0.0f;
		}
	}
*/
//  initial guess for inpainting.
/*
    for(i=0;i<m_vnum;i++)
		if(m_vertices[i].mask==0.0)
		{
			m_vertices[i].r = (double)(rand())/RAND_MAX;
			m_vertices[i].g = m_vertices[i].r;
			m_vertices[i].b = m_vertices[i].r;
		}*/
//  WriteData();//
	for(i=0;i<m_vnum;i++)
	{
		crgb.r = m_vertices[i].r;
        crgb.g = m_vertices[i].g;
        crgb.b = m_vertices[i].b;
		m_tempRGB.push_back(crgb);
	}
/*
	for(i=0;i<m_vnum;i++)
	{
		m_tempContour.push_back(m_vertices[i].contour);
	}
*/
	for(i=0;i<m_vnum;i++)
	{//initialize the normals of vertices.
		m_vertices[i].normal_x = 0;
		m_vertices[i].normal_y = 0;
		m_vertices[i].normal_z = 0;
	}
	FindNormals();
//	AddGaussianNoise2VertexNormal(0.2);

//
//  First of all we compute the mesh information indicating its quality.
//	GetMeshInformation();
//
	BuildOneDisks();
//  BuildCCDualVerticesForTriangles();//this is to construct cc dual.
	BuildBaryCentersForTriangles();//this is to construct bc dual.
	BuildTPPIBG();//to compute the gradient of basis functions of primal-primal interpolation.
//  WriteTPPIBG();
	BuildONEDISKS();
	BuildONEDISKSANISOTROPIC();
	ComputeBCDualArea();
	BuildOnerings();
	BuildTrianglesArea();

/*  firstly blur the image to test the contour detection algorithm for blurred images.
	
	double temp;
	temp = m_tStep;
	m_tStep = 0.0001;
	ImplicitL2Diff();// blur the image.
	m_tStep = temp;
*/
	BuildEdgeIndicator(2);

	m_WeightedCurveLength.push_back(ComputeWeightedCurveLength());
	m_tSteps.push_back(m_tStep);
	m_WeightedCurveLength_Forobservation.push_back(ComputeWeightedCurveLength());
	m_tSteps_Forobservation.push_back(m_tStep);

//  if one uses the time step adaptive GCF, then the non-adaptive GCF should be called once before as following.
//    ImplicitGeodesicCurvatureFlow();

//  if one uses the time step adaptive WGCF, then the non-adaptive WGCF should be called once before as following.
//  ImplicitGeodesicWeightedCurvatureFlow();

//	ImplicitL2DiffDirecData();//blending the normals.
//	BuildPrincipleDirections();

	bool hasOA = HasObtuseAngle();

	return true;
}

void TriangularMesh::UnifyData()
{//归一化物体坐标，使得原点在物体中心，
 //物体大小为[-0.5f,0.5f]*[-0.5f,0.5f]*[-0.5f,0.5f]
	double minx,maxx,miny,maxy,minz,maxz;
	
	minx=maxx=m_vertices[0].x;
	miny=maxy=m_vertices[0].y;
	minz=maxz=m_vertices[0].z;
	
	int i;

	for(i=1;i<m_vnum;i++) {
		if(minx>m_vertices[i].x) minx=m_vertices[i].x;
		if(maxx<m_vertices[i].x) maxx=m_vertices[i].x;
		if(miny>m_vertices[i].y) miny=m_vertices[i].y;
		if(maxy<m_vertices[i].y) maxy=m_vertices[i].y;
		if(minz>m_vertices[i].z) minz=m_vertices[i].z;
		if(maxz<m_vertices[i].z) maxz=m_vertices[i].z;
	}
	
	double midx=0.5*(minx+maxx), midy=0.5*(miny+maxy), midz=0.5*(minz+maxz);
	for(i=0;i<m_vnum;i++) {
		m_vertices[i].x-=midx;
		m_vertices[i].y-=midy;
		m_vertices[i].z-=midz;
	}
	
	double widx=maxx-minx, widy=maxy-miny, widz=maxz-minz;
	if(widz<widy) widz=widy;
	if(widz<widx) widz=widx;
	for(i=0;i<m_vnum;i++) m_vertices[i]/=widz;
}

void TriangularMesh::Draw()
{
//
	switch(m_ObjectDrawing){
	case 0:
		DrawMesh();
		break;
	case 1:
	    DrawContour_MulRegionLabelling();
		break;
	case 2:
		DrawMesh();
		DrawContour_MulRegionLabelling();
		break;
	default:
		DrawMesh();
		break;
	}
//
}

void TriangularMesh::FindNormals()
{
	VECTOR3D trinormal,v1,v2;
	int i;
	for(i=0;i<m_trinum;i++)
	{
		v1.x = m_vertices[m_triangles[i].ver1].x-m_vertices[m_triangles[i].ver0].x;
		v1.y = m_vertices[m_triangles[i].ver1].y-m_vertices[m_triangles[i].ver0].y;
		v1.z = m_vertices[m_triangles[i].ver1].z-m_vertices[m_triangles[i].ver0].z;
		v2.x = m_vertices[m_triangles[i].ver2].x-m_vertices[m_triangles[i].ver0].x;
		v2.y = m_vertices[m_triangles[i].ver2].y-m_vertices[m_triangles[i].ver0].y;
		v2.z = m_vertices[m_triangles[i].ver2].z-m_vertices[m_triangles[i].ver0].z;
		trinormal = CrossProduct(v1,v2);
		m_vertices[m_triangles[i].ver0].normal_x+= trinormal.x;
		m_vertices[m_triangles[i].ver0].normal_y+= trinormal.y;
		m_vertices[m_triangles[i].ver0].normal_z+= trinormal.z;
		m_vertices[m_triangles[i].ver1].normal_x+= trinormal.x;
		m_vertices[m_triangles[i].ver1].normal_y+= trinormal.y;
		m_vertices[m_triangles[i].ver1].normal_z+= trinormal.z;
		m_vertices[m_triangles[i].ver2].normal_x+= trinormal.x;
		m_vertices[m_triangles[i].ver2].normal_y+= trinormal.y;
		m_vertices[m_triangles[i].ver2].normal_z+= trinormal.z;
	}
	
	for(i=0;i<m_vnum;i++)
	{
		v1.x = m_vertices[i].normal_x;
		v1.y = m_vertices[i].normal_y;
		v1.z = m_vertices[i].normal_z;
		NormalizeVector3D(&v1);
		m_vertices[i].normal_x = v1.x;
		m_vertices[i].normal_y = v1.y;
		m_vertices[i].normal_z = v1.z;
	}
}

void TriangularMesh::NormalizeVector3D(VECTOR3D *v)
{
	double l;
	l = sqrt((*v).x*(*v).x+(*v).y*(*v).y+(*v).z*(*v).z);
	if(l>0)
	{
		(*v).x/= l;
		(*v).y/= l;
		(*v).z/= l;
	}
}

bool TriangularMesh::WriteNormalsToFile()
{
	ofstream ofp("horse_normals.txt");
	//ofstream ofp("venushead_normals.txt");
	int i;
	ofp<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
	{
		ofp<<m_vertices[i].normal_x<<"   "<<m_vertices[i].normal_y<<"   "<<m_vertices[i].normal_z<<"\n";
	}
	ofp.close();
	return true;
}

void TriangularMesh::BuildOneDisks()
{
	int i,j,k,tmpindex;
	EDGE e;
	int v0;
	vecInts ind;
	ind.clear();
	m_onedisks.clear();
	for(i=0;i<m_vnum;i++)
	{
		m_onedisks.push_back(ind);
	}
    for(i=0;i<m_trinum;i++)
	{
		m_onedisks[m_triangles[i].ver0].push_back(i);
        m_onedisks[m_triangles[i].ver1].push_back(i);
		m_onedisks[m_triangles[i].ver2].push_back(i);
	}
/*	
	ofstream ofp1("onedisks.txt");
	ofp1<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
	{
		for(j=0;j<m_onedisks[i].size();j++)
		{
		    ofp1<<m_onedisks[i][j]<<",";
		}
		ofp1<<"\n";
	}
	ofp1.close();
*/	
	// one disks for all vertices are built, but they may be neither counterclockwise nor clockwise.
	// in the following we adjust them to be clockwise or counterclockwise order.
	for(i=0;i<m_vnum;i++)
	{
/****    bubble-sort like adjustment for each one-disk.               ********
 ****	 but at last we need tell if the one disk is a circle,             ***
 ****    this helps to judge if the vertex i is on the boundary of the mesh.**
 */
		v0 = FindAnotherVertex(m_triangles[m_onedisks[i][0]],i);
		e.ver0 = i;e.ver1 = v0;
		for(j=1;j<m_onedisks[i].size();j++)
		{
			for(k=j;k<m_onedisks[i].size();k++)
			{
				if(IsAdjointTriangle2(m_triangles[m_onedisks[i][j-1]],e,m_triangles[m_onedisks[i][k]]))
				{//there exists an really adjoint triangle.
					tmpindex = m_onedisks[i][k];
					m_onedisks[i][k] = m_onedisks[i][j];
					m_onedisks[i][j] = tmpindex;
					e = FindAdjointCoTriEdge(m_triangles[m_onedisks[i][j]],e,e.ver0);
					break;
				}
			}
		}
		//e = FindAdjointCoTriEdge(m_triangles[m_onedisks[i][m_onedisks[i].size()-1]],e,e.ver0);//this line is not performed for the last triangle in the above loop.
		if(!IsAdjointTriangle2(m_triangles[m_onedisks[i][0]],e,m_triangles[
			m_onedisks[i][m_onedisks[i].size()-1]]))//boundary vertex.
			m_vertices[i].IsOnBoundary = true;
		else
			m_vertices[i].IsOnBoundary = false;
	}
/*
	ofstream ofp2("orderedonedisks.txt");
	ofp2<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
	{
		for(j=0;j<m_onedisks[i].size();j++)
		{
		    ofp2<<m_onedisks[i][j]<<",";
		}
//		if(m_vertices[i].IsOnBoundary)
//			ofp2<<" On Boundary!";
		ofp2<<"\n";
	}
	ofp2.close();
*/
}

bool TriangularMesh::IsSharingEdge(TRIANGLE t1, TRIANGLE t2)
{
	// not realized yet.
	return false;
}

bool TriangularMesh::HasThisEdge(TRIANGLE t, EDGE e)
{
	if((e.ver0==t.ver0)&&(e.ver1==t.ver1||e.ver1==t.ver2))
		return true;
	if((e.ver0==t.ver1)&&(e.ver1==t.ver0||e.ver1==t.ver2))
		return true;
	if((e.ver0==t.ver2)&&(e.ver1==t.ver0||e.ver1==t.ver1))
		return true;
	return false;
}

int TriangularMesh::FindAnotherVertex(TRIANGLE t, int v0)
{//here v0 is one of the 3 vertices of t;
 //if not, return -1.
	if(t.ver0==v0)
		return t.ver1;
	if(t.ver1==v0)
		return t.ver0;
	if(t.ver2==v0)
		return t.ver0;
	return -1;
}

int TriangularMesh::FindOppositeVertex(TRIANGLE t, EDGE e)
{// to find the vertex oppositing the given edge e.
	if((e.ver0==t.ver0&&e.ver1==t.ver1)||(e.ver0==t.ver1&&e.ver1==t.ver0))
		return t.ver2;
	if((e.ver0==t.ver1&&e.ver1==t.ver2)||(e.ver0==t.ver2&&e.ver1==t.ver1))
		return t.ver0;
	if((e.ver0==t.ver0&&e.ver1==t.ver2)||(e.ver0==t.ver2&&e.ver1==t.ver0))
		return t.ver1;
	return -1;
}

EDGE TriangularMesh::FindAdjointCoTriEdge(TRIANGLE t, EDGE e, int v0)
{
	EDGE ret;
	ret.ver0 = v0;
	ret.ver1 = FindOppositeVertex(t,e);
	return ret;
}

bool TriangularMesh::HasThisEdge2(TRIANGLE t, EDGE e)
{//judge under the condition that e.ver0 is a vertex of t.
	if(e.ver1==t.ver0||e.ver1==t.ver1||e.ver1==t.ver2)
		return true;
	return false;
}

bool TriangularMesh::IsAdjointTriangle(TRIANGLE t, EDGE e, TRIANGLE at)
{//calling HasThisEdge(at,e); not applied up to now.
	return false;
}

bool TriangularMesh::IsAdjointTriangle2(TRIANGLE t, EDGE e, TRIANGLE at)
{//calling HasThisEdge2(at,e); return the really adjoint triangle of t which shares the edge e.
	if(HasThisEdge2(at,e))
	{
		if(FindOppositeVertex(t,e)!=FindOppositeVertex(at,e))
			return true;
		else return false;
	}
	return false;
}

VERTEX3D TriangularMesh::GetCCDual(TRIANGLE t)
{
/* 
 *    return the dual of triangle t.
 *    compute it through edges (ver0ver1) and (ver0ver2).
 */
	VERTEX3D dv;
	VECTOR3D v1,v2,v3,v4;
	VERTEX3D c1,c2;
	double asp;
	c1.x = (m_vertices[t.ver1].x+m_vertices[t.ver0].x)/2;
	c1.y = (m_vertices[t.ver1].y+m_vertices[t.ver0].y)/2;
	c1.z = (m_vertices[t.ver1].z+m_vertices[t.ver0].z)/2;
	c2.x = (m_vertices[t.ver1].x+m_vertices[t.ver0].x)/2;
	c2.y = (m_vertices[t.ver1].y+m_vertices[t.ver0].y)/2;
	c2.z = (m_vertices[t.ver1].z+m_vertices[t.ver0].z)/2;
	v1.x = m_vertices[t.ver1].x-m_vertices[t.ver0].x;
	v1.y = m_vertices[t.ver1].y-m_vertices[t.ver0].y;
	v1.z = m_vertices[t.ver1].z-m_vertices[t.ver0].z;
	v2.x = m_vertices[t.ver2].x-m_vertices[t.ver0].x;
	v2.y = m_vertices[t.ver2].y-m_vertices[t.ver0].y;
	v2.z = m_vertices[t.ver2].z-m_vertices[t.ver0].z;
	v3 = CrossProduct(v1,v2);//the normal vector of this triangle,also of its underlying plane.
	v1 = CrossProduct(v1,v3);
	v2 = CrossProduct(v2,v3);
	NormalizeVector3D(&v1);
	NormalizeVector3D(&v2);
	v3 = v2;
	//so dv=c1+u*v1;dv=c2+v*v2; find (u,v).
	//so dv*v1=c1*v1+0;dv*v1=c2*v1+v*v2*v1;
	//   0 = (c2-c1)*v1+v v2*v1;
	v4.x = c2.x-c1.x;
	v4.y = c2.y-c1.y;
	v4.z = c2.z-c1.z;
    v2 = CrossProduct(v1,v2);
	v1 = CrossProduct(v4,v1);
	if(v2.x!=0)
	    asp = v1.x/v2.x;
	else if(v2.y!=0)
		asp = v1.y/v2.y;
	else
		asp = v1.z/v2.z;
    dv.x = c2.x+asp*v3.x;
	dv.y = c2.y+asp*v3.y;
	dv.z = c2.z+asp*v3.z;
	return dv;
}

double TriangularMesh::Norm(VECTOR3D v)
{
	return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

double TriangularMesh::NormSquare(VECTOR3D v)
{
	return v.x*v.x+v.y*v.y+v.z*v.z;
}

double TriangularMesh::Norm(VECTOR3D v,double alpha)
{
	return sqrt(v.x*v.x+v.y*v.y+v.z*v.z+alpha*alpha);
}

void TriangularMesh::BuildCCDualVerticesForTriangles()
{
	m_trianglesBC.clear();
	int i;
	VERTEX3D dv;
	for(i=0;i<m_trinum;i++)
	{
		dv = GetCCDual(m_triangles[i]);
		m_trianglesBC.push_back(dv);
	}
}

VERTEX3D TriangularMesh::GetBaryCenter(EDGE e)
{
	EDGEWITHCOOR ewc;
	ewc.ver0 = GetVertex3DFromPoint3D(m_vertices[e.ver0]);
	ewc.ver1 = GetVertex3DFromPoint3D(m_vertices[e.ver1]);
	return GetBaryCenter(ewc);
}

VERTEX3D TriangularMesh::GetBaryCenter(EDGEWITHCOOR ewc)
{
	VERTEX3D bc;
	bc.x = (ewc.ver0.x+ewc.ver1.x)/2;
	bc.y = (ewc.ver0.y+ewc.ver1.y)/2;
	bc.z = (ewc.ver0.z+ewc.ver1.z)/2;
	return bc;
}

VERTEX3D TriangularMesh::GetVertex3DFromPoint3D(POINT3d p)
{
	VERTEX3D v;
	v.x = p.x;
	v.y = p.y;
	v.z = p.z;
	return v;
}

COLORRGB TriangularMesh::GetColorFromPoint3D(POINT3d p)
{
	COLORRGB c;
	c.r = p.r;
	c.g = p.g;
	c.b = p.b;
	return c;
}

VERTEX3D TriangularMesh::GetBaryCenter(TRIANGLE t)
{
	TRIANGLEWITHCOOR twc;
	twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
	twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
	twc.ver2 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
	return GetBaryCenter(twc);
}

VERTEX3D TriangularMesh::GetBaryCenter(TRIANGLEWITHCOOR twc)
{
	VERTEX3D bc;
	bc.x = (twc.ver0.x+twc.ver1.x+twc.ver2.x)/3;
	bc.y = (twc.ver0.y+twc.ver1.y+twc.ver2.y)/3;
	bc.z = (twc.ver0.z+twc.ver1.z+twc.ver2.z)/3;
	return bc;
}

void TriangularMesh::BuildBaryCentersForTriangles()
{
	m_trianglesBC.clear();
	int i;
	VERTEX3D bc;
	for(i=0;i<m_trinum;i++)
	{
		bc = GetBaryCenter(m_triangles[i]);
		m_trianglesBC.push_back(bc);
	}
}

void TriangularMesh::FindPPIBGradient(int ver, TRIANGLE t, VECTOR3D *g, double *h)
{//to compute the gradient of the primal-primal interpolation function.
	double alpha;  double epsilon = 1.0e-12;
	VECTOR3D v1,v2;
	VERTEX3D pf;
    TRIANGLEWITHCOOR twc;
	twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
	twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
	twc.ver2 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
	
	if(ver==t.ver0)
	{
		v1.x = twc.ver1.x-twc.ver2.x;
		v1.y = twc.ver1.y-twc.ver2.y;
		v1.z = twc.ver1.z-twc.ver2.z;
		v2.x = twc.ver0.x-twc.ver2.x;
		v2.y = twc.ver0.y-twc.ver2.y;
		v2.z = twc.ver0.z-twc.ver2.z;
		alpha = DotProduct(v2,v1)/(DotProduct(v1,v1)+epsilon);
        //pf = alpha*twc.ver1+(1-alpha)*twc.ver2;
		pf.x = alpha*twc.ver1.x+(1-alpha)*twc.ver2.x;
		pf.y = alpha*twc.ver1.y+(1-alpha)*twc.ver2.y;
		pf.z = alpha*twc.ver1.z+(1-alpha)*twc.ver2.z;
		(*g).x = twc.ver0.x-pf.x;
		(*g).y = twc.ver0.y-pf.y;
		(*g).z = twc.ver0.z-pf.z;
		*h = sqrt(DotProduct(*g,*g))+epsilon;//this is the height corresponding to the edge opposite ver.
		(*g).x/=*h;
		(*g).y/=*h;
		(*g).z/=*h;
		(*g).x/=*h;
		(*g).y/=*h;
		(*g).z/=*h;
	}
	else if(ver==t.ver1)
	{
		v1.x = twc.ver2.x-twc.ver0.x;
		v1.y = twc.ver2.y-twc.ver0.y;
		v1.z = twc.ver2.z-twc.ver0.z;
		v2.x = twc.ver1.x-twc.ver0.x;
		v2.y = twc.ver1.y-twc.ver0.y;
		v2.z = twc.ver1.z-twc.ver0.z;
		alpha = DotProduct(v2,v1)/(DotProduct(v1,v1)+epsilon);
        pf.x = alpha*twc.ver2.x+(1-alpha)*twc.ver0.x;
		pf.y = alpha*twc.ver2.y+(1-alpha)*twc.ver0.y;
		pf.z = alpha*twc.ver2.z+(1-alpha)*twc.ver0.z;
		(*g).x = twc.ver1.x-pf.x;
		(*g).y = twc.ver1.y-pf.y;
		(*g).z = twc.ver1.z-pf.z;
		*h = sqrt(DotProduct(*g,*g))+epsilon;//this is the height corresponding to the edge opposite ver.
		(*g).x/=*h;
		(*g).y/=*h;
		(*g).z/=*h;
		(*g).x/=*h;
		(*g).y/=*h;
		(*g).z/=*h;
	}
	else
	{
		v1.x = twc.ver0.x-twc.ver1.x;
		v1.y = twc.ver0.y-twc.ver1.y;
		v1.z = twc.ver0.z-twc.ver1.z;
		v2.x = twc.ver2.x-twc.ver1.x;
		v2.y = twc.ver2.y-twc.ver1.y;
		v2.z = twc.ver2.z-twc.ver1.z;
		alpha = DotProduct(v2,v1)/(DotProduct(v1,v1)+epsilon);
        pf.x = alpha*twc.ver0.x+(1-alpha)*twc.ver1.x;
		pf.y = alpha*twc.ver0.y+(1-alpha)*twc.ver1.y;
		pf.z = alpha*twc.ver0.z+(1-alpha)*twc.ver1.z;
		(*g).x = twc.ver2.x-pf.x;
		(*g).y = twc.ver2.y-pf.y;
		(*g).z = twc.ver2.z-pf.z;
		*h = sqrt(DotProduct(*g,*g))+epsilon;//this is the height corresponding to the edge opposite ver.
		(*g).x/=*h;
		(*g).y/=*h;
		(*g).z/=*h;
		(*g).x/=*h;
		(*g).y/=*h;
		(*g).z/=*h;
	}
}

void TriangularMesh::FindPPIGradient(TRIANGLE t, VECTOR3D *g)
{
	double h;
	VECTOR3D g0,g1,g2;
	FindPPIBGradient(t.ver0,t,&g0,&h);
	FindPPIBGradient(t.ver1,t,&g1,&h);
	FindPPIBGradient(t.ver2,t,&g2,&h);
	(*g).x = m_vertices[t.ver0].r*g0.x+m_vertices[t.ver1].r*g1.x+m_vertices[t.ver2].r*g2.x;
	(*g).y = m_vertices[t.ver0].r*g0.y+m_vertices[t.ver1].r*g1.y+m_vertices[t.ver2].r*g2.y;
	(*g).z = m_vertices[t.ver0].r*g0.z+m_vertices[t.ver1].r*g1.z+m_vertices[t.ver2].r*g2.z;
}

void TriangularMesh::BuildTPPIBG()
{
	double h;
	m_TPPIBG.clear();
	TRIANGLEPPIBGRADIENT tppibg;
	int i;
	for(i=0;i<m_trinum;i++)
	{
		FindPPIBGradient(m_triangles[i].ver0,m_triangles[i],&(tppibg.v0),&h);
		FindPPIBGradient(m_triangles[i].ver1,m_triangles[i],&(tppibg.v1),&h);
		FindPPIBGradient(m_triangles[i].ver2,m_triangles[i],&(tppibg.v2),&h);
		m_TPPIBG.push_back(tppibg);
	}
}

VECTOR3D TriangularMesh::FindPPIGradient(int it)
{//from the index of the triangle to compute the ppi gradient of that.
	VECTOR3D g;
	g.x = m_vertices[m_triangles[it].ver0].r*m_TPPIBG[it].v0.x+m_vertices[m_triangles[it].ver1].r*m_TPPIBG[it].v1.x+
		m_vertices[m_triangles[it].ver2].r*m_TPPIBG[it].v2.x;
	g.y = m_vertices[m_triangles[it].ver0].r*m_TPPIBG[it].v0.y+m_vertices[m_triangles[it].ver1].r*m_TPPIBG[it].v1.y+
		m_vertices[m_triangles[it].ver2].r*m_TPPIBG[it].v2.y;
	g.z = m_vertices[m_triangles[it].ver0].r*m_TPPIBG[it].v0.z+m_vertices[m_triangles[it].ver1].r*m_TPPIBG[it].v1.z+
		m_vertices[m_triangles[it].ver2].r*m_TPPIBG[it].v2.z;
	return g;
}

void TriangularMesh::WriteTPPIBG()
{
	ofstream ofp("horse_TPPIBG.txt");
    //ofstream ofp("venushead_TPPIBG.txt");
	int i;
	ofp<<m_trinum<<"\n";
	for(i=0;i<m_trinum;i++)
	{
		ofp<<m_TPPIBG[i].v0.x<<"   "<<m_TPPIBG[i].v0.y<<"   "<<m_TPPIBG[i].v0.z<<"\n";
		ofp<<m_TPPIBG[i].v1.x<<"   "<<m_TPPIBG[i].v1.y<<"   "<<m_TPPIBG[i].v1.z<<"\n";
		ofp<<m_TPPIBG[i].v2.x<<"   "<<m_TPPIBG[i].v2.y<<"   "<<m_TPPIBG[i].v2.z<<"\n";
	}
	ofp.close();
}

VERTEX3D TriangularMesh::GetBaryCenter(VERTEX3D v)
{
	return v;
}

VECTOR3D TriangularMesh::GetBCDNormal(VERTEX3D v, EDGEWITHCOOR ewc, TRIANGLEWITHCOOR twc)
{
//to compute the normal vector of the boudary of the barycentric dual mesh at (v,ewc,twc).
	VECTOR3D n,v1,v2;
	float alpha;
    VERTEX3D bv,bewc,btwc;
	bv = GetBaryCenter(v);
	bewc = GetBaryCenter(ewc);
	btwc = GetBaryCenter(twc);
	v1 = GetVector3Dfrom2Vertices(btwc,bewc);
	v2 = GetVector3Dfrom2Vertices(btwc,bv);
	alpha = DotProduct(v2,v1)/DotProduct(v1,v1);
	n.x = alpha*v1.x-v2.x;
	n.y = alpha*v1.y-v2.y;
	n.z = alpha*v1.z-v2.z;
	NormalizeVector3D(&n);//normalization.
	return n;
}

VECTOR3D TriangularMesh::GetBCDNormal(int v,EDGE e,TRIANGLE t)
{
	VERTEX3D vv;
	EDGEWITHCOOR ewc;
	TRIANGLEWITHCOOR twc;
	vv = GetVertex3DFromPoint3D(m_vertices[v]);
	ewc.ver0 = GetVertex3DFromPoint3D(m_vertices[e.ver0]);
	ewc.ver1 = GetVertex3DFromPoint3D(m_vertices[e.ver1]);
	twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
	twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
	twc.ver2 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
	return GetBCDNormal(vv,ewc,twc);
}

VECTOR3D TriangularMesh::GetVector3Dfrom2Vertices(VERTEX3D begin, VERTEX3D end)
{
	VECTOR3D v;
	v.x = end.x-begin.x;
	v.y = end.y-begin.y;
	v.z = end.z-begin.z;
	return v;
}

void TriangularMesh::LBDiffusion()
{
	int i,j,time;
	double ts;//time step.
	double div;
	EDGE e;
	TRIANGLE t;
	//ts = 0.0005; //this ts is for venushead.
	ts = 0.0005; //this ts is for horse.
	for(time=0;time<m_nTime;time++)
	{
		for(i=0;i<m_vnum;i++)
		{
			div = 0;
			if(!m_vertices[i].IsOnBoundary)
			{//handle the inner points of the mesh.
				e.ver0 = i;
				for(j=0;j<m_onedisks[i].size();j++)
				{
					t = m_triangles[m_onedisks[i][j]];
					if(t.ver0==i)
					{
						e.ver1 = t.ver1;
						div+= DotProduct(FindPPIGradient(m_onedisks[i][j]),GetBCDNormal(i,e,t))*Norm(
							GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
						e.ver1 = t.ver2;
						div+= DotProduct(FindPPIGradient(m_onedisks[i][j]),GetBCDNormal(i,e,t))*Norm(
							GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
					}
					else if(t.ver1==i)
					{
						e.ver1 = t.ver0;
						div+= DotProduct(FindPPIGradient(m_onedisks[i][j]),GetBCDNormal(i,e,t))*Norm(
							GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
						e.ver1 = t.ver2;
						div+= DotProduct(FindPPIGradient(m_onedisks[i][j]),GetBCDNormal(i,e,t))*Norm(
							GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
					}
					else
					{
						e.ver1 = t.ver1;
						div+= DotProduct(FindPPIGradient(m_onedisks[i][j]),GetBCDNormal(i,e,t))*Norm(
							GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
						e.ver1 = t.ver0;
						div+= DotProduct(FindPPIGradient(m_onedisks[i][j]),GetBCDNormal(i,e,t))*Norm(
							GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
					}
				}
				m_tempRGB[i].r = m_vertices[i].r+ts*div;
			}
			else
			{
			}
		}
		for(i=0;i<m_vnum;i++)
		{
			m_vertices[i].r = m_tempRGB[i].r;
			m_vertices[i].g = m_vertices[i].r;
			m_vertices[i].b = m_vertices[i].r;
		}
	}
}

void TriangularMesh::BuildONEDISKS()
{
//与buildonedisks不同的是，这里构造的ONEDISKS已经带有lb算子中需要的顶点组合系数。
//此构建过程基于已经构建好的onedisks数据。
	int i,j;
	DISKTRIANGLE dt;
	DISK d;
	m_ONEDISKS.clear();
    EDGE e;
	TRIANGLE t;
	for(i=0;i<m_onedisks.size();i++)
	{
		d.clear();
		e.ver0 = i;
		for(j=0;j<m_onedisks[i].size();j++)
		{
			t = m_triangles[m_onedisks[i][j]];
			dt.itri = m_onedisks[i][j];
			dt.ver0co = 0;
			dt.ver1co = 0;
			dt.ver2co = 0;
			if(t.ver0==i)
			{
				e.ver1 = t.ver1;
				dt.ver1co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver2co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				e.ver1 = t.ver2;
				dt.ver1co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver2co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver0co = -dt.ver1co-dt.ver2co;
			}
			else if(t.ver1==i)
			{
				e.ver1 = t.ver0;
				dt.ver0co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver2co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				e.ver1 = t.ver2;
				dt.ver0co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver2co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver1co = -dt.ver0co-dt.ver2co;
			}
			else
			{
				e.ver1 = t.ver1;
				dt.ver1co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver0co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				e.ver1 = t.ver0;
				dt.ver1co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver0co+= DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dt.ver2co = -dt.ver1co-dt.ver0co;
			}
			d.push_back(dt);
		}
		m_ONEDISKS.push_back(d);
	}
/*
	ofstream ofp2("orderedonedisks2.txt");
	ofp2<<m_ONEDISKS.size()<<"\n";
	for(i=0;i<m_ONEDISKS.size();i++)
	{
		ofp2<<i<<"\n";
		for(j=0;j<m_ONEDISKS[i].size();j++)
		{
			t = m_triangles[m_ONEDISKS[i][j].itri];
		    ofp2<<"   "<<t.ver0<<":"<<m_ONEDISKS[i][j].ver0co<<";  "
				<<"   "<<t.ver1<<":"<<m_ONEDISKS[i][j].ver1co<<";  "
				<<"   "<<t.ver2<<":"<<m_ONEDISKS[i][j].ver2co<<"\n";
		}
		ofp2<<"\n";
	}
	ofp2.close();
*/
}

void TriangularMesh::ComputeBCDualArea()
{
//to compute the area of the bc dual of each vertex of the mesh.
	int i,j;
	double area;
	VECTOR3D v1,v2;
	EDGE e;
	TRIANGLE t;
	for(i=0;i<m_vnum;i++)
	{
		area = 0;
		e.ver0 = i;
		for(j=0;j<m_ONEDISKS[i].size();j++)
		{
			t = m_triangles[m_ONEDISKS[i][j].itri];
			if(t.ver0==i)
			{
				e.ver1 = t.ver1;
			    v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(e));
				v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(t));
				area+= Norm(CrossProduct(v1,v2))/2.0;
				e.ver1 = t.ver2;
				v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(e));
				v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(t));
				area+= Norm(CrossProduct(v1,v2))/2.0;
			}
			else if(t.ver1==i)
			{
				e.ver1 = t.ver0;
			    v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(e));
				v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(t));
				area+= Norm(CrossProduct(v1,v2))/2.0;
				e.ver1 = t.ver2;
				v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(e));
				v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(t));
				area+= Norm(CrossProduct(v1,v2))/2.0;
			}
			else
			{
				e.ver1 = t.ver1;
			    v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(e));
				v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(t));
				area+= Norm(CrossProduct(v1,v2))/2.0;
				e.ver1 = t.ver0;
				v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(e));
				v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[i]),GetBaryCenter(t));
				area+= Norm(CrossProduct(v1,v2))/2.0;
			}
		}
		m_vertices[i].BCDArea = area;
	}
/*
	ofstream ofp("venushead_dualareas.txt");
	ofp<<m_vnum<<"\n";
	double max,min;
	max = m_vertices[0].BCDArea;
	min = m_vertices[0].BCDArea;
	for(i=0;i<m_vnum;i++)
	{
		ofp<<m_vertices[i].BCDArea<<"\n";
		if(max<m_vertices[i].BCDArea)
			max = m_vertices[i].BCDArea;
		if(min>m_vertices[i].BCDArea)
			min = m_vertices[i].BCDArea;
	}
	ofp<<"max: "<<max<<"  min: "<<min;
	ofp.close();
*/
}

void TriangularMesh::LBDiffusion2()
{
	int i,j,time;
	double ts;//time step.
	double div;
	TRIANGLE t;
	ts = 0.0000000001; //this ts is for venushead.
	//ts = 0.00000005; //this ts is for horse.
	for(time=0;time<m_nTime;time++)
	{
		for(i=0;i<m_vnum;i++)
		{
			div = 0;
			if(!m_vertices[i].IsOnBoundary)
			{//handle the inner points of the mesh.
				for(j=0;j<m_ONEDISKS[i].size();j++)
				{
					t = m_triangles[m_ONEDISKS[i][j].itri];
					div+= m_vertices[t.ver0].r*m_ONEDISKS[i][j].ver0co+
						m_vertices[t.ver1].r*m_ONEDISKS[i][j].ver1co+
						m_vertices[t.ver2].r*m_ONEDISKS[i][j].ver2co;
				}	
				div/= m_vertices[i].BCDArea;
				m_tempRGB[i].r = m_vertices[i].r+ts*div;
			}
			else
			{
			}
		}
		for(i=0;i<m_vnum;i++)
		{
			m_vertices[i].r = m_tempRGB[i].r;
		}
	}
	for(i=0;i<m_vnum;i++)
	{
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}
}

void TriangularMesh::BuildOnerings()
{//还未考虑对偶网格的控制体面积。修正在后面的对角阵中。
	m_ONERINGS.clear();
	int i,j,k,ods;
	ONERINGVERTEX orv;
	ONERING or;
	DISKTRIANGLE dt,dtn;
	TRIANGLE t,tn;
	for(i=0;i<m_ONEDISKS.size();i++)
	{
		or.clear();
		ods = m_ONEDISKS[i].size();
		for(j=0;j<ods;j++)
		{//将onedisk中的所有三角形的非i 顶点都唯一加入到i 的onering.
			dt = m_ONEDISKS[i][j];
			t = m_triangles[dt.itri];
			if(t.ver0!=i)
			{
				orv.iver = t.ver0;
				orv.co = 0;
				if(!HasOrvInOnering(or,orv))
					or.push_back(orv);
			}
			if(t.ver1!=i)
			{
				orv.iver = t.ver1;
				orv.co = 0;
				if(!HasOrvInOnering(or,orv))
					or.push_back(orv);
			}
			if(t.ver2!=i)
			{
				orv.iver = t.ver2;
				orv.co = 0;
				if(!HasOrvInOnering(or,orv))
					or.push_back(orv);
			}
		}
		for(j=0;j<ods;j++)
		{//进一步确定每个顶点的lb 算子系数。
			dt = m_ONEDISKS[i][j];//第i 个节点的ONEDISKS的第j 个DISK三角形。
			t = m_triangles[dt.itri];//对应的三角片。
			for(k=0;k<or.size();k++)
			{
				if(or[k].iver==t.ver0)
					or[k].co+= dt.ver0co;
				if(or[k].iver==t.ver1)
                    or[k].co+= dt.ver1co;
				if(or[k].iver==t.ver2)
                    or[k].co+= dt.ver2co;
			}
		}
		m_ONERINGS.push_back(or);
	}
/*	
	ofstream ofp2("onerings.txt");
	ofp2<<m_ONERINGS.size()<<"\n";
	for(i=0;i<m_ONERINGS.size();i++)
	{
		ofp2<<i<<":\n";
		for(j=0;j<m_ONERINGS[i].size();j++)
		{
		    ofp2<<"   "<<m_ONERINGS[i][j].iver<<":"<<m_ONERINGS[i][j].co<<"\n";
		}
		ofp2<<"\n";
	}
	ofp2.close();
*/	
}

bool TriangularMesh::HasOrvInOnering(ONERING or, ONERINGVERTEX orv)
{
	int i;
	for(i=0;i<or.size();i++)
	{
		if(or[i].iver==orv.iver)
			break;
	}
	if(i==or.size())
		return false;
	return true;
}

void TriangularMesh::nr_dsprsax(unsigned long n, double x[], double b[], double tstep)
{
//参考 nrutil 的实现。
//用x 右乘该mesh 的lb算子系数矩阵，将结果保存于b 中。
//参数tstep 用于从lb算子系数矩阵确定最终的隐格式系数矩阵。
//这里的实现为l2 denoising算法；
//修改这里的程序，可得到各种其它算法。
//不同的模型算法对应的系数矩阵A不一样，因此此函数也要作相应更改。
	switch(m_ProcessingMethod)
	{
	case 1:
		nr_dsprsax_L2Diff(n, x, b, tstep);
		break;
	case 2:
		nr_dsprsax_TVDiff(n, x, b, tstep);
		break;
	case 3:
		nr_dsprsax_L2Denois(n, x, b, tstep);
		break;
	case 4:
		nr_dsprsax_TVDenois(n, x, b, tstep);
		break;
	case 5:
		nr_dsprsax_L2Inpaint(n, x, b, tstep);
		break;
	case 6:
		nr_dsprsax_TVInpaint(n, x, b, tstep);
		break;
	case 7:
		nr_dsprsax_L2RDTextur(n, x, b, tstep);
		break;
	case 8:
		nr_dsprsax_TVRDTextur(n, x, b, tstep);
		break;
	case 9:
		nr_dsprsax_L2DiffDirectionalData(n, x, b, tstep);
		break;
	case 10:
		nr_dsprsax_TVDiffDirectionalData(n, x, b, tstep);
		break;
	case 11:
		nr_dsprsax_AnisotropicRDTextur(n, x, b, tstep);
		break;
    case 12:
		nr_dsprsax_GeodesicCurvatureFlow(n, x, b, tstep);
		break;
	case 13:
		nr_dsprsax_GeodesicWeightedCurvatureFlow(n, x, b, tstep);
		break;
    case 14:
		// the same as the curve evolution under geodesic curvature flow.
		nr_dsprsax_GeodesicCurvatureFlow(n, x, b, tstep);
		break;
	case 16:
		nr_dsprsax_AnisotropicDAntialiasing(n, x, b, tstep);
		break;
	case 17:
		nr_dsprsax_NonlinearDiff(n, x, b, tstep);
		break;
	default://L2 diffusion as the default.
		nr_dsprsax_L2Diff(n, x, b, tstep);
		break;
	}
}

void TriangularMesh::nr_dsprsax_L2Diff(unsigned long n,double x[],double b[],double tstep)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]+= (-tstep)*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= ((1+0)*m_vertices[i].BCDArea-tstep*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::nr_dsprsax_L2DiffDirectionalData(unsigned long n,double x[],double b[],double tstep)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]+= (-tstep)*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= ((1+0)*m_vertices[i].BCDArea-tstep*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::nr_dsprsax_L2Denois(unsigned long n,double x[],double b[],double tstep)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]+= (-tstep)*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= ((1+m_Lambda*tstep)*m_vertices[i].BCDArea-tstep*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::nr_dsprsax_L2Inpaint(unsigned long n,double x[],double b[],double tstep)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]+= (-tstep)*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= ((1+m_Lambda*m_vertices[i].mask*tstep)*m_vertices[i].BCDArea-tstep*diag)*x[i]; //对角项,inpainting mask。
	}
}

void TriangularMesh::nr_dsprsax_TVDiff(unsigned long n,double x[],double b[],double tstep)
{
/*  for TV, it's a nonlinear case.
 *  the vector field is not grad(f) yet,but grad(f)/|grad(f)|.
 *  to get a linear system, semi-implicit scheme is used here,
 *  that is, we compute the 1/|grad(f)| by the values at the current
 *  time.
 */
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			b[i]+= q;
		}
		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_NonlinearDiff(unsigned long n,double x[],double b[],double tstep)
{
/*  A nonlinear case.
 *  the vector field is not grad(f) yet,but g(|grad(f)|)grad(f).
 *  to get a linear system, semi-implicit scheme is used here,
 *  that is, we compute g by the values at the current
 *  time.
 */
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			q*= gHeatCoeff(Norm(m_triangles[dt.itri].grad));
			b[i]+= q;
		}
		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_GeodesicCurvatureFlow(unsigned long n,double x[],double b[],double tstep)
{
/*  for geodesic curvature flow, a semi-implicit scheme is used.
 *  that is, for the nonlinear part we use the data at the current time step.
 *  a tip: the absolute gradient outside of div is approximated by average.
 */
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	double abgrad;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		abgrad = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			b[i]+= q;
			abgrad+= Norm(m_triangles[dt.itri].grad)*m_trianglesArea[dt.itri]/3;
		}
		abgrad/= m_vertices[i].BCDArea;
		b[i]*= abgrad;

		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_GeodesicWeightedCurvatureFlow(unsigned long n,double x[],double b[],double tstep)
{
/*  for geodesic curvature flow, a semi-implicit scheme is used.
 *  that is, for the nonlinear part we use the data at the current time step.
 *  a tip: the absolute gradient outside of div is approximated by average.
 */
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	double abgrad;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		abgrad = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			q*= m_triangles[dt.itri].edgeindicator;// added for edge indicator.
			b[i]+= q;
			abgrad+= Norm(m_triangles[dt.itri].grad)*m_trianglesArea[dt.itri]/3;
		}
		abgrad/= m_vertices[i].BCDArea;
		b[i]*= abgrad;

		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_TVDiffDirectionalData(unsigned long n,double x[],double b[],double tstep)
{
/*  for TV, it's a nonlinear case.
 */
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	VECTOR3D gd;//
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			gd.x = Norm(m_triangles[dt.itri].grad);
			gd.y = Norm(m_triangles[dt.itri].grad1);
			gd.z = Norm(m_triangles[dt.itri].grad2);//gradient of directional data.
			q/= Norm(gd,m_alpha);
			b[i]+= q;
		}
		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_TVDenois(unsigned long n,double x[],double b[],double tstep)
{
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			b[i]+= q;
		}
		b[i]+= (1+m_Lambda*tstep)*m_vertices[i].BCDArea*x[i];//
	}
}

void TriangularMesh::nr_dsprsax_TVInpaint(unsigned long n,double x[],double b[],double tstep)
{
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0];
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1];
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2];
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			b[i]+= q;
		}
		b[i]+= (1+m_Lambda*m_vertices[i].mask*tstep)*m_vertices[i].BCDArea*x[i];//
	}
}

void TriangularMesh::nr_dsprsax_L2RDTextur(unsigned long n,double x[],double b[],double tstep)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]+= m_RDTDiffRate*(-tstep)*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= ((1+0)*m_vertices[i].BCDArea-m_RDTDiffRate*tstep*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::nr_dsprsax_AnisotropicDAntialiasing(unsigned long n,double x[],double b[],double tstep)
{
	DISKANISOTROPIC da;
	DISKTRIANGLEANISOTROPIC dta;
	unsigned long i,j;
	double q;
	double s;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		da = m_ONEDISKSANISOTROPIC[i];//每个顶点对应的ONEDISKANISOTROPIC.
		for(j=0;j<da.size();j++)
		{//对该顶点的ONEDISKANISOTROPIC中的所有三角形循环。
			q = 0;
			dta = da[j];
            s = Norm(m_triangles[dta.itri].grad);
			q+= -tstep*(dta.d1ver0co*g1HeatCoeff(s)+dta.d2ver0co*g2HeatCoeff(s))*x[m_triangles[dta.itri].ver0];
			q+= -tstep*(dta.d1ver1co*g1HeatCoeff(s)+dta.d2ver1co*g2HeatCoeff(s))*x[m_triangles[dta.itri].ver1];
			q+= -tstep*(dta.d1ver2co*g1HeatCoeff(s)+dta.d2ver2co*g2HeatCoeff(s))*x[m_triangles[dta.itri].ver2];
			b[i]+= q;
		}
		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_AnisotropicRDTextur(unsigned long n,double x[],double b[],double tstep)
{
	DISKANISOTROPIC da;
	DISKTRIANGLEANISOTROPIC dta;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		da = m_ONEDISKSANISOTROPIC[i];//每个顶点对应的ONEDISKANISOTROPIC.
		for(j=0;j<da.size();j++)
		{//对该顶点的ONEDISKANISOTROPIC中的所有三角形循环。
			q = 0;
			dta = da[j];
			q+= -tstep*(dta.d1ver0co+dta.d2ver0co)*x[m_triangles[dta.itri].ver0];
			q+= -tstep*(dta.d1ver1co+dta.d2ver1co)*x[m_triangles[dta.itri].ver1];
			q+= -tstep*(dta.d1ver2co+dta.d2ver2co)*x[m_triangles[dta.itri].ver2];
			b[i]+= q*m_RDTDiffRate;
		}
		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprsax_TVRDTextur(unsigned long n,double x[],double b[],double tstep)
{
}

void TriangularMesh::nr_atimes(unsigned long n, double x[], double r[], double tstep, int itrnsp)
{
	if(itrnsp)
		nr_dsprstx(n,x,r,tstep);//转置乘。
	else
	    nr_dsprsax(n,x,r,tstep);
}

void TriangularMesh::nr_dsprstx(unsigned long n, double x[], double b[], double tstep)
{
    switch(m_ProcessingMethod)
	{
	case 1:
		nr_dsprsax_L2Diff(n, x, b, tstep);
		break;
	case 2:
		nr_dsprsax_TVDiff(n, x, b, tstep);
		break;
	case 3:
		nr_dsprsax_L2Denois(n, x, b, tstep);
		break;
	case 4:
		nr_dsprsax_TVDenois(n, x, b, tstep);
		break;
	case 5:
		nr_dsprsax_L2Inpaint(n, x, b, tstep);
		break;
	case 6:
		nr_dsprsax_TVInpaint(n, x, b, tstep);
		break;
	case 7:
		nr_dsprsax_L2RDTextur(n, x, b, tstep);
		break;
	case 8:
		nr_dsprsax_TVRDTextur(n, x, b, tstep);
		break;
	case 9:
		nr_dsprsax_L2DiffDirectionalData(n, x, b, tstep);
		break;
	case 10:
		nr_dsprsax_TVDiffDirectionalData(n, x, b, tstep);
		break;
	case 11:
		nr_dsprsax_AnisotropicRDTextur(n, x, b, tstep);
		break;
    case 12:
		// the coefficient matrix is not symmetric.
		nr_dsprstx_GeodesicCurvatureFlow(n, x, b, tstep);
		break;
	case 13:
		// the coefficient matrix is not symmetric.
		nr_dsprstx_WeightedGeodesicCurvatureFlow(n, x, b, tstep);
		break;
    case 14:
		// the same as the curve evolution under geodesic curvature flow.
		nr_dsprstx_GeodesicCurvatureFlow(n, x, b, tstep);
		break;
	case 16:
		nr_dsprsax_AnisotropicDAntialiasing(n, x, b, tstep);
		break;
	case 17:
		nr_dsprsax_NonlinearDiff(n, x, b, tstep);
		break;
	default:
		//L2 diffusion as the default.
		nr_dsprsax_L2Diff(n, x, b, tstep);
		break;
	}
}

void TriangularMesh::nr_dsprstx_GeodesicCurvatureFlow(unsigned long n, double x[], double b[], double tstep)
{
// ! the only difference between the coefficient matrix and its transpose is the averaged
// ! gradient term.
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			dt = d[j];
			q = 0;
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0]*m_vertices[m_triangles[dt.itri].ver0].averagedAbGradient;
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1]*m_vertices[m_triangles[dt.itri].ver1].averagedAbGradient;
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2]*m_vertices[m_triangles[dt.itri].ver2].averagedAbGradient;
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			b[i]+= q;
		}

		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_dsprstx_WeightedGeodesicCurvatureFlow(unsigned long n, double x[], double b[], double tstep)
{
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double q;
	for(i=0;i<m_vnum;i++)
	{
		b[i] = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			q = 0;
			dt = d[j];
			q+= -tstep*dt.ver0co*x[m_triangles[dt.itri].ver0]*m_vertices[m_triangles[dt.itri].ver0].averagedAbGradient;
			q+= -tstep*dt.ver1co*x[m_triangles[dt.itri].ver1]*m_vertices[m_triangles[dt.itri].ver1].averagedAbGradient;
			q+= -tstep*dt.ver2co*x[m_triangles[dt.itri].ver2]*m_vertices[m_triangles[dt.itri].ver2].averagedAbGradient;
			q/= Norm(m_triangles[dt.itri].grad,m_alpha);
			q*= m_triangles[dt.itri].edgeindicator;// added for edge indicator.
			b[i]+= q;
		}

		b[i]+= m_vertices[i].BCDArea*x[i];
	}
}

void TriangularMesh::nr_asolve(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
//由于采用对角阵作为preconditioner,所以参数itrnsp 暂时不用。
//不同的模型算法对应的系数矩阵A不一样，因此次函数也要做相应改变。
	switch(m_ProcessingMethod)
	{
	case 1:
		nr_asolve_L2Diff(n, b, x, tstep, itrnsp);
		break;
	case 2:
		nr_asolve_TVDiff(n, b, x, tstep, itrnsp);
		break;
	case 3:
		nr_asolve_L2Denois(n, b, x, tstep, itrnsp);
		break;
	case 4:
		nr_asolve_TVDenois(n, b, x, tstep, itrnsp);
		break;
	case 5:
		nr_asolve_L2Inpaint(n, b, x, tstep, itrnsp);
		break;
	case 6:
		nr_asolve_TVInpaint(n, b, x, tstep, itrnsp);
		break;
	case 7:
		nr_asolve_L2RDTextur(n, b, x, tstep, itrnsp);
		break;
	case 8:
		nr_asolve_TVRDTextur(n, b, x, tstep, itrnsp);
		break;
	case 9:
		nr_asolve_L2DiffDirectionalData(n, b, x, tstep, itrnsp);
		break;
	case 10:
		nr_asolve_TVDiffDirectionalData(n, b, x, tstep, itrnsp);
		break;
	case 11:
		nr_asolve_AnisotropicRDTextur(n, b, x, tstep, itrnsp);
		break;
    case 12:
		nr_asolve_GeodesicCurvatureFlow(n, b, x, tstep, itrnsp);
		break;
    case 13:
		nr_asolve_GeodesicWeightedCurvatureFlow(n, b, x, tstep, itrnsp);
		break;
    case 14:
		nr_asolve_GeodesicCurvatureFlow(n, b, x, tstep, itrnsp);
		break;
	case 16:
		nr_asolve_AnisotropicDAntialiasing(n, b, x, tstep, itrnsp);
		break;
	case 17:
		nr_asolve_NonlinearDiff(n, b, x, tstep, itrnsp);
		break;
	default://L2 diffusion as the default.
		nr_asolve_L2Diff(n, b, x, tstep, itrnsp);
		break;
	}
}

void TriangularMesh::nr_asolve_L2Diff(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = (1+0)*m_vertices[i].BCDArea-tstep*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_L2DiffDirectionalData(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = (1+0)*m_vertices[i].BCDArea-tstep*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_L2Denois(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = (1+m_Lambda*tstep)*m_vertices[i].BCDArea-tstep*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_L2Inpaint(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = (1+m_Lambda*m_vertices[i].mask*tstep)*m_vertices[i].BCDArea-tstep*diag;//inpainting mask.
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_TVDiff(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
//关键计算对角阵元素。
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else
				diag+= -tstep*dt.ver2co/Norm(m_triangles[dt.itri].grad,m_alpha);
		}
		diag+= m_vertices[i].BCDArea;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_NonlinearDiff(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
//关键计算对角阵元素。
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co*gHeatCoeff(Norm(m_triangles[dt.itri].grad));
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co*gHeatCoeff(Norm(m_triangles[dt.itri].grad));
			else
				diag+= -tstep*dt.ver2co*gHeatCoeff(Norm(m_triangles[dt.itri].grad));
		}
		diag+= m_vertices[i].BCDArea;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_GeodesicCurvatureFlow(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
//关键计算对角阵元素。
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	double abgrad;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		abgrad = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{//遍历该顶点的1-disk三角形。
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else
				diag+= -tstep*dt.ver2co/Norm(m_triangles[dt.itri].grad,m_alpha);
			abgrad+= Norm(m_triangles[dt.itri].grad)*m_trianglesArea[dt.itri]/3;
		}
		abgrad/= m_vertices[i].BCDArea;
		diag*= abgrad;
		diag+= m_vertices[i].BCDArea;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_GeodesicWeightedCurvatureFlow(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
//关键计算对角阵元素。
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	double abgrad;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		abgrad = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{//遍历该顶点的1-disk三角形。
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co/Norm(m_triangles[dt.itri].grad,m_alpha)*m_triangles[dt.itri].edgeindicator;
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co/Norm(m_triangles[dt.itri].grad,m_alpha)*m_triangles[dt.itri].edgeindicator;
			else
				diag+= -tstep*dt.ver2co/Norm(m_triangles[dt.itri].grad,m_alpha)*m_triangles[dt.itri].edgeindicator;
			abgrad+= Norm(m_triangles[dt.itri].grad)*m_trianglesArea[dt.itri]/3;
		}
		abgrad/= m_vertices[i].BCDArea;
		diag*= abgrad;
		diag+= m_vertices[i].BCDArea;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_TVDiffDirectionalData(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
//关键计算对角阵元素。
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	VECTOR3D gd;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{
			dt = d[j];
			gd.x = Norm(m_triangles[dt.itri].grad);
			gd.y = Norm(m_triangles[dt.itri].grad1);
			gd.z = Norm(m_triangles[dt.itri].grad2);
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co/Norm(gd,m_alpha);
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co/Norm(gd,m_alpha);
			else
				diag+= -tstep*dt.ver2co/Norm(gd,m_alpha);
		}
		diag+= m_vertices[i].BCDArea;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_TVDenois(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else
				diag+= -tstep*dt.ver2co/Norm(m_triangles[dt.itri].grad,m_alpha);
		}
		diag+= (1+m_Lambda*tstep)*m_vertices[i].BCDArea;//
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_TVInpaint(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	DISK d;
	DISKTRIANGLE dt;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		d = m_ONEDISKS[i];
		for(j=0;j<d.size();j++)
		{
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				diag+= -tstep*dt.ver0co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else if(i==m_triangles[dt.itri].ver1)
				diag+= -tstep*dt.ver1co/Norm(m_triangles[dt.itri].grad,m_alpha);
			else
				diag+= -tstep*dt.ver2co/Norm(m_triangles[dt.itri].grad,m_alpha);
		}
		diag+= (1+m_Lambda*m_vertices[i].mask*tstep)*m_vertices[i].BCDArea;//
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_L2RDTextur(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = (1+0)*m_vertices[i].BCDArea-m_RDTDiffRate*tstep*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_AnisotropicDAntialiasing(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	DISKANISOTROPIC da;
	DISKTRIANGLEANISOTROPIC dta;
	unsigned long i,j;
	double diag,s;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		da = m_ONEDISKSANISOTROPIC[i];//每个顶点对应的ONEDISKANISOTROPIC.
		for(j=0;j<da.size();j++)
		{//对该顶点的ONEDISKANISOTROPIC中的所有三角形循环。
			dta = da[j];
            s = Norm(m_triangles[dta.itri].grad);
			if(m_triangles[dta.itri].ver0==i)
				diag+= -tstep*(dta.d1ver0co*g1HeatCoeff(s)+dta.d2ver0co*g2HeatCoeff(s));
			else if(m_triangles[dta.itri].ver1==i)
				diag+= -tstep*(dta.d1ver1co*g1HeatCoeff(s)+dta.d2ver1co*g2HeatCoeff(s));
			else
				diag+= -tstep*(dta.d1ver2co*g1HeatCoeff(s)+dta.d2ver2co*g2HeatCoeff(s));
		}
		diag = m_vertices[i].BCDArea + diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_AnisotropicRDTextur(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
	DISKANISOTROPIC da;
	DISKTRIANGLEANISOTROPIC dta;
	unsigned long i,j;
	double diag;
	for(i=0;i<m_vnum;i++)
	{
		diag = 0;
		da = m_ONEDISKSANISOTROPIC[i];//每个顶点对应的ONEDISKANISOTROPIC.
		for(j=0;j<da.size();j++)
		{//对该顶点的ONEDISKANISOTROPIC中的所有三角形循环。
			dta = da[j];
			if(m_triangles[dta.itri].ver0==i)
				diag+= -tstep*(dta.d1ver0co+dta.d2ver0co);
			else if(m_triangles[dta.itri].ver1==i)
				diag+= -tstep*(dta.d1ver1co+dta.d2ver1co);
			else
				diag+= -tstep*(dta.d1ver2co+dta.d2ver2co);
		}
		diag*= m_RDTDiffRate;
		diag = m_vertices[i].BCDArea + diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::nr_asolve_TVRDTextur(unsigned long n, double b[], double x[], double tstep, int itrnsp)
{
}

double TriangularMesh::nr_snrm(unsigned long n, double sx[], int itol)
{
//compute one of 2 norms for a vector sx[1..n],as signaled by itol. used by nr_linbcg.
	unsigned long i,imax;
	double ans;
	if(itol<=3)
	{
		ans = 0.0;
		for(i=0;i<n;i++)
			ans+= sx[i]*sx[i];
		return sqrt(ans);
	}
	else
	{
		imax = 0;
		for(i=0;i<n;i++)
			if(fabs(sx[i])>fabs(sx[imax]))
				imax = i;
		ans = fabs(sx[imax]);
		return ans;
	}
}

void TriangularMesh::nr_linbcg(unsigned long n, double b[], double x[], double tstep, int itol, double tol, int itmax, int *iter, double *err)
{
/*Solves Ax=b for x[1..n],given b[1..n],by the iterative preconditional biconjugate gradient method.
  on input x[1..n] should be set to an initial guess of the solution or all zeros;
  itol is 1,2,3,4, specialfying which convergence test is applied;
  itmax is allowed maximum number of iterations; tol is the desired convergence tolerance;
  tstep is the time step of the implicit scheme of the PDE.
  on output x[1..n] is set to the improved solution; iter is the number of iterations actually taken;
  err is the estimated error.
  the matrix A is referenced only through the user-supplied routines nr_atimes, which computes
  the product of A or its transpose on a vector; and nr_asolve, which solves for precondition, 
  whose coefficient matrix is usually the diagonal part of A.
 */
	unsigned long j;
    double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zminrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p = new double[n];
	pp = new double[n];
	r = new double[n];
	rr = new double[n];
	z = new double[n];
	zz = new double[n];

	//calculate initial residual.
	*iter = 0;
	nr_atimes(n,x,r,tstep,0);
	for(j=0;j<n;j++)
	{
		r[j] = b[j]-r[j];
		rr[j] = r[j];
	}
	/* nr_atimes(n,r,rr,tstep,0) */
	if(itol==1)
	{
		bnrm = nr_snrm(n,b,itol);
		nr_asolve(n,r,z,tstep,0);
	}
	else if(itol==2)
	{
		nr_asolve(n,b,z,tstep,0);
		bnrm = nr_snrm(n,z,itol);
		nr_asolve(n,r,z,tstep,0);
	}
	else if(itol==3 || itol==4)
	{
		nr_asolve(n,b,z,tstep,0);
		bnrm = nr_snrm(n,z,itol);
		nr_asolve(n,r,z,tstep,0);
		znrm = nr_snrm(n,z,itol);
	}
	else return;//
    while(*iter<=itmax)
	{// main loop.
		++(*iter);
		nr_asolve(n,rr,zz,tstep,1);
		for(bknum=0.0,j=0;j<n;j++)
			bknum+= z[j]*rr[j];
		/*  calculate coefficient bk and direction vectors p and pp  */
		if(*iter==1)
		{
			for(j=0;j<n;j++)
			{
				p[j] = z[j];
				pp[j] = zz[j];
			}
		}
		else
		{
			bk = bknum/bkden;
			for(j=0;j<n;j++)
			{
				p[j] = bk*p[j]+z[j];
				pp[j] = bk*pp[j]+zz[j];
			}
		}
		bkden = bknum;//calculate new coefficients ak, new x, .
		nr_atimes(n,p,z,tstep,0);
		for(akden=0.0,j=0;j<n;j++)
			akden+= z[j]*pp[j];
		ak = bknum/akden;
		nr_atimes(n,pp,zz,tstep,1);
		for(j=0;j<n;j++)
		{
			x[j]+= ak*p[j];
			r[j]-= ak*z[j];
			rr[j]-= ak*zz[j];
		}
		nr_asolve(n,r,z,tstep,0);//
		if(itol==1)
			*err = nr_snrm(n,r,itol)/bnrm;
		else if(itol==2)
			*err = nr_snrm(n,z,itol)/bnrm;
		else if(itol==3 || itol==4)
		{
			zminrm = znrm;
			znrm = nr_snrm(n,z,itol);
			if(fabs(zminrm-znrm)>EPS*znrm)
			{
				dxnrm = fabs(ak)*nr_snrm(n,p,itol);
				*err = znrm/fabs(zminrm-znrm)*dxnrm;
			}
			else
			{
				*err = znrm/bnrm;
				continue;
			}
			xnrm = nr_snrm(n,x,itol);
			if(*err<=0.5*xnrm) *err/= xnrm;
			else
			{
				*err = znrm/bnrm;
				continue;
			}
		}
		if(*err<=tol)
			break;
	}
	delete r;
	delete rr;
	delete p;
	delete pp;
	delete z;
	delete zz;
}

void TriangularMesh::ImplicitL2Denoising()
{
//这里的实现应该考虑到系数矩阵的对称性。
	SetProcessingMethod(3);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*tstep*m_tempRGB[i].r);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
    for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*tstep*m_tempRGB[i].g);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*tstep*m_tempRGB[i].b);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitAnisotropicDAntialiasing()
{
//for anisotropic diffusion antialiasing.
    SetProcessingMethod(16);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

	delete u;
	delete uu;

//	WriteData();
}

void TriangularMesh::ImplicitL2RDTextur()
{
//generating patterns on meshes via reaction diffusion equations.
//here 2 chemicals are considered.
//semi-implicit scheme is used here.
    SetProcessingMethod(7);
    //SetRDTextureParam(0.0001,0.00000625,0.1,16,12);//for venushead.
	//SetRDTextureParam(0.005,0.00000625,0.14,16,12);//for bunny21_2.
	//SetRDTextureParam(0.00005,0.00000625,0.14,16,12);//for armadillo11.
    //SetRDTextureParam(0.005,0.000625,0.14,9,3);//for rabbit_param_nofilehead.
	//SetRDTextureParam(0.005,0.00000625,0.14,9,3);//for knot.
	SetRDTextureParam(0.00005,0.00000625,0.14,9,3);//for teapot.

	ImplicitL2RDTTuring_Turing();
//  ImplicitL2RDTTuring_Brusselator();	
	WriteData();
}

void TriangularMesh::ImplicitAnisotropicRDTextur()
{
	SetProcessingMethod(11);
	SetRDTextureParam(0.005,0.00000625,0.14,14,4);
	
	//Turing RD model ...
	ImplicitL2RDTTuring_Turing();
}

void TriangularMesh::SetProcessingMethod(unsigned short pm)
{//不同的处理方法得到的系数矩阵不一样，因此根据这个参数选择
 //不同的Ax函数调用以及A_{-1}函数调用。
 //也意味着对不同的方法要编写相应的这两个函数。
	m_ProcessingMethod = pm;
}

double TriangularMesh::RDTexturing_TurkF(double a, double b,double rm)
{
//the reaction term for chemical a.
//actually Turing's model.
	return m_RDTp.ar*(m_RDTp.gr+rm-a*b);
}

double TriangularMesh::RDTexturing_TurkG(double a, double b,double rm)
{
//the reaction term for chemical b.
//actually Turing's model.
	return m_RDTp.ar*(a*b-b-m_RDTp.decay-rm);
}

void TriangularMesh::ImplicitL2Diff_Contour(double smoothtime)
{//   ! 对contour数据的L2smoothing过程，只迭代一步。
    SetProcessingMethod(1);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = smoothtime;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].contour;
		uu[i] = u[i];
	}
	for(time=0;time<1;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].contour = uu[i];
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitL2Diff()
{
    SetProcessingMethod(1);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
/*
    //then, the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}
*/
    for(i=0;i<n;i++)
	{
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitL2DiffDirecData()
{
    SetProcessingMethod(9);

	unsigned long n,i,j;
	double *u1,*uu1,*u2,*uu2,*u3,*uu3;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;
	DISK d;
	DISKTRIANGLE dt;
    double gd;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u1 = new double[n];
	uu1 = new double[n];
	u2 = new double[n];
	uu2 = new double[n];
	u3 = new double[n];
	uu3 = new double[n];
	//initialization...
	
	for(i=0;i<n;i++)
	{
		u1[i] = m_vertices[i].normal_x;
		u2[i] = m_vertices[i].normal_y;
		u3[i] = m_vertices[i].normal_z;
		uu1[i] = u1[i];
		uu2[i] = u2[i];
		uu3[i] = u3[i];
	}
	/*
	for(i=0;i<n;i++)
	{
		u1[i] = m_vertices[i].pde2.x;
		u2[i] = m_vertices[i].pde2.y;
		u3[i] = m_vertices[i].pde2.z;
		uu1[i] = u1[i];
		uu2[i] = u2[i];
		uu3[i] = u3[i];
	}
	*/
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu1);
		GetPPIGradients1(n,uu2);
		GetPPIGradients2(n,uu3);
	    for(i=0;i<n;i++)
		{
			gd = 0.0f;
			d = m_ONEDISKS[i];
            for(j=0;j<d.size();j++)
			{
				dt = d[j];
				gd+= (POWER(Norm(m_triangles[dt.itri].grad))+
					POWER(Norm(m_triangles[dt.itri].grad1))+
					POWER(Norm(m_triangles[dt.itri].grad2)))*
					m_trianglesArea[dt.itri]/3;
			}
	    	u1[i] = m_vertices[i].BCDArea*uu1[i]+uu1[i]*tstep*gd;
			u2[i] = m_vertices[i].BCDArea*uu2[i]+uu2[i]*tstep*gd;
			u3[i] = m_vertices[i].BCDArea*uu3[i]+uu3[i]*tstep*gd;
		}
	    nr_linbcg(n,u1,uu1,tstep,itol,tol,itmax,&iter,&err);
		nr_linbcg(n,u2,uu2,tstep,itol,tol,itmax,&iter,&err);
		nr_linbcg(n,u3,uu3,tstep,itol,tol,itmax,&iter,&err);
		for(i=0;i<n;i++)
			NormalizeVector3D(uu1+i,uu2+i,uu3+i);
	}
	
	for(i=0;i<n;i++)
	{
		m_vertices[i].normal_x = uu1[i];
		m_vertices[i].normal_y = uu2[i];
		m_vertices[i].normal_z = uu3[i];
	}
	/*
	for(i=0;i<n;i++)
	{
		m_vertices[i].pde2.x = uu1[i];
		m_vertices[i].pde2.y = uu2[i];
		m_vertices[i].pde2.z = uu3[i];
	}

	ofstream ofp("pde2normaldp.txt");
	ofp<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
		ofp<<(m_vertices[i].pde2.x*m_vertices[i].normal_x+
		m_vertices[i].pde2.y*m_vertices[i].normal_y+
		m_vertices[i].pde2.z*m_vertices[i].normal_z)<<"\n";
	ofp.close();*/

	delete u1;
	delete uu1;
	delete u2;
	delete uu2;
	delete u3;
	delete uu3;
}

void TriangularMesh::ImplicitL2Inpaint()
{
	SetProcessingMethod(5);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*m_vertices[i].mask*tstep*m_tempRGB[i].r);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
    for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*m_vertices[i].mask*tstep*m_tempRGB[i].g);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*m_vertices[i].mask*tstep*m_tempRGB[i].b);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitTVDiff()
{
	SetProcessingMethod(2);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitNonlinearDiff()
{
	SetProcessingMethod(17);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitGeodesicCurvatureFlow()
{
// in this routine, we realize the geodesic curvature flow of a curve
// representing as a zero level set of a function.

//	Reinitialize_Contour(0.0005);//求解level set方程之前重新初始化曲线表示函数。

	SetProcessingMethod(12);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	// the contour montion.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].contour;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);//各三角形内的uu数据的梯度。
		ComputeAveragedVertexAbGradient();//各顶点处的由uu数据决定的平均梯度模（调用此函数前必须调用上面的函数）。
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].contour = uu[i];
	}

	m_WeightedCurveLength.push_back(ComputeWeightedCurveLength());// to compute the wcl.
	m_tSteps.push_back(m_tStep);// m_tStep does not change in this routine.
	m_WeightedCurveLength_Forobservation.push_back(ComputeWeightedCurveLength());
	m_tSteps_Forobservation.push_back(m_tStep);

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitGeodesicCurvatureFlow_AdaptiveTimeSteps()
{
// in this routine, we realize the geodesic curvature flow of a curve
// representing as a zero level set of a function.
// Also the time steps are chosen adaptively to avoid missing some closed geodesics.

//	Reinitialize_Contour(0.0005);//求解level set方程之前重新初始化曲线表示函数。

	SetProcessingMethod(12);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	// the contour montion.
	for(i=0;i<n;i++)
	{
		m_tempContour[i] = m_vertices[i].contour;
		u[i] = m_vertices[i].contour;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);//各三角形内的uu数据的梯度。
		ComputeAveragedVertexAbGradient();//各顶点处的由uu数据决定的平均梯度模（调用此函数前必须调用上面的函数）。
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].contour = uu[i];
	}

	m_WeightedCurveLength_Forobservation.push_back(ComputeWeightedCurveLength());
	m_WeightedCurveLength.push_back(ComputeWeightedCurveLength());
    
	ResetContourTimestep_3();
	m_tSteps_Forobservation.push_back(m_tStep);
	m_tSteps.push_back(m_tStep);

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitScaleSpacebyGCF()
{
// in this routine, we realize the scale space construction of images by geodesic curvature flow
// of the grey levels of those images.
	SetProcessingMethod(14);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	// the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
		ComputeAveragedVertexAbGradient();//
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
/*
	// the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
		ComputeAveragedVertexAbGradient();//
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}

	// the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
		ComputeAveragedVertexAbGradient();//
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}
*/
    for(i=0;i<n;i++)
	{
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitGeodesicWeightedCurvatureFlow()
{
// in this routine, we realize the weighted geodesic curvature flow with an edge indicator of a curve
// representing as a zero level set of a function.
	SetProcessingMethod(13);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;
	//CTime begin,end;
	//CTimeSpan ts;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//begin = CTime::GetCurrentTime();
	// the contour montion.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].contour;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
		ComputeAveragedVertexAbGradient();//
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].contour = uu[i];
	}

	//end = CTime::GetCurrentTime();
	//ts = end - begin;

	m_WeightedCurveLength.push_back(ComputeWeightedCurveLength());
	m_tSteps.push_back(m_tStep);
	m_WeightedCurveLength_Forobservation.push_back(ComputeWeightedCurveLength());
	m_tSteps_Forobservation.push_back(m_tStep);
	//m_CPUtime.push_back(ts.GetTotalSeconds());

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitGeodesicWeightedCurvatureFlow_AdaptiveTimeSteps()
{
// in this routine, we realize the weighted geodesic curvature flow with an edge indicator of a curve
// representing as a zero level set of a function.
	SetProcessingMethod(13);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	// the contour montion.
	for(i=0;i<n;i++)
	{
		m_tempContour[i] = m_vertices[i].contour;
		u[i] = m_vertices[i].contour;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
		ComputeAveragedVertexAbGradient();//
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+0);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].contour = uu[i];
	}

	m_WeightedCurveLength_Forobservation.push_back(ComputeWeightedCurveLength());
	m_WeightedCurveLength.push_back(ComputeWeightedCurveLength());
    
	ResetContourTimestep_3();
	m_tSteps_Forobservation.push_back(m_tStep);
	m_tSteps.push_back(m_tStep);

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitTVDiffDirecData()
{
	SetProcessingMethod(10);

	unsigned long n,i,j;
	double *u1,*uu1,*u2,*uu2,*u3,*uu3;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;
	DISK d;
	DISKTRIANGLE dt;
    double gd;
	VECTOR3D gdd;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u1 = new double[n];
	uu1 = new double[n];
	u2 = new double[n];
	uu2 = new double[n];
	u3 = new double[n];
	uu3 = new double[n];
	//initialization...
	for(i=0;i<n;i++)
	{
		u1[i] = m_vertices[i].normal_x;
		u2[i] = m_vertices[i].normal_y;
		u3[i] = m_vertices[i].normal_z;
		uu1[i] = u1[i];
		uu2[i] = u2[i];
		uu3[i] = u3[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu1);
		GetPPIGradients1(n,uu2);
		GetPPIGradients2(n,uu3);
	    for(i=0;i<n;i++)
		{
			gd = 0.0f;
			d = m_ONEDISKS[i];
            for(j=0;j<d.size();j++)
			{
				dt = d[j];
				gdd.x = Norm(m_triangles[dt.itri].grad);
				gdd.y = Norm(m_triangles[dt.itri].grad1);
				gdd.z = Norm(m_triangles[dt.itri].grad2);
				gd+= Norm(gdd)*m_trianglesArea[dt.itri]/3;
			}
	    	u1[i] = m_vertices[i].BCDArea*uu1[i]+uu1[i]*tstep*gd;
			u2[i] = m_vertices[i].BCDArea*uu2[i]+uu2[i]*tstep*gd;
			u3[i] = m_vertices[i].BCDArea*uu3[i]+uu3[i]*tstep*gd;
		}
	    nr_linbcg(n,u1,uu1,tstep,itol,tol,itmax,&iter,&err);
		nr_linbcg(n,u2,uu2,tstep,itol,tol,itmax,&iter,&err);
		nr_linbcg(n,u3,uu3,tstep,itol,tol,itmax,&iter,&err);
		for(i=0;i<n;i++)
			NormalizeVector3D(uu1+i,uu2+i,uu3+i);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].normal_x = uu1[i];
		m_vertices[i].normal_y = uu2[i];
		m_vertices[i].normal_z = uu3[i];
	}

	delete u1;
	delete uu1;
	delete u2;
	delete uu2;
	delete u3;
	delete uu3;
}

void TriangularMesh::ImplicitTVDenois()
{
	SetProcessingMethod(4);

    unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*tstep*m_tempRGB[i].r);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*tstep*m_tempRGB[i].g);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*tstep*m_tempRGB[i].b);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

    WriteData();

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitTVInpaint()
{
	SetProcessingMethod(6);

	unsigned long n,i;
	double *u,*uu;
	int itol,itmax,iter;
    double tol,err,tstep;
	int time;

	n = m_vnum;
	tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[n];
	uu = new double[n];
	//firstly, the r component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].r;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*m_vertices[i].mask*tstep*m_tempRGB[i].r);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = uu[i];
	}
    //then, the g component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].g;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*m_vertices[i].mask*tstep*m_tempRGB[i].g);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = uu[i];
	}
	//at last, the b component.
	for(i=0;i<n;i++)
	{
		u[i] = m_vertices[i].b;
		uu[i] = u[i];
	}
	for(time=0;time<m_nTime;time++)
	{
		GetPPIGradients(n,uu);
	    for(i=0;i<n;i++)
	    	u[i] = m_vertices[i].BCDArea*(uu[i]+m_Lambda*m_vertices[i].mask*tstep*m_tempRGB[i].b);  
	    nr_linbcg(n,u,uu,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].b = uu[i];
	}

	delete u;
	delete uu;
}

void TriangularMesh::ImplicitTVRDTextur()
{
	SetProcessingMethod(8);
}

void TriangularMesh::SetRDTextureParam(double Da, double Db, double ar, double gr, double decay)
{
	m_RDTp.Da = Da;
	m_RDTp.Db = Db;
	m_RDTp.ar = ar;
	m_RDTp.gr = gr;
	m_RDTp.decay = decay;
}

void TriangularMesh::WriteData()
{
	unsigned long i;
//	ofstream ofp("Mesh_onlyContour_.txt");
//	ofstream ofp("UnifiedMesh_.txt");
	ofstream ofp("ColoredMesh_.txt");
//	ofstream ofp("RDtexturedMesh_.txt");
/*	ofp<<"Da: "<<m_RDTp.Da<<";\n";
	ofp<<"Db: "<<m_RDTp.Db<<";\n";
	ofp<<"reaction rate: "<<m_RDTp.ar<<";\n";
	ofp<<"growth rate: "<<m_RDTp.gr<<";\n";
	ofp<<"decay rate:  "<<m_RDTp.decay<<";\n";
	*/
	ofp<<m_vnum<<"   "<<m_trinum<<"\n";
	for(i=0;i<m_vnum;i++)
	{
		ofp<<m_vertices[i].x<<"   "<<m_vertices[i].y<<"   "<<m_vertices[i].z<<"   "
			<<m_vertices[i].r<<"   "<<m_vertices[i].g<<"   "<<m_vertices[i].b<<"   \n";
//			<<m_vertices[i].mask<<"\n";// for inpainting.
//            <<m_vertices[i].contour<<"\n";// for edge detection.
	}
	for(i=0;i<m_trinum;i++)
	{
		ofp<<m_triangles[i].ver0<<"   "<<m_triangles[i].ver1<<"   "<<m_triangles[i].ver2<<"\n";
	}
	ofp.close();
}

void TriangularMesh::AddUniformNoise(float mag)
{
	unsigned long i;
	double rd;
	for(i=0;i<m_vnum;i++)
	{
		rd = mag*double(rand())/RAND_MAX-mag/2.0;
		m_vertices[i].r+= rd;
//		rd = mag*double(rand())/RAND_MAX-mag/2.0;
		m_vertices[i].g+= rd;
//		rd = mag*double(rand())/RAND_MAX-mag/2.0;
		m_vertices[i].b+= rd;
	}
}


double TriangularMesh::RDTexturing_BrusselatorF(double a, double b, double s, double alpha, double beta)
{
//Brusselator's reaction term.
    return s*(alpha-(1+beta)*a+a*a*b);
}

double TriangularMesh::RDTexturing_BrusselatorG(double a, double b, double s, double alpha, double beta)
{
//Brusselator's reaction term.
    return s*(beta*a-a*a*b);
}

void TriangularMesh::ImplicitL2RDTTuring_Turing()
{
//Turing's model with Turing's reaction term: only 2 chemicals.
	unsigned long n,i;
	int time;
	double *a,*aa,*b,*bb,*random;

	int itol,itmax,iter;
    double tol,err,tstep;

    tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	n = m_vnum;
	a = new double[n];
	aa = new double[n];
	b = new double[n];
	bb = new double[n];
	random = new double[n];
	float zmin,zmax;
	zmin = m_vertices[0].y;
	zmax = zmin;
	for(i=0;i<n;i++)
	{
		if(m_vertices[i].z<zmin) zmin=m_vertices[i].z;
		if(m_vertices[i].z>zmax) zmax=m_vertices[i].z;
	}
	for(i=0;i<n;i++)
	{//initialization.
		a[i]=4.0;aa[i]=4.0;
		b[i]=4.0;bb[i]=4.0;
		random[i]=0.1*double(rand()-RAND_MAX/2.0)*2/RAND_MAX;
		//random[i] = (m_vertices[i].z-(zmin+zmax)/2)/(zmax-zmin)*0.35+0.1;
	}
	for(time=0;time<m_nTime;time++)
	{
		for(i=0;i<n;i++)
		{
	    	a[i] = m_vertices[i].BCDArea*(aa[i]+tstep*RDTexturing_TurkF(aa[i],bb[i],random[i]));
			b[i] = m_vertices[i].BCDArea*(bb[i]+tstep*RDTexturing_TurkG(aa[i],bb[i],0));
		}
		m_RDTDiffRate = m_RDTp.Da;
	    nr_linbcg(n,a,aa,tstep,itol,tol,itmax,&iter,&err);
		m_RDTDiffRate = m_RDTp.Db;
		nr_linbcg(n,b,bb,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = bb[i];
		m_vertices[i].g = bb[i];
		m_vertices[i].b = bb[i];
	}
	float rmin,rmax;
	rmin = m_vertices[0].r;
	rmax = m_vertices[0].r;
	for(i=0;i<n;i++)
	{
		if(m_vertices[i].r<rmin) rmin=m_vertices[i].r;
		if(m_vertices[i].r>rmax) rmax=m_vertices[i].r;
	}
	for(i=0;i<n;i++)
		m_vertices[i].r-= rmin;
	for(i=0;i<n;i++)
		m_vertices[i].r/= (rmax-rmin);
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = 0.4;
		m_vertices[i].b = 0.8;
	}
/*
	ofstream ofp("venushead_text.txt");
	ofp<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
	{
		ofp<<aa[i]<<"   "<<m_vertices[i].r<<"\n";
	}
	ofp.close();
*/
	delete a;
	delete aa;
	delete b;
	delete bb;
	delete random;
}

void TriangularMesh::ImplicitL2RDTTuring_Brusselator()
{
//Turing's model with Brusselator's reaction term.
	unsigned long n,i;
	int time;
	double *a,*aa,*b,*bb,*random1,*random2;

	int itol,itmax,iter;
    double tol,err,tstep;

    tstep = m_tStep;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	n = m_vnum;
	a = new double[n];
	aa = new double[n];
	b = new double[n];
	bb = new double[n];
	random1 = new double[n];
	random2 = new double[n];
	float zmin,zmax;
	zmin = m_vertices[0].y;
	zmax = zmin;
	for(i=0;i<n;i++)
	{
		if(m_vertices[i].z<zmin) zmin=m_vertices[i].z;
		if(m_vertices[i].z>zmax) zmax=m_vertices[i].z;
	}
	for(i=0;i<n;i++)
	{//initialization.
		random1[i]=0.01*double(rand()-RAND_MAX/2.0)*2/RAND_MAX;
		random2[i]=0.01*double(rand()-RAND_MAX/2.0)*2/RAND_MAX;
		//random1[i] = (m_vertices[i].z-(zmin+zmax)/2)/(zmax-zmin)*0.35+0.1;
		a[i]=3.0*(1+random1[i]);aa[i]=3.0*(1+random1[i]);
		b[i]=3.0*(1+random1[i]);bb[i]=3.0*(1+random1[i]);
	}
	for(time=0;time<m_nTime;time++)
	{
		for(i=0;i<n;i++)
		{
	    	a[i] = m_vertices[i].BCDArea*(aa[i]+tstep*RDTexturing_BrusselatorF(aa[i],bb[i],
				m_RDTp.ar+(m_vertices[i].z-(zmin+zmax)/2)*2/(zmax-zmin)*0.2,
				m_RDTp.decay*(1+random1[i]),m_RDTp.gr*(1+random2[i])));
			b[i] = m_vertices[i].BCDArea*(bb[i]+tstep*RDTexturing_BrusselatorG(aa[i],bb[i],
				m_RDTp.ar+(m_vertices[i].z-(zmin+zmax)/2)*2/(zmax-zmin)*0.2,
				m_RDTp.decay*(1+random1[i]),m_RDTp.gr*(1+random2[i])));
		}
		m_RDTDiffRate = m_RDTp.Da;
	    nr_linbcg(n,a,aa,tstep,itol,tol,itmax,&iter,&err);
		m_RDTDiffRate = m_RDTp.Db;
		nr_linbcg(n,b,bb,tstep,itol,tol,itmax,&iter,&err);
	}
	for(i=0;i<n;i++)
	{
		m_vertices[i].r = bb[i];
		m_vertices[i].g = bb[i];
		m_vertices[i].b = bb[i];
	}
	float rmin,rmax;
	rmin = m_vertices[0].r;
	rmax = m_vertices[0].r;
	for(i=0;i<n;i++)
	{
		if(m_vertices[i].r<rmin) rmin=m_vertices[i].r;
		if(m_vertices[i].r>rmax) rmax=m_vertices[i].r;
	}
	for(i=0;i<n;i++)
		m_vertices[i].r-= rmin;
	for(i=0;i<n;i++)
		m_vertices[i].r/= (rmax-rmin);
	for(i=0;i<n;i++)
	{
		m_vertices[i].g = 0.5;
		m_vertices[i].b = 1;
	}
/*
	ofstream ofp("venushead_text.txt");
	ofp<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
	{
		ofp<<aa[i]<<"   "<<m_vertices[i].r<<"\n";
	}
	ofp.close();
*/
	delete a;
	delete aa;
	delete b;
	delete bb;
	delete random1;
	delete random2;
}

void TriangularMesh::ImplicitL2RDTTuring_GrayScott()
{

}

void TriangularMesh::ImplicitL2RDTTuring_GiererMein()
{

}

void TriangularMesh::ImplicitL2RDTMeinhardtStrips_Mein()
{
//realize the meinhardt's strips formating system which includes five chemicals.
}

void TriangularMesh::FindPPIGradient(TRIANGLE t, VECTOR3D *g, unsigned int dataType)
{
/*  在构建了mesh的PPIBGradient数组并保存之后，此函数不再使用。
 *  dataType:指示对何种数据在该三角形内求梯度，如颜色分量的r分量？g分量？b分量？
 *  dataType=1: r 分量；
 *  dataType=2: g 分量；
 *  dataType=3: b 分量；
 *  dataType=4: 或许是contour 数据；
 */
	double h;
	VECTOR3D g0,g1,g2;
	FindPPIBGradient(t.ver0,t,&g0,&h);
	FindPPIBGradient(t.ver1,t,&g1,&h);
	FindPPIBGradient(t.ver2,t,&g2,&h);
	switch(dataType)
	{
	case 1:
	    (*g).x = m_vertices[t.ver0].r*g0.x+m_vertices[t.ver1].r*g1.x+m_vertices[t.ver2].r*g2.x;
	    (*g).y = m_vertices[t.ver0].r*g0.y+m_vertices[t.ver1].r*g1.y+m_vertices[t.ver2].r*g2.y;
	    (*g).z = m_vertices[t.ver0].r*g0.z+m_vertices[t.ver1].r*g1.z+m_vertices[t.ver2].r*g2.z;
		break;
	case 2:
		(*g).x = m_vertices[t.ver0].g*g0.x+m_vertices[t.ver1].g*g1.x+m_vertices[t.ver2].g*g2.x;
	    (*g).y = m_vertices[t.ver0].g*g0.y+m_vertices[t.ver1].g*g1.y+m_vertices[t.ver2].g*g2.y;
	    (*g).z = m_vertices[t.ver0].g*g0.z+m_vertices[t.ver1].g*g1.z+m_vertices[t.ver2].g*g2.z;
		break;
	case 3:
		(*g).x = m_vertices[t.ver0].b*g0.x+m_vertices[t.ver1].b*g1.x+m_vertices[t.ver2].b*g2.x;
	    (*g).y = m_vertices[t.ver0].b*g0.y+m_vertices[t.ver1].b*g1.y+m_vertices[t.ver2].b*g2.y;
	    (*g).z = m_vertices[t.ver0].b*g0.z+m_vertices[t.ver1].b*g1.z+m_vertices[t.ver2].b*g2.z;
		break;
	case 4:
		break;
	default:
		(*g).x = 0;
		(*g).y = 0;
		(*g).z = 0;
		break;
	}
}

VECTOR3D TriangularMesh::FindPPIGradient(int it, unsigned int dataType)
{
/*  from the index of the triangle to compute the ppi gradient of that.
 *  dataType:指示对何种数据在该三角形内求梯度，如颜色分量的r分量？g分量？b分量？
 *  dataType=1: r 分量；
 *  dataType=2: g 分量；
 *  dataType=3: b 分量；;
 *  dataType=4: 或许是contour 数据；
 */
	VECTOR3D g;
	switch(dataType)
	{
	case 1:
	    g.x = m_vertices[m_triangles[it].ver0].r*m_TPPIBG[it].v0.x+m_vertices[m_triangles[it].ver1].r*m_TPPIBG[it].v1.x+
		    m_vertices[m_triangles[it].ver2].r*m_TPPIBG[it].v2.x;
	    g.y = m_vertices[m_triangles[it].ver0].r*m_TPPIBG[it].v0.y+m_vertices[m_triangles[it].ver1].r*m_TPPIBG[it].v1.y+
		    m_vertices[m_triangles[it].ver2].r*m_TPPIBG[it].v2.y;
	    g.z = m_vertices[m_triangles[it].ver0].r*m_TPPIBG[it].v0.z+m_vertices[m_triangles[it].ver1].r*m_TPPIBG[it].v1.z+
		    m_vertices[m_triangles[it].ver2].r*m_TPPIBG[it].v2.z;
		break;
	case 2:
		g.x = m_vertices[m_triangles[it].ver0].g*m_TPPIBG[it].v0.x+m_vertices[m_triangles[it].ver1].g*m_TPPIBG[it].v1.x+
		    m_vertices[m_triangles[it].ver2].g*m_TPPIBG[it].v2.x;
	    g.y = m_vertices[m_triangles[it].ver0].g*m_TPPIBG[it].v0.y+m_vertices[m_triangles[it].ver1].g*m_TPPIBG[it].v1.y+
		    m_vertices[m_triangles[it].ver2].g*m_TPPIBG[it].v2.y;
	    g.z = m_vertices[m_triangles[it].ver0].g*m_TPPIBG[it].v0.z+m_vertices[m_triangles[it].ver1].g*m_TPPIBG[it].v1.z+
		    m_vertices[m_triangles[it].ver2].g*m_TPPIBG[it].v2.z;
		break;
	case 3:
		g.x = m_vertices[m_triangles[it].ver0].b*m_TPPIBG[it].v0.x+m_vertices[m_triangles[it].ver1].b*m_TPPIBG[it].v1.x+
		    m_vertices[m_triangles[it].ver2].b*m_TPPIBG[it].v2.x;
	    g.y = m_vertices[m_triangles[it].ver0].b*m_TPPIBG[it].v0.y+m_vertices[m_triangles[it].ver1].b*m_TPPIBG[it].v1.y+
		    m_vertices[m_triangles[it].ver2].b*m_TPPIBG[it].v2.y;
	    g.z = m_vertices[m_triangles[it].ver0].b*m_TPPIBG[it].v0.z+m_vertices[m_triangles[it].ver1].b*m_TPPIBG[it].v1.z+
		    m_vertices[m_triangles[it].ver2].b*m_TPPIBG[it].v2.z;
		break;
	case 4:
		break;
	default:
		g.x = 0;
		g.y = 0;
		g.z = 0;
		break;
	}
	return g;
}

void TriangularMesh::GetPPIGradients(long n, double q[])
{
//be sure that n is the number of the vertices of the whole mesh.
//and the index of q is the same as the index of m_vertices.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		m_triangles[i].grad.x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		m_triangles[i].grad.y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		m_triangles[i].grad.z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

void TriangularMesh::GetPPIGradients1(long n, double q[])
{
//be sure that n is the number of the vertices of the whole mesh.
//and the index of q is the same as the index of m_vertices.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		m_triangles[i].grad1.x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		m_triangles[i].grad1.y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		m_triangles[i].grad1.z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

void TriangularMesh::GetPPIGradients2(long n, double q[])
{
//be sure that n is the number of the vertices of the whole mesh.
//and the index of q is the same as the index of m_vertices.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		m_triangles[i].grad2.x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		m_triangles[i].grad2.y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		m_triangles[i].grad2.z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

void TriangularMesh::BuildTrianglesArea()
{
	m_trianglesArea.clear();
	int i;
	double area;
	TRIANGLE t;
	VERTEX3D v0,v1,v2;
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
		v0 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
		v1 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
		v2 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
		area = Norm(CrossProduct(GetVector3Dfrom2Vertices(v0,v1),
			GetVector3Dfrom2Vertices(v0,v2)));
		m_trianglesArea.push_back(area);
	}
}

void TriangularMesh::NormalizeVector3D(double *xx, double *yy, double *zz)
{
	VECTOR3D v;
	v.x = *xx;
	v.y = *yy;
	v.z = *zz;
	NormalizeVector3D(&v);
	*xx = v.x;
	*yy = v.y;
	*zz = v.z;
}

void TriangularMesh::GetPrincipleDirection(int iver, VECTOR3D *e1, VECTOR3D *e2)
{
//compute the 2 unit principle directions of vertex indexed by iver.
	if(iver<0||iver>m_vnum)
	{
		(*e1).x = 0.0f;
		(*e1).y = 0.0f;
		(*e1).z = 0.0f;
		(*e2).x = 0.0f;
		(*e2).y = 0.0f;
		(*e2).z = 0.0f;
		return;
	}
	FRAME f;
    ONERING or;
	MATRIX3BY3 m;
	VECTOR3D ve,dir;
	double alpha,beta;
	double a,b,c,r1,r2,r3;
	double weight = 1.0;
	double hN,max;
	int i;
	or = m_ONERINGS[iver];
	f.p = m_vertices[iver];
	f.e3.x = m_vertices[iver].normal_x;
    f.e3.y = m_vertices[iver].normal_y;
	f.e3.z = m_vertices[iver].normal_z;//e3 stands for the normal direction.
	NormalizeVector3D(&f.e3);//
    ve = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[iver]),GetVertex3DFromPoint3D(m_vertices[or[0].iver]));
	f.e1 = GetCoProjection(ve,f.e3);
	NormalizeVector3D(&f.e1);//单位化。
	f.e2 = CrossProduct(f.e3,f.e1);
	NormalizeVector3D(&f.e2);//
	//up to now, a frame at current vertex is constructed.
	m.a11=0;m.a12=0;m.a13=0;
	m.a21=0;m.a22=0;m.a23=0;
	m.a31=0;m.a32=0;m.a33=0;
	r1=0;r2=0;r3=0;//
	max = 0;
	for(i=0;i<or.size();i++)
	{
		ve = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[iver]),GetVertex3DFromPoint3D(m_vertices[or[i].iver]));
		if(max<Norm(ve))
			max = Norm(ve);
	}
	for(i=0;i<or.size();i++)
	{
		ve = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[iver]),GetVertex3DFromPoint3D(m_vertices[or[i].iver]));
		hN = GetProjection(ve,f.e3);//ve沿f.e3的分量。
		dir = GetCoProjection(ve,f.e3);
		alpha = DotProduct(dir,f.e1);
		beta  = DotProduct(dir,f.e2);
		weight = exp(-Norm(ve)/max);//
		m.a11+= weight*POWER4(alpha);
		m.a12+= weight*beta*POWER3(alpha);
		m.a13+= weight*POWER(alpha)*POWER(beta);
		r1   += weight*hN*POWER(alpha);
		m.a21+= weight*POWER3(alpha)*beta;
		m.a22+= weight*POWER(beta)*POWER(alpha);
		m.a23+= weight*alpha*POWER3(beta);
		r2   += weight*hN*alpha*beta;
		m.a31+= weight*POWER(alpha)*POWER(beta);
		m.a32+= weight*alpha*POWER3(beta);
		m.a33+= weight*POWER4(beta);
		r3   += weight*hN*POWER(beta);
	}
	//then solve for the curvature tensor.
	Cramer3(m,r1,r2,r3,&a,&b,&c);//
	if(a==0&b==0&c==0)
	{//this case is not considered in detail, a plane?
		(*e1).x = 0.0;
		(*e1).y = 0.0;
		(*e1).z = 0.0;
		(*e2).x = 0.0;
		(*e2).y = 0.0;
		(*e2).z = 0.0;
	}
	else
	{
		if(b!=0)
		{
			double lamda=c-a-sqrt(POWER(a-c)+b*b);
			(*e1).x = b*f.e1.x+lamda*f.e2.x;
			(*e1).y = b*f.e1.y+lamda*f.e2.y;
			(*e1).z = b*f.e1.z+lamda*f.e2.z;
			(*e2).x = lamda*f.e1.x-b*f.e2.x;
			(*e2).y = lamda*f.e1.y-b*f.e2.y;
			(*e2).z = lamda*f.e1.z-b*f.e2.z;
		}
		else
		{
			if(a<c)
			{
				*e1 = f.e1;//
				*e2 = f.e2;
			}
			else if(a>c)
			{
				*e1 = f.e2;//
				*e2 = f.e1;
			}
			else
			{
				(*e1).x = 0;
				(*e1).y = 0;
				(*e1).z = 0;
				(*e2).x = 0;
				(*e2).y = 0;
				(*e2).z = 0;
			}
		}
	}
	if(Norm(*e1)>0)
		NormalizeVector3D(e1);
	if(Norm(*e2)>0)
		NormalizeVector3D(e2);
}

double TriangularMesh::GetProjection(VECTOR3D v, VECTOR3D r)
{
//compute the coordinate of vector v with the reference vector r.
//r is a unit vector.
	return DotProduct(v,r);
}

VECTOR3D TriangularMesh::GetCoProjection(VECTOR3D v, VECTOR3D r)
{
	double coor;
	VECTOR3D cop;
	coor = GetProjection(v,r);
	cop.x = v.x-coor*r.x;
	cop.y = v.y-coor*r.y;
	cop.z = v.z-coor*r.z;
	return cop;
}

void TriangularMesh::Cramer3(MATRIX3BY3 m, double r[], double x[])
{//不完善的cramer law.
	MATRIX3BY3 tmp;
	if(Determinant3(m)==0)
	{
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
	}
	else
	{
		tmp = m;
		tmp.a11 = r[0];
		tmp.a21 = r[1];
		tmp.a31 = r[2];
		x[0] = Determinant3(tmp)/Determinant3(m);
		tmp = m;
		tmp.a12 = r[0];
		tmp.a22 = r[1];
		tmp.a32 = r[2];
		x[1] = Determinant3(tmp)/Determinant3(m);
		tmp = m;
		tmp.a13 = r[0];
		tmp.a23 = r[1];
		tmp.a33 = r[2];
		x[2] = Determinant3(tmp)/Determinant3(m);
	}
}

double TriangularMesh::Determinant3(MATRIX3BY3 m)
{
/*        a11  a12  a13
          a21  a22  a23
		  a31  a32  a33
 */
	return m.a11*(m.a22*m.a33-m.a32*m.a23)
		-m.a12*(m.a21*m.a33-m.a31*m.a23)
		+m.a13*(m.a21*m.a32-m.a31*m.a22);
}

void TriangularMesh::Cramer3(MATRIX3BY3 m, double r1, double r2, double r3, double *x1, double *x2, double *x3)
{//不完善的cramer law.
	MATRIX3BY3 tmp;
	if(Determinant3(m)==0)
	{
		*x1 = 0;
		*x2 = 0;
		*x3 = 0;
	}
	else
	{
		tmp = m;
		tmp.a11 = r1;
		tmp.a21 = r2;
		tmp.a31 = r3;
		*x1 = Determinant3(tmp)/Determinant3(m);
		tmp = m;
		tmp.a12 = r1;
		tmp.a22 = r2;
		tmp.a32 = r3;
		*x2 = Determinant3(tmp)/Determinant3(m);
		tmp = m;
		tmp.a13 = r1;
		tmp.a23 = r2;
		tmp.a33 = r3;
		*x3 = Determinant3(tmp)/Determinant3(m);
	}
}

void TriangularMesh::BuildPrincipleDirections()
{
	int i;
	for(i=0;i<m_vnum;i++)
		GetPrincipleDirection(i,&(m_vertices[i].pde1),&(m_vertices[i].pde2));
	ofstream ofp("pde1.txt");
	ofp<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
		ofp<<m_vertices[i].pde1.x<<"   "<<m_vertices[i].pde1.y<<"   "<<
		m_vertices[i].pde1.z<<"   \n";
	ofp.close();

	ofstream ofp1("pde2normaldp.txt");
	ofp1<<m_vnum<<"\n";
	for(i=0;i<m_vnum;i++)
		ofp1<<(m_vertices[i].pde2.x*m_vertices[i].normal_x+
		m_vertices[i].pde2.y*m_vertices[i].normal_y+
		m_vertices[i].pde2.z*m_vertices[i].normal_z)<<"\n";
	ofp1.close();
}

void TriangularMesh::WriteDataToPovray()
{
	WriteMeshToPovray();
//	WriteContourToPovray();
//	WriteWeightedCurveLength();
//	WriteWeightedCurveLength_Forobservation();
}

void TriangularMesh::WriteWeightedCurveLength()
{
	int i=0;
	vector<double>::iterator dite1, dite2, dite3;
	ofstream ofp("WeightedCurveLength_Timesteps_.txt");
	for(dite1=m_WeightedCurveLength.begin(),dite2=m_tSteps.begin(),dite3=m_CPUtime.begin();
	dite1<m_WeightedCurveLength.end(),dite2<m_tSteps.end(),dite3<m_CPUtime.end();dite1++,dite2++,dite3++)
	{
		ofp<<i<<":   "<<*dite1<<";   "<<*dite2<<";   "<<*dite3<<"\n";
		i++;
	}
	ofp.close();
}

void TriangularMesh::WriteContourToPovray()
{
	int i;
	POINT3d tmp1,tmp2,tmp3, linend1,linend2;
	ofstream ofp("ContourPr_.txt");
	ofp<<"#declare Contour=\n";
	ofp<<"union{\n";
	for(i=0;i<m_trinum;i++)
	{
	    tmp1 = m_vertices[m_triangles[i].ver0];
	    tmp2 = m_vertices[m_triangles[i].ver1];
	    tmp3 = m_vertices[m_triangles[i].ver2];
	    if(ContouringATriangle(tmp1,tmp2,tmp3,&linend1,&linend2))
		{
		    ofp<<"cylinder{ <"<<linend1.x<<","<<linend1.y<<","<<linend1.z<<" >, "
				<<"<"<<linend2.x<<","<<linend2.y<<","<<linend2.z<<" >,"
				<<"rad_cyl \n"
				<<"pigment {ContourPig} \n"
				<<"finish {DefaultFinish}"
				<<"}\n";
		}   
	}
	ofp<<"}";
	ofp.close();
}

void TriangularMesh::WriteMeshToPovray()
{
	int i;
	POINT3d p;
	TRIANGLE t;
	ofstream ofp("MeshPr_.txt");
//	ofp<<"object{\n";
	ofp<<"#declare Mesh=\n";
	ofp<<"union{\n";
	ofp<<"mesh2{\n";
    /* Begin vertex_vectors */
	ofp<<"     vertex_vectors{\n";
	ofp<<"        "<<m_vnum<<",\n       ";//顶点数目.
	for(i=0;i<m_vnum;i++)
	{
		p = m_vertices[i];
		ofp<<" <"<<p.x<<","<<p.y<<","<<p.z<<">,";
		if(i%3 == 2)
		    ofp<<"\n        ";
	}
	ofp<<"\n}\n";
    /* End vertex_vectors */
    /* Begin vertex_normals */
	ofp<<"     normal_vectors{\n";
	ofp<<"        "<<m_vnum<<",\n       ";//顶点数目.
	for(i=0;i<m_vnum;i++)
	{
		p = m_vertices[i];
		ofp<<" <"<<p.normal_x<<","<<p.normal_y<<","<<p.normal_z<<">,";
		if(i%3 == 2)
		    ofp<<"\n        ";
	}
	ofp<<"\n}\n";
	/* End vertex_normals */
    /* Begin texture_list */
	ofp<<"     texture_list{\n";
    ofp<<"        "<<m_vnum<<",\n        ";//顶点数目.
	for(i=0;i<m_vnum;i++)
	{
		p = m_vertices[i];
		ofp<<"texture{pigment{rgb <"<<p.r<<","<<p.g<<","<<p.b<<">}}";// for colored meshes.
//		ofp<<"texture{pigment{MeshPig}}";// for noncolored meshes.
		ofp<<"\n        ";
	}
	ofp<<"}\n";
    /* End texture_list */
    /* Begin face_indices */
	ofp<<"face_indices{\n";
	ofp<<"        "<<m_trinum<<",\n         ";//三角片数目.
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
		ofp<<" <"<<t.ver0<<","<<t.ver1<<","<<t.ver2<<">,"
			<<t.ver0<<","<<t.ver1<<","<<t.ver2<<",   ";
		if(i%2 == 1)
    		ofp<<"\n             ";
	}
	ofp<<"}\n";
    /* End face_indices */
	ofp<<"}";
//	ofp<<"\n}";
	ofp<<"\n}";
	ofp.close();
}

void TriangularMesh::BuildONEDISKSANISOTROPIC()
{
	float lambda1,lambda2;//eigenvalues of diffusion tensor along 2 eigendirections.
	VECTOR3D trinormal,e1,e2;//normal direction of current triangle and 2 eigendirections of the diffusion tensor stricted in this triangle.
    int i,j;
	float alpha;
	DISKTRIANGLEANISOTROPIC dta;
	DISKANISOTROPIC da;
	EDGE e;
	TRIANGLE t;
	m_ONEDISKSANISOTROPIC.clear();
	for(i=0;i<m_onedisks.size();i++)
	{
		da.clear();
		e.ver0 = i;
		for(j=0;j<m_onedisks[i].size();j++)
		{
			t = m_triangles[m_onedisks[i][j]];
			Design_e1e2(t, &e1, &e2);
			//e1,e2 are 2 orthonormal basis vectors of the underlying space of this triangle.
			lambda1 = 1;
			lambda2 = 1;
			//2 eigenvalues along the 2 eigendirections e1 and e2.
			dta.iDiskCent = i;
			dta.itri = m_onedisks[i][j];
			dta.d1ver0co = 0;
			dta.d1ver1co = 0;
			dta.d1ver2co = 0;
			dta.d2ver0co = 0;
			dta.d2ver1co = 0;
			dta.d2ver2co = 0;
			if(t.ver0==i)
			{
				e.ver1 = t.ver1;
				dta.d1ver0co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver1co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver2co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver0co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver1co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver2co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				e.ver1 = t.ver2;
				dta.d1ver0co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver1co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver2co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver0co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver1co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver2co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
			}
			else if(t.ver1==i)
			{
				e.ver1 = t.ver0;
				dta.d1ver0co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver1co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver2co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver0co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver1co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver2co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				e.ver1 = t.ver2;
				dta.d1ver0co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver1co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver2co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver0co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver1co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver2co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
			}
			else
			{
				e.ver1 = t.ver0;
				dta.d1ver0co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver1co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver2co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver0co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver1co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver2co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				e.ver1 = t.ver1;
				dta.d1ver0co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver1co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d1ver2co+= lambda1*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e1)*DotProduct(e1,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver0co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v0,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver1co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v1,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
				dta.d2ver2co+= lambda2*DotProduct(m_TPPIBG[m_onedisks[i][j]].v2,e2)*DotProduct(e2,GetBCDNormal(i,e,t))*Norm(
					GetVector3Dfrom2Vertices(GetBaryCenter(e),GetBaryCenter(t)));
			}
			da.push_back(dta);
		}
		m_ONEDISKSANISOTROPIC.push_back(da);
	}
/*
	ofstream ofp2("onedisksaniso.txt");
	ofp2<<m_ONEDISKSANISOTROPIC.size()<<"\n";
	for(i=0;i<m_ONEDISKSANISOTROPIC.size();i++)
	{
		ofp2<<i<<"\n";
		for(j=0;j<m_ONEDISKSANISOTROPIC[i].size();j++)
		{
			t = m_triangles[m_ONEDISKSANISOTROPIC[i][j].itri];
		    ofp2<<"   "<<t.ver0<<":"<<m_ONEDISKSANISOTROPIC[i][j].d1ver0co<<";  "
				<<"   "<<t.ver1<<":"<<m_ONEDISKSANISOTROPIC[i][j].d1ver1co<<";  "
				<<"   "<<t.ver2<<":"<<m_ONEDISKSANISOTROPIC[i][j].d1ver2co<<";  "
				<<"   "<<t.ver0<<":"<<m_ONEDISKSANISOTROPIC[i][j].d2ver0co<<";  "
				<<"   "<<t.ver1<<":"<<m_ONEDISKSANISOTROPIC[i][j].d2ver1co<<";  "
				<<"   "<<t.ver2<<":"<<m_ONEDISKSANISOTROPIC[i][j].d2ver2co<<"\n";
		}
		ofp2<<"\n";
	}
	ofp2.close();
*/
}

void TriangularMesh::DrawContour()
{
	//int i;
	//POINT3d tmp1,tmp2,tmp3, linend1,linend2;
	//glLineWidth(6.0);
	//glBegin(GL_LINES);
	//   for(i=0;i<m_trinum;i++)
	//   {
	//	   tmp1 = m_vertices[m_triangles[i].ver0];
	//	   tmp2 = m_vertices[m_triangles[i].ver1];
	//	   tmp3 = m_vertices[m_triangles[i].ver2];
	//	   if(ContouringATriangle(tmp1,tmp2,tmp3,&linend1,&linend2))
	//	   {
	//		   glColor3f(0.0f,1.0f,0.0f);
	//	       glVertex3d(linend1.x,linend1.y,linend1.z);
	//		   glColor3f(0.0f,1.0f,0.0f);
	//		   glVertex3d(linend2.x,linend2.y,linend2.z);
	//	   }
	//   }
	//glEnd();
	//glLineWidth(1.0);//reset the linewidth to the default value.
}

void TriangularMesh::DrawContour_MulRegionLabelling()
{
	//int i;
	//POINT3d tmp1,tmp2,tmp3;
	//VERTEX3D bct,bce1,bce2,bce3;
	//glLineWidth(3.0);
	//glBegin(GL_LINES);
	//   for(i=0;i<m_trinum;i++)
	//   {
	//	   tmp1 = m_vertices[m_triangles[i].ver0];
	//	   tmp2 = m_vertices[m_triangles[i].ver1];
	//	   tmp3 = m_vertices[m_triangles[i].ver2];
	//	   bct.x = (tmp1.x+tmp2.x+tmp3.x)/3;
	//	   bct.y = (tmp1.y+tmp2.y+tmp3.y)/3;
	//	   bct.z = (tmp1.z+tmp2.z+tmp3.z)/3;
	//	   bce1.x = (tmp1.x+tmp2.x)/2; bce1.y = (tmp1.y+tmp2.y)/2; bce1.z = (tmp1.z+tmp2.z)/2;
	//	   bce2.x = (tmp2.x+tmp3.x)/2; bce2.y = (tmp2.y+tmp3.y)/2; bce2.z = (tmp2.z+tmp3.z)/2;
	//	   bce3.x = (tmp3.x+tmp1.x)/2; bce3.y = (tmp3.y+tmp1.y)/2; bce3.z = (tmp3.z+tmp1.z)/2;
	//	   if(tmp1.RegionId!=tmp2.RegionId)
	//	   {
	//		   glColor3f(0.0f,1.0f,0.0f);
	//		   glVertex3d(bct.x,bct.y,bct.z);
	//		   glVertex3d(bce1.x,bce1.y,bce1.z);
	//	   }
	//	   if(tmp2.RegionId!=tmp3.RegionId)
	//	   {
	//		   glColor3f(0.0f,1.0f,0.0f);
	//		   glVertex3d(bct.x,bct.y,bct.z);
	//		   glVertex3d(bce2.x,bce2.y,bce2.z);
	//	   }
	//	   if(tmp3.RegionId!=tmp1.RegionId)
	//	   {
	//		   glColor3f(0.0f,1.0f,0.0f);
	//		   glVertex3d(bct.x,bct.y,bct.z);
	//		   glVertex3d(bce3.x,bce3.y,bce3.z);
	//	   }
	//   }
	//glEnd();
	//glLineWidth(1.0);
}

void TriangularMesh::DrawMesh()
{
	//int i;
	//POINT3d tmp;
	//glBegin(GL_TRIANGLES);
	//   for(i=0;i<m_trinum;i++)
	//   {
	//	   tmp = m_vertices[m_triangles[i].ver0];
	//	   glColor3f(tmp.r,tmp.g,tmp.b);
	//	   glNormal3d(tmp.normal_x,tmp.normal_y,tmp.normal_z);
	//	   glVertex3d(tmp.x,tmp.y,tmp.z);
	//	   tmp = m_vertices[m_triangles[i].ver1];
	//	   glColor3f(tmp.r,tmp.g,tmp.b);
	//	   glNormal3d(tmp.normal_x,tmp.normal_y,tmp.normal_z);
	//	   glVertex3d(tmp.x,tmp.y,tmp.z);
	//	   tmp = m_vertices[m_triangles[i].ver2];
	//	   glColor3f(tmp.r,tmp.g,tmp.b);
	//	   glNormal3d(tmp.normal_x,tmp.normal_y,tmp.normal_z);
	//	   glVertex3d(tmp.x,tmp.y,tmp.z);
	//   }
	//glEnd();
}

bool TriangularMesh::ContouringATriangle(POINT3d tmp1, POINT3d tmp2, POINT3d tmp3, POINT3d *linend1, POINT3d *linend2)
{
	if(tmp1.contour>0 && tmp2.contour>0 && tmp3.contour>0)
		return false;
	if(tmp1.contour<0 && tmp2.contour<0 && tmp3.contour<0)
		return false;
	float ratio;
	if(tmp1.contour*tmp2.contour>0)
	{
		ratio = tmp3.contour/(tmp3.contour-tmp1.contour);
		(*linend1).x = (1-ratio)*tmp3.x+ratio*tmp1.x;
		(*linend1).y = (1-ratio)*tmp3.y+ratio*tmp1.y;
		(*linend1).z = (1-ratio)*tmp3.z+ratio*tmp1.z;

		ratio = tmp3.contour/(tmp3.contour-tmp2.contour);
		(*linend2).x = (1-ratio)*tmp3.x+ratio*tmp2.x;
		(*linend2).y = (1-ratio)*tmp3.y+ratio*tmp2.y;
		(*linend2).z = (1-ratio)*tmp3.z+ratio*tmp2.z;
		return true;
	}
	if(tmp2.contour*tmp3.contour>0)
	{
		ratio = tmp1.contour/(tmp1.contour-tmp3.contour);
		(*linend1).x = (1-ratio)*tmp1.x+ratio*tmp3.x;
		(*linend1).y = (1-ratio)*tmp1.y+ratio*tmp3.y;
		(*linend1).z = (1-ratio)*tmp1.z+ratio*tmp3.z;

		ratio = tmp1.contour/(tmp1.contour-tmp2.contour);
		(*linend2).x = (1-ratio)*tmp1.x+ratio*tmp2.x;
		(*linend2).y = (1-ratio)*tmp1.y+ratio*tmp2.y;
		(*linend2).z = (1-ratio)*tmp1.z+ratio*tmp2.z;
		return true;
	}
	if(tmp3.contour*tmp1.contour>0)
	{
		ratio = tmp2.contour/(tmp2.contour-tmp3.contour);
		(*linend1).x = (1-ratio)*tmp2.x+ratio*tmp3.x;
		(*linend1).y = (1-ratio)*tmp2.y+ratio*tmp3.y;
		(*linend1).z = (1-ratio)*tmp2.z+ratio*tmp3.z;

		ratio = tmp2.contour/(tmp2.contour-tmp1.contour);
		(*linend2).x = (1-ratio)*tmp2.x+ratio*tmp1.x;
		(*linend2).y = (1-ratio)*tmp2.y+ratio*tmp1.y;
		(*linend2).z = (1-ratio)*tmp2.z+ratio*tmp1.z;
		return true;
	}
}

void TriangularMesh::SetObjectDrawing(unsigned char od)
{
	m_ObjectDrawing = od;
}

void TriangularMesh::ShiftObjectDrawing()
{
	m_ObjectDrawing++;
	m_ObjectDrawing = m_ObjectDrawing/3;
}

void TriangularMesh::BuildEdgeIndicator(unsigned char type)
{
/* to build edge indicator of data of type type.
 * type = 1: no weight, EdgeIndicator = 1;
 * type = 2: gray image or averaged (r,g,b) of a color image;
 * type = 3: the r component of a color image;
 * type = 4: the g component of a color image;
 * type = 5: the b component of a color image;
 * type = 6: the l2 o grad;
 * type = 7: curvature.......
 * an edge indicator function g should be chosen.
 */
	long i;
	float k;
	double grey0,grey1,grey2;
	VECTOR3D gradr,gradg,gradb;
	switch(type){
	case 1:
		for(i=0;i<m_trinum;i++)
			m_triangles[i].edgeindicator = 1;
		break;
	case 2:
    	for(i=0;i<m_trinum;i++)
		{
		    grey0 = m_vertices[m_triangles[i].ver0].r+m_vertices[m_triangles[i].ver0].g+m_vertices[m_triangles[i].ver0].b;
			grey0/= 3;
			grey1 = m_vertices[m_triangles[i].ver1].r+m_vertices[m_triangles[i].ver1].g+m_vertices[m_triangles[i].ver1].b;
			grey1/= 3;
			grey2 = m_vertices[m_triangles[i].ver2].r+m_vertices[m_triangles[i].ver2].g+m_vertices[m_triangles[i].ver2].b;
			grey2/= 3;
			gradr.x = grey0*m_TPPIBG[i].v0.x+grey1*m_TPPIBG[i].v1.x+grey2*m_TPPIBG[i].v2.x;
			gradr.y = grey0*m_TPPIBG[i].v0.y+grey1*m_TPPIBG[i].v1.y+grey2*m_TPPIBG[i].v2.y;
			gradr.z = grey0*m_TPPIBG[i].v0.z+grey1*m_TPPIBG[i].v1.z+grey2*m_TPPIBG[i].v2.z;
			m_triangles[i].edgeindicator = 1/(1+10*POWER(Norm(gradr)));
		}
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
//		k = 7.0;  // for bunny.
		k = 7.0; // for horse.
		for(i=0;i<m_trinum;i++)
		{
		    grey0 = m_vertices[m_triangles[i].ver0].r;
			grey1 = m_vertices[m_triangles[i].ver1].r;
			grey2 = m_vertices[m_triangles[i].ver2].r;
			gradr.x = grey0*m_TPPIBG[i].v0.x+grey1*m_TPPIBG[i].v1.x+grey2*m_TPPIBG[i].v2.x;
			gradr.y = grey0*m_TPPIBG[i].v0.y+grey1*m_TPPIBG[i].v1.y+grey2*m_TPPIBG[i].v2.y;
			gradr.z = grey0*m_TPPIBG[i].v0.z+grey1*m_TPPIBG[i].v1.z+grey2*m_TPPIBG[i].v2.z;
			grey0 = m_vertices[m_triangles[i].ver0].g;
			grey1 = m_vertices[m_triangles[i].ver1].g;
			grey2 = m_vertices[m_triangles[i].ver2].g;
			gradg.x = grey0*m_TPPIBG[i].v0.x+grey1*m_TPPIBG[i].v1.x+grey2*m_TPPIBG[i].v2.x;
			gradg.y = grey0*m_TPPIBG[i].v0.y+grey1*m_TPPIBG[i].v1.y+grey2*m_TPPIBG[i].v2.y;
			gradg.z = grey0*m_TPPIBG[i].v0.z+grey1*m_TPPIBG[i].v1.z+grey2*m_TPPIBG[i].v2.z;
			grey0 = m_vertices[m_triangles[i].ver0].b;
			grey1 = m_vertices[m_triangles[i].ver1].b;
			grey2 = m_vertices[m_triangles[i].ver2].b;
			gradb.x = grey0*m_TPPIBG[i].v0.x+grey1*m_TPPIBG[i].v1.x+grey2*m_TPPIBG[i].v2.x;
			gradb.y = grey0*m_TPPIBG[i].v0.y+grey1*m_TPPIBG[i].v1.y+grey2*m_TPPIBG[i].v2.y;
			gradb.z = grey0*m_TPPIBG[i].v0.z+grey1*m_TPPIBG[i].v1.z+grey2*m_TPPIBG[i].v2.z;

			m_triangles[i].edgeindicator = 1/(1+k*POWER(Norm(gradr))+k*POWER(Norm(gradg))+k*POWER(Norm(gradb)));
		}
		break;
	case 7:
		break;
	default:
		break;
	}
}

void TriangularMesh::ComputeAveragedVertexAbGradient()
{
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double abgrad;
	for(i=0;i<m_vnum;i++)
	{
		abgrad = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			dt = d[j];
			abgrad+= Norm(m_triangles[dt.itri].grad)*m_trianglesArea[dt.itri]/3;
		}
		abgrad/= m_vertices[i].BCDArea;
		m_vertices[i].averagedAbGradient = abgrad;
	}
}

double TriangularMesh::ComputeWeightedCurveLength()
{
	double len=0;
	int i;
	POINT3d tmp1,tmp2,tmp3, linend1,linend2;
	for(i=0;i<m_trinum;i++)
	{
	    tmp1 = m_vertices[m_triangles[i].ver0];
		tmp2 = m_vertices[m_triangles[i].ver1];
		tmp3 = m_vertices[m_triangles[i].ver2];
		if(ContouringATriangle(tmp1,tmp2,tmp3,&linend1,&linend2))
		{
		    len+=Norm(GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(linend1),
				GetVertex3DFromPoint3D(linend2)))*m_triangles[i].edgeindicator;
		}
	}
	return len;
}

void TriangularMesh::Reinitialize_Contour(double time)
{
	int i;

	for(i=0;i<m_vnum;i++)
	{
		m_vertices[i].r = 0.5;
		m_vertices[i].g = 0.5;
		m_vertices[i].b = 0.5;
	}

	for(i=0;i<m_vnum;i++)
	{
		if(m_vertices[i].contour<0)
			m_vertices[i].contour = -1;
		else if(m_vertices[i].contour>0)
			m_vertices[i].contour = 1;
	}

	ImplicitL2Diff_Contour(time);

	for(i=0;i<m_vnum;i++)
		if(fabs(m_vertices[i].contour)<0.9f)
		{
			m_vertices[i].r = 0;
			m_vertices[i].g = 0;
			m_vertices[i].b = 1;
		}
}

void TriangularMesh::FMM_LocalSolver_Acute(TRIANGLEWITHCOOR twc, double u0, double u1, double *u2, double F)
{
/*  The local solver of fast marching methods for acute angle.
 *  The twc is the given triangle with vertices ver0, ver1, and ver2.
 *  The u0, u1, u2 are corresponding solution at three vertices.
 *  Therein u0 and u1, u2 are all given, but u2 is to be updated in this function.
 *  F is as defined in the Eikonal equation.
 *  The angle at ver2 should be acute.
 *  Also u1 > u0. 
 */
    
	VECTOR3D v02, v12;
	v02 = GetVector3Dfrom2Vertices(twc.ver2,twc.ver0);
	v12 = GetVector3Dfrom2Vertices(twc.ver2,twc.ver1);

    double u = u1-u0;

	// the coefficients for the quadratic equation.
	double a2, a1, a0;
	a2 = DotProduct(v02,v02) + DotProduct(v12,v12) - 2*DotProduct(v02,v12);
    a1 = 2*u*(DotProduct(v02,v12) - DotProduct(v02,v02));
	a0 = u*u*DotProduct(v02,v02)
		- F*F*(DotProduct(v12,v12)*DotProduct(v02,v02)-DotProduct(v02,v12)*DotProduct(v02,v12));
    // solve the quadratic a2*t_2 + a1*t + a0 = 0.
	double delta;
	if(delta=a1*a1-4*a2*a0 < 0)
	{
		*u2 = MIN(*u2, MIN(F*sqrt(DotProduct(v02,v02))+u0, F*sqrt(DotProduct(v12,v12))+u1));
		return;
	}

	double t1, t2;
	t1 = (-a1+sqrt(delta))/2/a2;
	t2 = (-a1-sqrt(delta))/2/a2;
	bool b1, b2;
	b1 = t1>u && DotProduct(v02,v12)<DotProduct(v02,v02)*(t1-u)/t1
		&& DotProduct(v02,v12)*(t1-u)/t1<DotProduct(v12,v12);
	b2 = t2>u && DotProduct(v02,v12)<DotProduct(v02,v02)*(t2-u)/t2
		     && DotProduct(v02,v12)*(t2-u)/t2<DotProduct(v12,v12);
	if(!b1 && !b2)
	{
        *u2 = MIN(*u2, MIN(F*sqrt(DotProduct(v02,v02))+u0, F*sqrt(DotProduct(v12,v12))+u1));
		return;
	}
	if(b1)
        *u2 = MIN(*u2, t1+u0);
	if(b2)
        *u2 = MIN(*u2, t2+u0);
}

void TriangularMesh::FMM_LocalSolver(int vertex, int triangle, double F)
{
// The localsolver for "vertex" in "triangle".
	TRIANGLEWITHCOOR twc;
	double u0, u1, u2;

	TRIANGLE t;
	t = m_triangles[triangle];
	// to construct twc and u0, u1 for FMM_LocalSolver_Acute().
	twc.ver2 = GetVertex3DFromPoint3D(m_vertices[vertex]);
	u2 = INFINITY;
	if(vertex==t.ver0)
	{
		if(m_vertices[t.ver1].contour>m_vertices[t.ver2].contour)
		{
			twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
			u1 = m_vertices[t.ver1].contour;
			twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
			u0 = m_vertices[t.ver2].contour;
		}
		else
		{
			twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
			u1 = m_vertices[t.ver2].contour;
			twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
			u0 = m_vertices[t.ver1].contour;
		}
	}
	if(vertex==t.ver1)
	{
		if(m_vertices[t.ver0].contour>m_vertices[t.ver2].contour)
		{
			twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
			u1 = m_vertices[t.ver0].contour;
			twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
			u0 = m_vertices[t.ver2].contour;
		}
		else
		{
			twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver2]);
			u1 = m_vertices[t.ver2].contour;
			twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
			u0 = m_vertices[t.ver0].contour;
		}
	}
	if(vertex==t.ver2)
	{
		if(m_vertices[t.ver0].contour>m_vertices[t.ver1].contour)
		{
			twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
			u1 = m_vertices[t.ver0].contour;
			twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
			u0 = m_vertices[t.ver1].contour;
		}
		else
		{
			twc.ver1 = GetVertex3DFromPoint3D(m_vertices[t.ver1]);
			u1 = m_vertices[t.ver1].contour;
			twc.ver0 = GetVertex3DFromPoint3D(m_vertices[t.ver0]);
			u0 = m_vertices[t.ver0].contour;
		}
	}
	VECTOR3D v02, v12;
	v02 = GetVector3Dfrom2Vertices(twc.ver2,twc.ver0);
	v12 = GetVector3Dfrom2Vertices(twc.ver2,twc.ver1);
	if(DotProduct(v02,v12)>=0)
		FMM_LocalSolver_Acute(twc,u0,u1,&u2,F);
	// else, recursively unfolding the neighbour triangles.
}

void TriangularMesh::Unfolding(TRIANGLEWITHCOOR s, TRIANGLEWITHCOOR *t)
{
/* to unfold triangle t to the underlying plane of triangle s.
 * We assume that s.ver1 = t.ver1 and s.ver2 = t.ver2.
 * Then we want to find new t.ver0 on the underlying plane of s.
 */ 
	VERTEX3D os, ot;
	ot = GetFootInATriangle(*t,0);
	os = GetFootInATriangle(s,0);
	VECTOR3D v;
	v = GetVector3Dfrom2Vertices(s.ver0,os);
	NormalizeVector3D(&v);
    double len;
	len = LengthVector3D(GetVector3Dfrom2Vertices(ot,(*t).ver0));

	Translate(&ot,v,len);
    (*t).ver0 = ot;
}

VERTEX3D TriangularMesh::GetFootInATriangle(TRIANGLEWITHCOOR twc, int ver)
{
/* to compute the foot of vertex ver on its opposite edge in triangle twc.
 *  ver == 0: twc.ver0;
 *  ver == 1: twc.ver1;
 *  ver == 2: twc.ver2.
 */
	VERTEX3D o;
	VECTOR3D v1,v2;
	double alpha;
	if(ver<=0)
	{// twc.ver0;
		v1 = GetVector3Dfrom2Vertices(twc.ver1,twc.ver2);
	    v2 = GetVector3Dfrom2Vertices(twc.ver1,twc.ver0);
        alpha = DotProduct(v1,v2)/DotProduct(v1,v1);
        o.x = alpha*twc.ver2.x+(1-alpha)*twc.ver1.x;
	    o.y = alpha*twc.ver2.y+(1-alpha)*twc.ver1.y;
	    o.z = alpha*twc.ver2.z+(1-alpha)*twc.ver1.z;
	}
	else if(ver==1)
	{// twc.ver1;
		v1 = GetVector3Dfrom2Vertices(twc.ver0,twc.ver2);
	    v2 = GetVector3Dfrom2Vertices(twc.ver0,twc.ver1);
        alpha = DotProduct(v1,v2)/DotProduct(v1,v1);
        o.x = alpha*twc.ver2.x+(1-alpha)*twc.ver0.x;
	    o.y = alpha*twc.ver2.y+(1-alpha)*twc.ver0.y;
	    o.z = alpha*twc.ver2.z+(1-alpha)*twc.ver0.z;
	}
	else
	{// twc.ver2;
		v1 = GetVector3Dfrom2Vertices(twc.ver0,twc.ver1);
	    v2 = GetVector3Dfrom2Vertices(twc.ver0,twc.ver2);
        alpha = DotProduct(v1,v2)/DotProduct(v1,v1);
        o.x = alpha*twc.ver1.x+(1-alpha)*twc.ver0.x;
	    o.y = alpha*twc.ver1.y+(1-alpha)*twc.ver0.y;
	    o.z = alpha*twc.ver1.z+(1-alpha)*twc.ver0.z;
	}
	return o;
}

double TriangularMesh::LengthVector3D(VECTOR3D v)
{
	return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

void TriangularMesh::Translate(VERTEX3D *v, VECTOR3D d, double len)
{
	(*v).x+= len*d.x;
	(*v).y+= len*d.y;
	(*v).z+= len*d.z;
}

void TriangularMesh::SetMarchingSteps(int ms)
{
	m_nTime = ms;
}

void TriangularMesh::Design_e1e2(TRIANGLE t, VECTOR3D *e1, VECTOR3D *e2)
{// to design the orthonormal constant vector fields e1 and e2 of the underlying plane of the triangle t.

	VECTOR3D trinormal, tmpe1, tmpe2;
	float alpha;
    //since the fields to be designed are in the tangent space, we firstly calculate its normal vector.
	tmpe1.x = m_vertices[t.ver1].x-m_vertices[t.ver0].x;
	tmpe1.y = m_vertices[t.ver1].y-m_vertices[t.ver0].y;
	tmpe1.z = m_vertices[t.ver1].z-m_vertices[t.ver0].z;
	tmpe2.x = m_vertices[t.ver2].x-m_vertices[t.ver0].x;
	tmpe2.y = m_vertices[t.ver2].y-m_vertices[t.ver0].y;
	tmpe2.z = m_vertices[t.ver2].z-m_vertices[t.ver0].z;
	trinormal = CrossProduct(tmpe1,tmpe2);
	NormalizeVector3D(&trinormal);

	tmpe1.x=0; tmpe1.y=0; tmpe1.z=1;// the z+ direction.
	alpha = DotProduct(tmpe1,trinormal);
	tmpe1.x =tmpe1.x-alpha*trinormal.x;
	tmpe1.y =tmpe1.y-alpha*trinormal.y;
	tmpe1.z =tmpe1.z-alpha*trinormal.z;
	NormalizeVector3D(&tmpe1);//
	tmpe2 = CrossProduct(trinormal,tmpe1);
	NormalizeVector3D(&tmpe2);//
	*e1 = tmpe1;
	*e2 = tmpe2;
}

void TriangularMesh::ResetContourTimestep()
{
// to reset time step for geodesic curvature flow using the weighted curve length information.
    double wcl1,wcl2,wcl3;
	vector<double>::iterator dite;

    dite=m_WeightedCurveLength.end();
	dite--;
	wcl1 = * (dite);
	dite--;
	wcl2 = * (dite);
	dite--;
	wcl3 = * (dite);

    if(BackwardsContour())
	{
		m_tStep /= (1+(wcl1-wcl2)/(wcl3-wcl2)*2);
		return;
	}
	else
		m_tStep *= 1+(wcl2-wcl1)/(wcl3-wcl2)/2;
}

bool TriangularMesh::BackwardsContour()
{
// to make a one-step backwards for the contour data if necessory.
	double wcl1,wcl2;
	vector<double>::iterator dite;
	dite=m_WeightedCurveLength.end();
	dite--;
	wcl1 = * (dite);
	dite--;
	wcl2 = * (dite);
	if(wcl1>wcl2)
	{
		for(int i=0;i<m_vnum;i++)
			m_vertices[i].contour = m_tempContour[i];
		m_WeightedCurveLength.pop_back();
		m_tSteps.pop_back();

		return true;
	}
	return false;
}

void TriangularMesh::WriteWeightedCurveLength_Forobservation()
{
	int i=0;
	vector<double>::iterator dite1, dite2;
	ofstream ofp("WeightedCurveLength_Timesteps__Forobservation.txt");
	for(dite1=m_WeightedCurveLength_Forobservation.begin(),dite2=m_tSteps_Forobservation.begin();
	dite1<m_WeightedCurveLength_Forobservation.end(),dite2<m_tSteps_Forobservation.end();dite1++,dite2++)
	{
		ofp<<i<<":   "<<*dite1<<";   "<<*dite2<<"\n";
		i++;
	}
	ofp.close();
}

void TriangularMesh::ResetContourTimestep_2()
{
	double wcl1,wcl2,wcl3;
	double dt2,dt3;
	vector<double>::iterator dite,ditedt;

    dite=m_WeightedCurveLength.end();
	dite--;
	wcl1 = * (dite);
	dite--;
	wcl2 = * (dite);
	dite--;
	wcl3 = * (dite);
    ditedt = m_tSteps.end();
	ditedt--;
	dt2 = *(ditedt);
	ditedt--;
	dt3 = *(ditedt);

    if(BackwardsContour())
	{
		m_tStep /= 1+(wcl1-wcl2)/(dt2)/((wcl3-wcl2)/dt3)*2;
	}
	else
		m_tStep *= 1+(wcl2-wcl1)/(dt2)/((wcl3-wcl2)/dt3)/2;
}

void TriangularMesh::ResetContourTimestep_3()
{
	double wcl1,wcl2,wcl3;
	double dt2,dt3;
	vector<double>::iterator dite,ditedt;

    dite=m_WeightedCurveLength.end();
	dite--;
	wcl1 = * (dite);
	dite--;
	wcl2 = * (dite);
	dite--;
	wcl3 = * (dite);
    ditedt = m_tSteps.end();
	ditedt--;
	dt2 = *(ditedt);
	ditedt--;
	dt3 = *(ditedt);

    if(BackwardsContour())
	{
		m_tStep = 0.5*(dt2+dt3*(wcl1-wcl2)/(wcl2-wcl3));
		if(m_tStep<0)
			m_tStep = dt2/2;
	}
	else
	{
		double dwcldt1, dwcldt2;
		dwcldt1 = (wcl1-wcl2)/dt2;
		dwcldt2 = (wcl2-wcl3)/dt3;
		if(dwcldt1<=dwcldt2)
    		m_tStep = dt2*((wcl2-wcl1)/(dt2)/((wcl3-wcl2)/dt3));
		else
			m_tStep = dt2*(0-dwcldt1)/(dwcldt1-dwcldt2)/2;
	}
}

void TriangularMesh::SetInitialContourFunction()
{
	int i;
	double smin;

	for(i=0;i<m_vnum;i++)
	{
		m_vertices[i].contour = 1.0;
//		m_vertices[i].contour = m_vertices[i].x-m_vertices[i].z+1.5*m_vertices[i].y+0.03;// for bunny.
//      m_vertices[i].contour = 0.08*m_vertices[i].x+m_vertices[i].z-0.16;// for horse.
//      m_vertices[i].contour = 0.25*m_vertices[i].x+m_vertices[i].y-0.16;
//      m_vertices[i].contour = m_vertices[i].x+0.4*m_vertices[i].z-0.1*m_vertices[i].y-0.13;
//      m_vertices[i].contour = POWER(m_vertices[i].x)+3*POWER(m_vertices[i].y)-0.46*0.46;//-(m_vertices[i].z+0.125);
//        m_vertices[i].contour = m_vertices[i].x+m_vertices[i].y;
/*
        smin = SignMin2(m_vertices[i].z-0.45, m_vertices[i].z+0.45);
		m_vertices[i].contour = smin;*/
	}
}

double TriangularMesh::SignMin2(double a, double b)
{
	return SIGN(a)*SIGN(b)*MIN(fabs(a),fabs(b));
}

void TriangularMesh::Contour2Intensity(double imin, double imax)
{
	double cmin, cmax;
	int i;
	cmin = 10000; cmax = -10000;
	for(i=0;i<m_vnum;i++)
	{
		if(m_vertices[i].contour<cmin)
			cmin = m_vertices[i].contour;
		if(m_vertices[i].contour>cmax)
			cmax = m_vertices[i].contour;
	}
	if(cmin==cmax)
		cmax = cmin + 1;
	for(i=0;i<m_vnum;i++)
	{
		m_vertices[i].r = (m_vertices[i].contour-cmin)/(cmax-cmin)*(imax-imin)+imin;
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}
}

double TriangularMesh::gHeatCoeff(double s)
{
	//return 	1/sqrt(s*s+m_alpha*m_alpha);// regularized TV Norm.
	//return 	1/(s*s+0.00001);// regularized BFB diffusion.
	//return 	1/(1+(s/1)*(s/1));// Perona-Malik model.
	return 	1/(s*s*s+0.001);
}

double TriangularMesh::g1HeatCoeff(double s)
{
	return 	1/sqrt(s*s+m_alpha*m_alpha);// regularized TV Norm.
	//return 	1/(s*s+0.00001);// regularized BFB diffusion.
	//return 0;
}

double TriangularMesh::g2HeatCoeff(double s)
{
	//return 	1/sqrt(s*s+m_alpha*m_alpha);// regularized TV Norm.
	//return 	1/(s*s+0.00001);// regularized BFB diffusion.
	return 0.1;
}

double TriangularMesh::GetSmallestAngleRatio()
{
/*
 *******************************************************************
 * to compute the smallest angle ratio of this mesh surface.
 * the ratio is calculated triangle by triangle.
 * in each triangle, the ratio is defined as the ratio of
 * the smallest angle and the largest angle.
 * this quantity indicates the quality of the triangular mesh.
 *******************************************************************
 */
	double ratio;
	ratio = 1; //initialization.
	int i;
	TRIANGLE t;
	POINT3d pA,pB,pC;
	VECTOR3D v1,v2;
	double angleA,angleB,angleC;
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
        // the 3 vertices of the triangle ABC.
		pA = m_vertices[t.ver0];
		pB = m_vertices[t.ver1];
		pC = m_vertices[t.ver2];
        // to compute the angles of ABC one by one.
        v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pB));
		v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pC));
        angleA = acos(DotProduct(v1,v2)/Norm(v1)/Norm(v2));
		v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pB),GetVertex3DFromPoint3D(pA));
		v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pB),GetVertex3DFromPoint3D(pC));
        angleB = acos(DotProduct(v1,v2)/Norm(v1)/Norm(v2));
		v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pC),GetVertex3DFromPoint3D(pA));
		v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pC),GetVertex3DFromPoint3D(pB));
        angleC = acos(DotProduct(v1,v2)/Norm(v1)/Norm(v2));
		// to compute the ratio of smallest and largest angle; and then update the value.
		ratio = MIN(ratio,MIN3(angleA,angleB,angleC)/MAX3(angleA,angleB,angleC));
	}

    return ratio;
}

double TriangularMesh::GetSmallestEdgeLengthRatio()
{
/*
 *******************************************************************
 * to compute the smallest edge length ratio of this mesh surface.
 * the ratio is calculated triangle by triangle.
 * in each triangle, the ratio is defined as the ratio of
 * the shortest edge length and the longest edge length.
 * this quantity indicates the quality of the triangular mesh.
 *******************************************************************
 */
	double ratio;
	ratio = 1; //initialization.
	int i;
	TRIANGLE t;
	POINT3d pA,pB,pC;
	VECTOR3D eA,eB,eC;
	double lengtheA,lengtheB,lengtheC;
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
        // the 3 vertices of the triangle ABC.
		pA = m_vertices[t.ver0];
		pB = m_vertices[t.ver1];
		pC = m_vertices[t.ver2];
        // to compute the lengths of edges one by one.
        eC = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pB));
		lengtheC = Norm(eC);
		eB = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pC));
        lengtheB = Norm(eB);
        eA = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pB),GetVertex3DFromPoint3D(pC));
        lengtheA = Norm(eA);
        // to compute the ratio between the shortest and the longest lengths; and then update the value.
		ratio = MIN(ratio,MIN3(lengtheA,lengtheB,lengtheC)/MAX3(lengtheA,lengtheB,lengtheC));
	}

    return ratio;
}

double TriangularMesh::GetSmallestEdgeLength()
{
/*
 *******************************************************************
 * to compute the smallest edge length of this mesh surface.
 * this quantity indicates the quality of the triangular mesh.
 *******************************************************************
 */
	double length;
	length = INFINITY; //initialization.
	int i;
	TRIANGLE t;
	POINT3d pA,pB,pC;
	VECTOR3D eA,eB,eC;
	double lengtheA,lengtheB,lengtheC;
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
        // the 3 vertices of the triangle ABC.
		pA = m_vertices[t.ver0];
		pB = m_vertices[t.ver1];
		pC = m_vertices[t.ver2];
        // to compute the lengths of edges one by one.
        eC = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pB));
		lengtheC = Norm(eC);
		eB = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pC));
        lengtheB = Norm(eB);
        eA = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pB),GetVertex3DFromPoint3D(pC));
        lengtheA = Norm(eA);
        // to compute the ratio between the shortest and the longest lengths; and then update the value.
		length = MIN(length,MIN3(lengtheA,lengtheB,lengtheC));
	}

    return length;
}

double TriangularMesh::GetSmallestLargestTriAreaRatio()
{
/*
 ***************************************************************
 * to compute the ratio of areas of the smallest and largest 
 * triangles. this is a global mesh quality measure.
 ***************************************************************
 */
	double ratio;
	int i;
	TRIANGLE t;
	POINT3d pA,pB,pC;
	VECTOR3D eA,eB,eC;
	double lengtheA,lengtheB,lengtheC;
	double semipre;
    double smallarea,largearea,area;
	smallarea = INFINITY;
	largearea = -1;
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
        // the 3 vertices of the triangle ABC.
		pA = m_vertices[t.ver0];
		pB = m_vertices[t.ver1];
		pC = m_vertices[t.ver2];
        // to compute the lengths of edges one by one.
        eC = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pB));
		lengtheC = Norm(eC);
		eB = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pC));
        lengtheB = Norm(eB);
        eA = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pB),GetVertex3DFromPoint3D(pC));
        lengtheA = Norm(eA);
        // to compute the area; and then update the smallest and largest area.
        semipre = (lengtheA+lengtheB+lengtheC)/2;
		area = sqrt(semipre*(semipre-lengtheA)*(semipre-lengtheB)*(semipre-lengtheC));
		smallarea = MIN(smallarea,area);
		largearea = MAX(largearea,area);
	}
	ratio = smallarea/largearea;

	return ratio;
}

double TriangularMesh::GetLargestEdgeLength()
{
/*
 *******************************************************************
 * to compute the largest edge length of this mesh surface.
 * this quantity indicates the quality of the triangular mesh.
 *******************************************************************
 */
	double length;
	length = 0; //initialization.
	int i;
	TRIANGLE t;
	POINT3d pA,pB,pC;
	VECTOR3D eA,eB,eC;
	double lengtheA,lengtheB,lengtheC;
	for(i=0;i<m_trinum;i++)
	{
		t = m_triangles[i];
        // the 3 vertices of the triangle ABC.
		pA = m_vertices[t.ver0];
		pB = m_vertices[t.ver1];
		pC = m_vertices[t.ver2];
        // to compute the lengths of edges one by one.
        eC = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pB));
		lengtheC = Norm(eC);
		eB = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pA),GetVertex3DFromPoint3D(pC));
        lengtheB = Norm(eB);
        eA = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(pB),GetVertex3DFromPoint3D(pC));
        lengtheA = Norm(eA);
        // to compute the longest edge of this triangle; and then update the longest edge of the whole mesh.
		length = MAX(length,MAX3(lengtheA,lengtheB,lengtheC));
	}

    return length;
}

void TriangularMesh::GetMeshInformation()
{
/*
 *  to compute the parameters describing the quality of the mesh.	
 *  The 1st and 2nd parameters indicate the quality of triangles, which are two local descriptors.
 *  Basically they are equivalent to each other.
 *  The 3rd parameter indicate the uniformility of the triangulation.
 *  The latter two parameters describe the mesh size.
 */	
	m_SmallestAngleRatio = GetSmallestAngleRatio(); // min_{tri} (min angle of tri / max angle of tri).
	m_SmallestEdgeLengthRatio = GetSmallestEdgeLengthRatio(); // min_{tri} (shorest edge of tri / longest edge of tri).
	m_SmallestLargestTriAreaRatio = GetSmallestLargestTriAreaRatio(); // min (area of all tri) / max (area of all tri).
	m_SmallestEdgeLength = GetSmallestEdgeLength(); // minimal length of all the edges of the mesh.
	m_LargestEdgeLength = GetLargestEdgeLength(); // maximal length of all the edges of the mesh.
	m_SmallestLargestEdgeLengthRatio = m_SmallestEdgeLength/m_LargestEdgeLength;

	ofstream ofp("MeshInformation_.txt");
	ofp<<"Mesh information:\n"
		<<"ver_number: "<<m_vnum<<"; tri_number: "<<m_trinum
		<<"; longest_edge: "<<m_LargestEdgeLength<<"; shortest_edge: "
		<<m_SmallestEdgeLength<<"; smallest_largest_area_ratio: "
		<<m_SmallestLargestTriAreaRatio<<"; smallest_edgelengthratio: "
		<<m_SmallestEdgeLengthRatio<<"; smallest_angleratio: "
		<<m_SmallestAngleRatio<<"\n";
	ofp.close();
}


void TriangularMesh::ALM_TVL2Denoising()
{
/*    
	Solves TVL2 denoising problem using ALM.

    The problem:
	    min_u regParam\sum|\nabla u|s_tau + 1/2 \|u-f\|_2.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   outTole  : the tolerance for outer iteration;
 */
	
    unsigned long nver,ntri,i;
	double *u,*uold,*b;
	VECTOR3D *p,*lambda,*w;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[nver];
	b = new double[nver];
	p = new VECTOR3D[ntri];
	lambda = new VECTOR3D[ntri];
	w = new VECTOR3D[ntri];
	uold = new double[nver];

    unsigned int innerL,l;
	double regParam, fidParam, penParam, outTole, stoppingCond;
	double oldObjFunVal,newObjFunVal;
/*   Parameters setting  */
	innerL = 1;
	regParam = 1;
	fidParam = 5000;
	penParam = 0.01;
	outTole = 0.0000001;

	numc::RowMatSym<double> alphaMinusrLap(nver,nver);
	ONERING or;
	double diag;
	unsigned long j;
	for(i=0;i<nver;i++)
	{
		or = m_ONERINGS[i];
		diag = 0;
		for(j=0;j<or.size();j++)
		{
			alphaMinusrLap(i,or[j].iver) = -penParam*or[j].co;
			diag-= or[j].co;
		}
		alphaMinusrLap(i,i) = fidParam*m_vertices[i].BCDArea - penParam*diag;
	}
    numc::SparseSolver solver;
    solver.getMatA() = alphaMinusrLap;
    solver.getMatA().mMtype = numc::CSRMatrix<double>::RealSymmIndef;
    solver.init();

/*   Initialization   */
	for(i=0;i<nver;i++)
		u[i] = 0;
	for(i=0;i<ntri;i++)
	{
		p[i].x=0; lambda[i].x=0;
		p[i].y=0; lambda[i].y=0;
		p[i].z=0; lambda[i].z=0;
	}
/*	newObjFunVal = 0;
	for(i=0;i<ntri;i++)
	{
		newObjFunVal+= Norm(p[i])*m_trianglesArea[i];
	}
	for(i=0;i<nver;i++)
	{
		newObjFunVal+= fidParam/2*POWER(u[i]-m_tempRGB[i].r)*m_vertices[i].BCDArea;
	}
*/
/*   Iteration        */
	do
	{
//		oldObjFunVal = newObjFunVal;
		for(i=0;i<nver;i++)
			uold[i] = u[i];

		/*  Inner iteration             */
		for(l=0;l<innerL;l++)
		{
		/*  Solve the u-sub problem     */
		   for(i=0;i<ntri;i++)
		   {
		    	p[i].x = lambda[i].x+penParam*p[i].x;
			    p[i].y = lambda[i].y+penParam*p[i].y;
			    p[i].z = lambda[i].z+penParam*p[i].z;
		   }
		   GetDivergence(ntri,p);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].r-m_vertices[i].divergence);
		   solver.solve(b,u);
//         ALM_TVL2Denoising_usub_nr_linbcg(nver,b,u,fidParam,penParam,itol,tol,itmax,&iter,&err);
		   /*  Projection of u             /
		   for(i=0;i<nver;i++)
		   {
			   if(u[i]<0) u[i] = 0;
			   if(u[i]>1) u[i] = 1;
		   }*/
        /*  Solve the p-sub problem     */
           GetPPIGradients(nver,u);
		   for(i=0;i<ntri;i++)
		   {
		    	w[i].x = m_triangles[i].grad.x - lambda[i].x/penParam;
		    	w[i].y = m_triangles[i].grad.y - lambda[i].y/penParam;
			    w[i].z = m_triangles[i].grad.z - lambda[i].z/penParam;
		   }
		   for(i=0;i<ntri;i++)
		   {
			    if(Norm(w[i])<=regParam/penParam)
				{
				     p[i].x = 0;
				     p[i].y = 0;
				     p[i].z = 0;
				}
			    else
				{
				     p[i].x = (1-regParam/(penParam*Norm(w[i])))*w[i].x;
				     p[i].y = (1-regParam/(penParam*Norm(w[i])))*w[i].y;
				     p[i].z = (1-regParam/(penParam*Norm(w[i])))*w[i].z;
				}
		   } 
		}
        /*   Update Lagrange multipliers        */
		//GetPPIGradients(nver,u);
		for(i=0;i<ntri;i++)
		{
            lambda[i].x+= penParam*(p[i].x-m_triangles[i].grad.x);
			lambda[i].y+= penParam*(p[i].y-m_triangles[i].grad.y);
			lambda[i].z+= penParam*(p[i].z-m_triangles[i].grad.z);
		}
		/*   Compute the stopping condition     */
        stoppingCond = 0;
		VECTOR3D temp;
/*		newObjFunVal = 0;
	    for(i=0;i<ntri;i++)
		{
		    newObjFunVal+= Norm(p[i])*m_trianglesArea[i];
		}
	    for(i=0;i<nver;i++)
		{
		    newObjFunVal+= fidParam/2*POWER(u[i]-m_tempRGB[i].r)*m_vertices[i].BCDArea;
		}

		for(i=0;i<ntri;i++)
		{
			temp.x = p[i].x-m_triangles[i].grad.x;
			temp.y = p[i].y-m_triangles[i].grad.y;
			temp.z = p[i].z-m_triangles[i].grad.z;
			stoppingCond+= m_trianglesArea[i]*(POWER(temp.x)+POWER(temp.y)+POWER(temp.z));
		}
*/
//		stoppingCond+= fabs(newObjFunVal-oldObjFunVal);
		for(i=0;i<nver;i++)
			stoppingCond+= (u[i]-uold[i]) * (u[i]-uold[i]) * m_vertices[i].BCDArea;
	}
	while(stoppingCond>outTole);

	for(i=0;i<nver;i++)
		m_vertices[i].r = u[i];
	for(i=0;i<nver;i++)
	{
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}

//    WriteData();

	delete u;
	delete b;
	delete p;
	delete w;
	delete lambda;
	delete uold;
}

void TriangularMesh::ALM_TVgL1()
{
/*    
	Solves TVgL1 problem using ALM.

    The problem:
	    min_u \sum g|\nabla u|s_tau + fidParam \|u-f\|.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   outTole  : the tolerance for outer iteration;
 */
	
    unsigned long nver,ntri,i;
	double *u,*uold,*b;
	VECTOR3D *p,*lambda_p,*w_p;
	double *z,*lambda_z,*w_z;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	u = new double[nver];
	b = new double[nver];
	p = new VECTOR3D[ntri];
	lambda_p = new VECTOR3D[ntri];
	w_p = new VECTOR3D[ntri];
	z = new double[nver];
	lambda_z = new double[nver];
	w_z = new double[nver];
	uold = new double[nver];

    unsigned int innerL,l;
	double fidParam, penParam_p, penParam_z, outTole, stoppingCond;
/*   Parameters setting  */
	innerL = 1;
	fidParam = 100;
	penParam_p = 0.0001;
	penParam_z = 1000;
	outTole = 0.00000001;

/*   Initialization   */
	for(i=0;i<nver;i++)
		u[i] = 0;
	for(i=0;i<ntri;i++)
	{
		p[i].x=0; lambda_p[i].x=0;
		p[i].y=0; lambda_p[i].y=0;
		p[i].z=0; lambda_p[i].z=0;
	}
	for(i=0;i<nver;i++)
	{
		z[i] = 0; lambda_z[i] = 0;
	}

/*   Iteration        */
	do
	{
		for(i=0;i<nver;i++)
			uold[i] = u[i];

		/*  Inner iteration             */
		for(l=0;l<innerL;l++)
		{
		/*  Solve the u-sub problem     */
		   for(i=0;i<ntri;i++)
		   {
		    	p[i].x = lambda_p[i].x+penParam_p*p[i].x;
			    p[i].y = lambda_p[i].y+penParam_p*p[i].y;
			    p[i].z = lambda_p[i].z+penParam_p*p[i].z;
		   }
		   GetDivergence(ntri,p);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(penParam_z*m_tempRGB[i].r+lambda_z[i]+penParam_z*z[i]-m_vertices[i].divergence);
           ALM_TVgL1_usub_nr_linbcg(nver,b,u,penParam_p,penParam_z,itol,tol,itmax,&iter,&err);
        /*  Solve the p-sub problem     */
           GetPPIGradients(nver,u);
		   for(i=0;i<ntri;i++)
		   {
		    	w_p[i].x = m_triangles[i].grad.x - lambda_p[i].x/penParam_p;
		    	w_p[i].y = m_triangles[i].grad.y - lambda_p[i].y/penParam_p;
			    w_p[i].z = m_triangles[i].grad.z - lambda_p[i].z/penParam_p;
		   }
		   for(i=0;i<ntri;i++)
		   {
			    if(Norm(w_p[i])<=1/penParam_p)
				{
				     p[i].x = 0;
				     p[i].y = 0;
				     p[i].z = 0;
				}
			    else
				{
				     p[i].x = (1-1/(penParam_p*Norm(w_p[i])))*w_p[i].x;
				     p[i].y = (1-1/(penParam_p*Norm(w_p[i])))*w_p[i].y;
				     p[i].z = (1-1/(penParam_p*Norm(w_p[i])))*w_p[i].z;
				}
		   } 
       /*   Solve the z-sub problem     */
		   for(i=0;i<nver;i++)
			   w_z[i] = u[i]-m_vertices[i].r-lambda_z[i]/penParam_z;
		   for(i=0;i<nver;i++)
		   {
			   if(fabs(w_z[i])<=fidParam/penParam_z)
				   z[i] = 0;
			   else
				   z[i] = (1-fidParam/(penParam_z*fabs(w_z[i])))*w_z[i];
		   }
		}
        /*   Update Lagrange multipliers        */
		//GetPPIGradients(nver,u);
		for(i=0;i<ntri;i++)
		{
            lambda_p[i].x+= penParam_p*(p[i].x-m_triangles[i].grad.x);
			lambda_p[i].y+= penParam_p*(p[i].y-m_triangles[i].grad.y);
			lambda_p[i].z+= penParam_p*(p[i].z-m_triangles[i].grad.z);
		}
		for(i=0;i<nver;i++)
			lambda_z[i]+= penParam_z*(z[i]-u[i]+m_vertices[i].r);

		/*   Compute the stopping condition     */
        stoppingCond = 0;
		for(i=0;i<nver;i++)
			stoppingCond+= (u[i]-uold[i]) * (u[i]-uold[i]) * m_vertices[i].BCDArea;
	}
	while(stoppingCond>outTole);

	for(i=0;i<nver;i++)
		m_vertices[i].r = u[i];
	for(i=0;i<nver;i++)
	{
		m_vertices[i].g = m_vertices[i].r;
		m_vertices[i].b = m_vertices[i].r;
	}

	delete u;
	delete b;
	delete p;
	delete w_p;
	delete lambda_p;
	delete z;
	delete w_z;
	delete lambda_z;
	delete uold;
}

void TriangularMesh::ALM_MulRegionLabelling()
{
/*    
	Solves Multi-Region Labelling problem using ALM.

    The problem:
	    min_u \sum |\nabla u|s_tau + fidParam <u,s>.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   medTole  : the tolerance for median iteration;
	   outTole  : the tolerance for outer iteration;
 */

	unsigned long nRegion; // the number of regions.
	long iRegion;
    unsigned long nver,ntri;
	long i;
    double *sDev;
	double *u,*uold,*b;
	double *utemp;
	VECTOR3D *p,*lambda_p,*w_p;
	VECTOR3D *ptemp;
	double *z,*lambda_z;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;

	nRegion = 3;
	
	sDev =  new double[nRegion*nver];
	u = new double[nRegion*nver];
	b = new double[nver];
	utemp = new double[nver];
	p = new VECTOR3D[nRegion*ntri];
	ptemp = new VECTOR3D[ntri];
	lambda_p = new VECTOR3D[nRegion*ntri];
	w_p = new VECTOR3D[nRegion*ntri];
	z = new double[nRegion*nver];
	lambda_z = new double[nRegion*nver];
	uold = new double[nRegion*nver];

    unsigned int innerL,il;
	unsigned int medL,ml;
	unsigned int nStoppingJudge,iStoppingJudge;
	double regParam, penParam_p, penParam_z, outTole, stoppingCond;
	double normMVec,normuVec;
/*   Parameters setting  */
	innerL = 1;
	medL = 50;// 10 for sphere.
	nStoppingJudge = 5;
	iStoppingJudge = 0;
	regParam = 0.0003;
	penParam_p = regParam*1.e-3;
	penParam_z = penParam_p*1.e6;
	outTole = 0.0001;
	stoppingCond = 100000;

//  for testing the algorithm .............
	int test_nIter=0;
	time_t test_begin, test_end;
	vecDoubles test_cpu_u;
	vecDoubles test_cpu_p;
	vecDoubles test_cpu_z;
	vecDoubles test_energy_TV;
	vecDoubles test_energy_fid;
	double energy_TV, energy_fid;
    test_cpu_u.clear();
	test_cpu_p.clear();
	test_cpu_z.clear();
	test_energy_TV.clear();
	test_energy_fid.clear();

	numc::RowMatSym<double> rzMinusrpLap(nver,nver);
	ONERING or;
	double diag;
	unsigned long j;
	for(i=0;i<nver;i++)
	{
		or = m_ONERINGS[i];
		diag = 0;
		for(j=0;j<or.size();j++)
		{
			rzMinusrpLap(i,or[j].iver) = -penParam_p*or[j].co;
			diag-= or[j].co;
		}
		rzMinusrpLap(i,i) = penParam_z*m_vertices[i].BCDArea - penParam_p*diag;
	}
    numc::SparseSolver solver;
    solver.getMatA() = rzMinusrpLap;
    solver.getMatA().mMtype = numc::CSRMatrix<double>::RealSymmIndef;
    solver.init();

/*   Initialization   */
	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
	    for(i=0;i<nver;i++)
		    u[iRegion*nver+i] = double(rand())/RAND_MAX;
	}
    ALM_MulRegionLabelling_Project2K(u,nRegion,nver);

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
	    for(i=0;i<ntri;i++)
		{
		    p[iRegion*ntri+i].x=0; lambda_p[iRegion*ntri+i].x=0;
		    p[iRegion*ntri+i].y=0; lambda_p[iRegion*ntri+i].y=0;
		    p[iRegion*ntri+i].z=0; lambda_p[iRegion*ntri+i].z=0;
		}
	}

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
	    for(i=0;i<nver;i++)
		{
		    z[iRegion*nver+i] = u[iRegion*nver+i]; lambda_z[iRegion*nver+i] = 0;
		}
	}
//  One more outer iteration ............   //
/*   Iteration        */
	do
	{
        //   for testing ................   //
		test_nIter++;

		for(i=0;i<nRegion*nver;i++)
			uold[i] = u[i];
		iStoppingJudge++;
		/*  update the s  */
//		ALM_MulRegionLabelling_sDev_grayimage(z,sDev,nRegion,nver);
		ALM_MulRegionLabelling_sDev_colorimage(z,sDev,nRegion,nver);

		/*   for new s, minimize the functional w.r.t. u  using ALM   */
        /*   median iteration, ALM for given s   */
	 	for(ml=0;ml<medL;ml++)
		{

    		/*  Inner iteration, ALM with given multipliers   */
		    for(il=0;il<innerL;il++)
			{
		        /*  Solve the u-sub problem     */
				test_begin = clock();

		        for(iRegion=0;iRegion<nRegion;iRegion++)
				{
				    for(i=0;i<ntri;i++)
					{
		    	        ptemp[i].x = lambda_p[iRegion*ntri+i].x+penParam_p*p[iRegion*ntri+i].x;
			            ptemp[i].y = lambda_p[iRegion*ntri+i].y+penParam_p*p[iRegion*ntri+i].y;
			            ptemp[i].z = lambda_p[iRegion*ntri+i].z+penParam_p*p[iRegion*ntri+i].z;
					}
		            GetDivergence(ntri,ptemp);
		            for(i=0;i<nver;i++)
		    	        b[i] = m_vertices[i].BCDArea*(lambda_z[iRegion*nver+i]+penParam_z*z[iRegion*nver+i]-m_vertices[i].divergence);
                    solver.solve(b, utemp);
//                  ALM_MulRegionLabelling_usub_nr_linbcg(nver,b,utemp,penParam_p,penParam_z,itol,tol,itmax,&iter,&err);
					GetPPIGradients(nver,utemp);
		            for(i=0;i<ntri;i++)
					{
		    	        w_p[iRegion*ntri+i].x = m_triangles[i].grad.x - lambda_p[iRegion*ntri+i].x/penParam_p;
		    	        w_p[iRegion*ntri+i].y = m_triangles[i].grad.y - lambda_p[iRegion*ntri+i].y/penParam_p;
			            w_p[iRegion*ntri+i].z = m_triangles[i].grad.z - lambda_p[iRegion*ntri+i].z/penParam_p;
					}
					for(i=0;i<nver;i++)
						u[iRegion*nver+i] = utemp[i];
				}

				test_end = clock();
				test_cpu_u.push_back(double(test_end - test_begin)/CLOCKS_PER_SEC);
                /*  Solve the p-sub problem     */
                test_begin = clock();

		        for(i=0;i<ntri;i++)
				{
					normMVec = 0;
					for(iRegion=0;iRegion<nRegion;iRegion++)
						//normMVec+= POWER(Norm(w_p[iRegion*ntri+i]));
						normMVec+= POWER(w_p[iRegion*ntri+i].x)+POWER(w_p[iRegion*ntri+i].y)+POWER(w_p[iRegion*ntri+i].z);
					normMVec = sqrt(normMVec);
			        if(normMVec<=regParam/penParam_p)
					{
				        for(iRegion=0;iRegion<nRegion;iRegion++)
						{
				             p[iRegion*ntri+i].x = 0;
				             p[iRegion*ntri+i].y = 0;
				             p[iRegion*ntri+i].z = 0;
						}
					}
			        else
					{
						for(iRegion=0;iRegion<nRegion;iRegion++)
						{
				            p[iRegion*ntri+i].x = (1-regParam/(penParam_p*normMVec))*w_p[iRegion*ntri+i].x;
				            p[iRegion*ntri+i].y = (1-regParam/(penParam_p*normMVec))*w_p[iRegion*ntri+i].y;
				            p[iRegion*ntri+i].z = (1-regParam/(penParam_p*normMVec))*w_p[iRegion*ntri+i].z;
						}
					}
				}
				
				test_end = clock();
				test_cpu_p.push_back(double(test_end - test_begin)/CLOCKS_PER_SEC);
                /*   Solve the z-sub problem     */
                test_begin = clock();

		        for(iRegion=0;iRegion<nRegion;iRegion++)
				{
				    for(i=0;i<nver;i++)
					{
						z[iRegion*nver+i] = u[iRegion*nver+i] - (sDev[iRegion*nver+i]+lambda_z[iRegion*nver+i])/penParam_z;
					}
				}
				ALM_MulRegionLabelling_Project2K(z,nRegion,nver);

				test_end = clock();
				test_cpu_z.push_back(double(test_end - test_begin)/CLOCKS_PER_SEC);
			}
            /*   Update Lagrange multipliers        */
			for(iRegion=0;iRegion<nRegion;iRegion++)
			{
				for(i=0;i<nver;i++)
					utemp[i] = u[iRegion*nver+i];
				GetPPIGradients(nver,utemp);
		        for(i=0;i<ntri;i++)
				{
                    lambda_p[iRegion*ntri+i].x+= penParam_p*(p[iRegion*ntri+i].x-m_triangles[i].grad.x);
			        lambda_p[iRegion*ntri+i].y+= penParam_p*(p[iRegion*ntri+i].y-m_triangles[i].grad.y);
			        lambda_p[iRegion*ntri+i].z+= penParam_p*(p[iRegion*ntri+i].z-m_triangles[i].grad.z);
				}
		        for(i=0;i<nver;i++)
			        lambda_z[iRegion*nver+i]+= penParam_z*(z[iRegion*nver+i]-u[iRegion*nver+i]);
			}
		}
		if(iStoppingJudge==nStoppingJudge)
		{		/*   Compute the stopping condition     */
            iStoppingJudge = 0;
            stoppingCond = 0;
		    for(i=0;i<nver;i++)
			{
			     normuVec = 0;
			     for(iRegion=0;iRegion<nRegion;iRegion++)
				 {
				     normuVec+= (u[iRegion*nver+i]-uold[iRegion*nver+i]) * (u[iRegion*nver+i]-uold[iRegion*nver+i]);
				 }
			     stoppingCond+= normuVec * m_vertices[i].BCDArea;
			}
		}
		//   Record the energy terms ..........   //
		energy_TV = 0;
		energy_fid = 0;
		for(i=0;i<nver;i++)
		{
			normuVec = 0;
			for(iRegion=0;iRegion<nRegion;iRegion++)
				normuVec+= z[iRegion*nver+i]*sDev[iRegion*nver+i];
			normuVec*= m_vertices[i].BCDArea;
			energy_fid+= normuVec;
		}
		test_energy_fid.push_back(energy_fid);
		for(i=0;i<ntri;i++)
		{
			normMVec = 0;
			for(iRegion=0;iRegion<nRegion;iRegion++)
				normMVec+= POWER(p[iRegion*ntri+i].x)+POWER(p[iRegion*ntri+i].y)+POWER(p[iRegion*ntri+i].z);
			normMVec = sqrt(normMVec);
			normMVec*= regParam;
			normMVec*= m_trianglesArea[i];
			energy_TV+= normMVec;
		}
		test_energy_TV.push_back(energy_TV);
	}
	while(stoppingCond>outTole);

	ALM_MulRegionLabelling_Binarization(u,nRegion,nver);
	for(i=0;i<nver;i++)
	{
		switch(m_vertices[i].RegionId)
		{
		case 0:
			m_vertices[i].r = 1;
			m_vertices[i].g = 0;
			m_vertices[i].b = 0;
			break;
		case 1:
			m_vertices[i].r = 0;
			m_vertices[i].g = 1;
			m_vertices[i].b = 0;
			break;
		case 2:
			m_vertices[i].r = 0;
			m_vertices[i].g = 0;
			m_vertices[i].b = 1;
			break;
		case 3:
			m_vertices[i].r = 1;
			m_vertices[i].g = 1;
			m_vertices[i].b = 0;
		}
	}

    //   for testing ..................   //
	ofstream ofp("test.txt");
	ofp<<" testing results: \n";
	ofp<<" outer iteration number: "<<test_nIter<<"\n";
	ofp<<" time for u-sub;    time for p-sub;    time for z-sub;\n";
	for(i=0;i<test_nIter*medL*innerL;i++)
	{
		ofp<<i<<":  "<<test_cpu_u[i]<<";   "<<test_cpu_p[i]<<";   "<<test_cpu_z[i]<<";\n";
	}
	ofp<<" fid energy  ;  TV energy;  total;\n";
	for(i=0;i<test_nIter;i++)
	{
		ofp<<i<<":  "<<test_energy_fid[i]<<";   "<<test_energy_TV[i]<<";   "<<test_energy_fid[i]+test_energy_TV[i]<<";\n";
	}
	ofp.close();

	delete u;
	delete b;
	delete p;
	delete w_p;
	delete lambda_p;
	delete z;
	delete lambda_z;
	delete uold;
	delete utemp;
	delete ptemp;
	delete sDev;
	test_cpu_u.clear();
	test_cpu_p.clear();
	test_cpu_z.clear();
	test_energy_fid.clear();
	test_energy_TV.clear();
}

void TriangularMesh::ALM_MulRegionLabelling_givenMean()
{
/*    
	Solves Multi-Region Labelling problem with given means using ALM.

    The problem:
	    min_u \sum |\nabla u|s_tau + fidParam <u,s>.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   outTole  : the tolerance for outer iteration;
 */

	unsigned long nRegion; // the number of regions.
	long iRegion;
    unsigned long nver,ntri;
	long i;
    double *sDev;
	double *u,*uold,*b;
	double *utemp;
	VECTOR3D *p,*lambda_p,*w_p;
	VECTOR3D *ptemp;
	double *z,*lambda_z;
	double *mean;
	double *meanr,*meang,*meanb;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 100;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;

	nRegion = 2;
	
	sDev =  new double[nRegion*nver];
	u = new double[nRegion*nver];
	b = new double[nver];
	utemp = new double[nver];
	p = new VECTOR3D[nRegion*ntri];
	ptemp = new VECTOR3D[ntri];
	lambda_p = new VECTOR3D[nRegion*ntri];
	w_p = new VECTOR3D[nRegion*ntri];
	z = new double[nRegion*nver];
	lambda_z = new double[nRegion*nver];
	uold = new double[nRegion*nver];

	mean = new double[nRegion];
    mean[0] = 0.1;
	mean[1] = 0.57;
//	mean[2] = 0.7;
//	mean[3] = 0.9;

//	meanr = new double[nRegion];
//	meang = new double[nRegion];
//	meanb = new double[nRegion];
//	meanr[0] = 0.9; meang[0] = 0.9; meanb[0] = 0.9;
//	meanr[1] = 0.501961; meang[1] = 0.25098; meanb[1] = 0.25098;
//	meanr[2] = 1; meang[2] = 0; meanb[2] = 0.501961;
//	meanr[3] = 0; meang[3] = 0; meanb[3] = 1;

	ALM_MulRegionLabelling_sDev_grayimage_givenMean(mean,sDev,nRegion,nver);
//	ALM_MulRegionLabelling_sDev_colorimage_givenMean(meanr,meang,meanb,sDev,nRegion,nver);

    unsigned int innerL,il;
	unsigned int nStoppingJudge,iStoppingJudge;
	double regParam, penParam_p, penParam_z, outTole, stoppingCond;
	double normMVec,normuVec;
/*   Parameters setting  */
	innerL = 1;
	nStoppingJudge = 5;
	iStoppingJudge = 0;
	regParam = 0.00003;
	penParam_p = regParam*1.e-3;
	penParam_z = penParam_p*1.e6;
	outTole = 0.00001;
	stoppingCond = 100000;

//  for testing the algorithm .............
	int test_nIter=0;
	time_t test_begin, test_end;
	double energy_TV, energy_fid;
	vecDoubles test_cpu_u;
	vecDoubles test_cpu_p;
	vecDoubles test_cpu_z;
	vecDoubles test_energy_TV;
	vecDoubles test_energy_fid;
    test_cpu_u.clear();
	test_cpu_p.clear();
	test_cpu_z.clear();
	test_energy_TV.clear();
	test_energy_fid.clear();

	numc::RowMatSym<double> rzMinusrpLap(nver,nver);
	ONERING or;
	double diag;
	unsigned long j;
	for(i=0;i<nver;i++)
	{
		or = m_ONERINGS[i];
		diag = 0;
		for(j=0;j<or.size();j++)
		{
			rzMinusrpLap(i,or[j].iver) = -penParam_p*or[j].co;
			diag-= or[j].co;
		}
		rzMinusrpLap(i,i) = penParam_z*m_vertices[i].BCDArea - penParam_p*diag;
	}
    numc::SparseSolver solver;
    solver.getMatA() = rzMinusrpLap;
    solver.getMatA().mMtype = numc::CSRMatrix<double>::RealSymmIndef;
    solver.init();

/*   Initialization   */
	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
	    for(i=0;i<nver;i++)
		    u[iRegion*nver+i] = double(rand())/RAND_MAX;
	}
    ALM_MulRegionLabelling_Project2K(u,nRegion,nver);

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
	    for(i=0;i<ntri;i++)
		{
		    p[iRegion*ntri+i].x=0; lambda_p[iRegion*ntri+i].x=0;
		    p[iRegion*ntri+i].y=0; lambda_p[iRegion*ntri+i].y=0;
		    p[iRegion*ntri+i].z=0; lambda_p[iRegion*ntri+i].z=0;
		}
	}

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
	    for(i=0;i<nver;i++)
		{
		    z[iRegion*nver+i] = u[iRegion*nver+i]; lambda_z[iRegion*nver+i] = 0;
		}
	}

/*   Iteration        */
	do
	{
        //   for testing ................   //
		test_nIter++;

		for(i=0;i<nRegion*nver;i++)
			uold[i] = u[i];
		iStoppingJudge++;

   		/*  Inner iteration, ALM with given multipliers   */
	    for(il=0;il<innerL;il++)
		{
	        /*  Solve the u-sub problem     */
			test_begin = clock();

			for(iRegion=0;iRegion<nRegion;iRegion++)
			{
			    for(i=0;i<ntri;i++)
				{
		            ptemp[i].x = lambda_p[iRegion*ntri+i].x+penParam_p*p[iRegion*ntri+i].x;
			        ptemp[i].y = lambda_p[iRegion*ntri+i].y+penParam_p*p[iRegion*ntri+i].y;
			        ptemp[i].z = lambda_p[iRegion*ntri+i].z+penParam_p*p[iRegion*ntri+i].z;
				}
		        GetDivergence(ntri,ptemp);
		        for(i=0;i<nver;i++)
		    	    b[i] = m_vertices[i].BCDArea*(lambda_z[iRegion*nver+i]+penParam_z*z[iRegion*nver+i]-m_vertices[i].divergence);
                solver.solve(b, utemp);
//              ALM_MulRegionLabelling_usub_nr_linbcg(nver,b,utemp,penParam_p,penParam_z,itol,tol,itmax,&iter,&err);
				GetPPIGradients(nver,utemp);
		        for(i=0;i<ntri;i++)
				{
		    	    w_p[iRegion*ntri+i].x = m_triangles[i].grad.x - lambda_p[iRegion*ntri+i].x/penParam_p;
		    	    w_p[iRegion*ntri+i].y = m_triangles[i].grad.y - lambda_p[iRegion*ntri+i].y/penParam_p;
			        w_p[iRegion*ntri+i].z = m_triangles[i].grad.z - lambda_p[iRegion*ntri+i].z/penParam_p;
				}
				for(i=0;i<nver;i++)
					u[iRegion*nver+i] = utemp[i];
			}

			test_end = clock();
			test_cpu_u.push_back(double(test_end - test_begin)/CLOCKS_PER_SEC);
            /*  Solve the p-sub problem     */
            test_begin = clock();

		    for(i=0;i<ntri;i++)
			{
				normMVec = 0;
				for(iRegion=0;iRegion<nRegion;iRegion++)
					//normMVec+= POWER(Norm(w_p[iRegion*ntri+i]));
					normMVec+= POWER(w_p[iRegion*ntri+i].x)+POWER(w_p[iRegion*ntri+i].y)+POWER(w_p[iRegion*ntri+i].z);
				normMVec = sqrt(normMVec);
			    if(normMVec<=regParam/penParam_p)
				{
				    for(iRegion=0;iRegion<nRegion;iRegion++)
					{
				         p[iRegion*ntri+i].x = 0;
				         p[iRegion*ntri+i].y = 0;
				         p[iRegion*ntri+i].z = 0;
					}
				}
			    else
				{
				    for(iRegion=0;iRegion<nRegion;iRegion++)
					{
				        p[iRegion*ntri+i].x = (1-regParam/(penParam_p*normMVec))*w_p[iRegion*ntri+i].x;
				        p[iRegion*ntri+i].y = (1-regParam/(penParam_p*normMVec))*w_p[iRegion*ntri+i].y;
				        p[iRegion*ntri+i].z = (1-regParam/(penParam_p*normMVec))*w_p[iRegion*ntri+i].z;
					}
				}
			}
				
			test_end = clock();
			test_cpu_p.push_back(double(test_end - test_begin)/CLOCKS_PER_SEC);
            /*   Solve the z-sub problem     */
            test_begin = clock();

	        for(iRegion=0;iRegion<nRegion;iRegion++)
			{
			    for(i=0;i<nver;i++)
				{
					z[iRegion*nver+i] = u[iRegion*nver+i] - (sDev[iRegion*nver+i]+lambda_z[iRegion*nver+i])/penParam_z;
				}
			}
			ALM_MulRegionLabelling_Project2K(z,nRegion,nver);

			test_end = clock();
			test_cpu_z.push_back(double(test_end - test_begin)/CLOCKS_PER_SEC);
		}
        /*   Update Lagrange multipliers        */
		for(iRegion=0;iRegion<nRegion;iRegion++)
		{
			for(i=0;i<nver;i++)
				utemp[i] = u[iRegion*nver+i];
			GetPPIGradients(nver,utemp);
		    for(i=0;i<ntri;i++)
			{
                lambda_p[iRegion*ntri+i].x+= penParam_p*(p[iRegion*ntri+i].x-m_triangles[i].grad.x);
			    lambda_p[iRegion*ntri+i].y+= penParam_p*(p[iRegion*ntri+i].y-m_triangles[i].grad.y);
			    lambda_p[iRegion*ntri+i].z+= penParam_p*(p[iRegion*ntri+i].z-m_triangles[i].grad.z);
			}
		    for(i=0;i<nver;i++)
			    lambda_z[iRegion*nver+i]+= penParam_z*(z[iRegion*nver+i]-u[iRegion*nver+i]);
		}

		if(iStoppingJudge==nStoppingJudge)
		{		/*   Compute the stopping condition     */
            iStoppingJudge = 0;
            stoppingCond = 0;
		    for(i=0;i<nver;i++)
			{
			     normuVec = 0;
			     for(iRegion=0;iRegion<nRegion;iRegion++)
				 {
				     normuVec+= (u[iRegion*nver+i]-uold[iRegion*nver+i]) * (u[iRegion*nver+i]-uold[iRegion*nver+i]);
				 }
			     stoppingCond+= normuVec * m_vertices[i].BCDArea;
			}
		}
				//   Record the energy terms ..........   //
		energy_TV = 0;
		energy_fid = 0;
		for(i=0;i<nver;i++)
		{
			normuVec = 0;
			for(iRegion=0;iRegion<nRegion;iRegion++)
				normuVec+= z[iRegion*nver+i]*sDev[iRegion*nver+i];
			normuVec*= m_vertices[i].BCDArea;
			energy_fid+= normuVec;
		}
		test_energy_fid.push_back(energy_fid);
		for(i=0;i<ntri;i++)
		{
			normMVec = 0;
			for(iRegion=0;iRegion<nRegion;iRegion++)
				normMVec+= POWER(p[iRegion*ntri+i].x)+POWER(p[iRegion*ntri+i].y)+POWER(p[iRegion*ntri+i].z);
			normMVec = sqrt(normMVec);
			normMVec*= regParam;
			normMVec*= m_trianglesArea[i];
			energy_TV+= normMVec;
		}
		test_energy_TV.push_back(energy_TV);
    }
	while(stoppingCond>outTole);
//    while(test_nIter<10);

	ALM_MulRegionLabelling_Binarization(u,nRegion,nver);
/*
	for(i=0;i<nver;i++)
	{
		switch(m_vertices[i].RegionId)
		{
		case 0:
			m_vertices[i].r = 1;
			m_vertices[i].g = 0;
			m_vertices[i].b = 0;
			break;
		case 1:
			m_vertices[i].r = 0;
			m_vertices[i].g = 1;
			m_vertices[i].b = 0;
			break;
		case 2:
			m_vertices[i].r = 0;
			m_vertices[i].g = 0;
			m_vertices[i].b = 1;
			break;
		case 3:
			m_vertices[i].r = 1;
			m_vertices[i].g = 1;
			m_vertices[i].b = 0;
		}

	}
*/
    //   to check the sDev ..........   //
    //   for testing ..................   //
	ofstream ofp("test.txt");
	ofp<<" testing results: \n";
	ofp<<" outer iteration number: "<<test_nIter<<"\n";
	ofp<<" time for u-sub;    time for p-sub;    time for z-sub;\n";
	for(i=0;i<test_nIter*innerL;i++)
	{
		ofp<<i<<":  "<<test_cpu_u[i]<<";   "<<test_cpu_p[i]<<";   "<<test_cpu_z[i]<<";\n";
	}
	ofp<<" fid energy  ;  TV energy;  total;\n";
	for(i=0;i<test_nIter;i++)
	{
		ofp<<i<<":  "<<test_energy_fid[i]<<";   "<<test_energy_TV[i]<<";   "<<test_energy_fid[i]+test_energy_TV[i]<<";\n";
	}
	ofp.close();

	delete u;
	delete b;
	delete p;
	delete w_p;
	delete lambda_p;
	delete z;
	delete lambda_z;
	delete uold;
	delete utemp;
	delete ptemp;
	delete sDev;
//	delete mean;
}

void TriangularMesh::ALM_vTVL2Denoising()
{
/*  Solves vectorial TV L2 denoising problem using ALM.

    The problem:
	    min_u regParam*\sum|\nabla u|s_tau + 1/2 \|u-f\|_2.

    Here u is a multi-channel image.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   outTole  : the tolerance for outer iteration;
 */
	
    unsigned long nver,ntri,i;
	double *ur,*ug,*ub,*b;
	double *urold,*ugold,*ubold;
	VECTOR3D *pr,*pg,*pb,*lambdar,*lambdag,*lambdab,wr,wg,wb,w;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 300;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	ur = new double[nver];
	ug = new double[nver];
	ub = new double[nver];
	urold = new double[nver];
	ugold = new double[nver];
	ubold = new double[nver];
	b = new double[nver];
	pr = new VECTOR3D[ntri];
	pg = new VECTOR3D[ntri];
	pb = new VECTOR3D[ntri];
	lambdar = new VECTOR3D[ntri];
	lambdag = new VECTOR3D[ntri];
	lambdab = new VECTOR3D[ntri];

    unsigned int innerL,l;
	double regParam, fidParam, penParam, outTole, stoppingCond;
/*   Parameters setting  */
	innerL = 1;
	fidParam = 1000;
	regParam = 1;
	penParam = 0.01;
	outTole = 0.0000001;

	numc::RowMatSym<double> alphaMinusrLap(nver,nver);
	ONERING or;
	double diag;
	unsigned long j;
	for(i=0;i<nver;i++)
	{
		or = m_ONERINGS[i];
		diag = 0;
		for(j=0;j<or.size();j++)
		{
			alphaMinusrLap(i,or[j].iver) = -penParam*or[j].co;
			diag-= or[j].co;
		}
		alphaMinusrLap(i,i) = fidParam*m_vertices[i].BCDArea - penParam*diag;
	}
    numc::SparseSolver solver;
    solver.getMatA() = alphaMinusrLap;
    solver.getMatA().mMtype = numc::CSRMatrix<double>::RealSymmIndef;
    solver.init();

/*   Initialization   */
	for(i=0;i<nver;i++)
	{	
		ur[i] = 0; ug[i] = 0; ub[i] = 0;
	};
	for(i=0;i<ntri;i++)
	{
		pr[i].x=0; lambdar[i].x=0;
		pr[i].y=0; lambdar[i].y=0;
		pr[i].z=0; lambdar[i].z=0;

		pg[i].x=0; lambdag[i].x=0;
		pg[i].y=0; lambdag[i].y=0;
		pg[i].z=0; lambdag[i].z=0;

		pb[i].x=0; lambdab[i].x=0;
		pb[i].y=0; lambdab[i].y=0;
		pb[i].z=0; lambdab[i].z=0;
	}
/*   Iteration        */
	do
	{
		for(i=0;i<nver;i++)
		{
			urold[i] = ur[i];
			ugold[i] = ug[i];
			ubold[i] = ub[i];
		}
		/*  Inner iteration             */
		for(l=0;l<innerL;l++)
		{
		/*  Solve the u-sub problem  (solve for ur,ug,ub).   */

		   for(i=0;i<ntri;i++)
		   {
		    	pr[i].x = lambdar[i].x+penParam*pr[i].x;
			    pr[i].y = lambdar[i].y+penParam*pr[i].y;
			    pr[i].z = lambdar[i].z+penParam*pr[i].z;
		   }
		   GetDivergence(ntri,pr);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].r-m_vertices[i].divergence);
		   solver.solve(b,ur);
//         ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,ur,fidParam,penParam,itol,tol,itmax,&iter,&err);

		   for(i=0;i<ntri;i++)
		   {
		    	pg[i].x = lambdag[i].x+penParam*pg[i].x;
			    pg[i].y = lambdag[i].y+penParam*pg[i].y;
			    pg[i].z = lambdag[i].z+penParam*pg[i].z;
		   }
		   GetDivergence(ntri,pg);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].g-m_vertices[i].divergence);
		   solver.solve(b,ug);
//         ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,ug,fidParam,penParam,itol,tol,itmax,&iter,&err);

		   for(i=0;i<ntri;i++)
		   {
		    	pb[i].x = lambdab[i].x+penParam*pb[i].x;
			    pb[i].y = lambdab[i].y+penParam*pb[i].y;
			    pb[i].z = lambdab[i].z+penParam*pb[i].z;
		   }
		   GetDivergence(ntri,pb);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].b-m_vertices[i].divergence);
		   solver.solve(b,ub);
//         ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,ub,fidParam,penParam,itol,tol,itmax,&iter,&err);

        /*  Solve the p-sub problem  (solve for pr,pg,pb).   */
           
		   GetPPIGradients(nver,ur);
		   GetPPIGradients1(nver,ug);
		   GetPPIGradients2(nver,ub);

		   for(i=0;i<ntri;i++)
		   {
			    wr.x = m_triangles[i].grad.x - lambdar[i].x/penParam;
		    	wr.y = m_triangles[i].grad.y - lambdar[i].y/penParam;
			    wr.z = m_triangles[i].grad.z - lambdar[i].z/penParam;

			    wg.x = m_triangles[i].grad1.x - lambdag[i].x/penParam;
		    	wg.y = m_triangles[i].grad1.y - lambdag[i].y/penParam;
			    wg.z = m_triangles[i].grad1.z - lambdag[i].z/penParam;

			    wb.x = m_triangles[i].grad2.x - lambdab[i].x/penParam;
		    	wb.y = m_triangles[i].grad2.y - lambdab[i].y/penParam;
			    wb.z = m_triangles[i].grad2.z - lambdab[i].z/penParam;

				if(DotProduct(wr,wr)+DotProduct(wg,wg)+DotProduct(wb,wb) <= POWER(regParam/penParam))
				{
				     pr[i].x = 0; pg[i].x = 0; pb[i].x = 0;
				     pr[i].y = 0; pg[i].y = 0; pb[i].y = 0;
				     pr[i].z = 0; pg[i].z = 0; pb[i].z = 0;
				}
			    else
				{
					 double tempNorm;
					 tempNorm = sqrt(DotProduct(wr,wr)+DotProduct(wg,wg)+DotProduct(wb,wb));

					 pr[i].x = (1-regParam/penParam/tempNorm)*wr.x; pg[i].x = (1-1/penParam/tempNorm)*wg.x; pb[i].x = (1-1/penParam/tempNorm)*wb.x;
				     pr[i].y = (1-regParam/penParam/tempNorm)*wr.y; pg[i].y = (1-1/penParam/tempNorm)*wg.y; pb[i].y = (1-1/penParam/tempNorm)*wb.y;
				     pr[i].z = (1-regParam/penParam/tempNorm)*wr.z; pg[i].z = (1-1/penParam/tempNorm)*wg.z; pb[i].z = (1-1/penParam/tempNorm)*wb.z;
				}
		   } 
		}
        /*   Update Lagrange multipliers  (lambdar,lambdag,lambdab).      */
//		GetPPIGradients(nver,ur);
//		GetPPIGradients1(nver,ug);
//		GetPPIGradients2(nver,ub);

		for(i=0;i<ntri;i++)
		{
            lambdar[i].x+= penParam*(pr[i].x-m_triangles[i].grad.x);
			lambdar[i].y+= penParam*(pr[i].y-m_triangles[i].grad.y);
			lambdar[i].z+= penParam*(pr[i].z-m_triangles[i].grad.z);

            lambdag[i].x+= penParam*(pg[i].x-m_triangles[i].grad1.x);
			lambdag[i].y+= penParam*(pg[i].y-m_triangles[i].grad1.y);
			lambdag[i].z+= penParam*(pg[i].z-m_triangles[i].grad1.z);

            lambdab[i].x+= penParam*(pb[i].x-m_triangles[i].grad2.x);
			lambdab[i].y+= penParam*(pb[i].y-m_triangles[i].grad2.y);
			lambdab[i].z+= penParam*(pb[i].z-m_triangles[i].grad2.z);
		}
		/*   Compute the stopping condition     */
        stoppingCond = 0;
/*		VECTOR3D tempr,tempg,tempb;
		for(i=0;i<ntri;i++)
		{
			tempr.x = pr[i].x-m_triangles[i].grad.x;
			tempr.y = pr[i].y-m_triangles[i].grad.y;
			tempr.z = pr[i].z-m_triangles[i].grad.z;
			
			tempg.x = pg[i].x-m_triangles[i].grad1.x;
			tempg.y = pg[i].y-m_triangles[i].grad1.y;
			tempg.z = pg[i].z-m_triangles[i].grad1.z;
			
			tempb.x = pb[i].x-m_triangles[i].grad2.x;
			tempb.y = pb[i].y-m_triangles[i].grad2.y;
			tempb.z = pb[i].z-m_triangles[i].grad2.z;

			stoppingCond+= m_trianglesArea[i]*(POWER(tempr.x)+POWER(tempr.y)+POWER(tempr.z)
				+POWER(tempg.x)+POWER(tempg.y)+POWER(tempg.z)
				+POWER(tempb.x)+POWER(tempb.y)+POWER(tempb.z));
		}
*/
		for(i=0;i<nver;i++)
			stoppingCond+= (POWER(ur[i]-urold[i])+POWER(ug[i]-ugold[i])+POWER(ub[i]-ubold[i])) * m_vertices[i].BCDArea;

	}
	while(stoppingCond>outTole);

	for(i=0;i<nver;i++)
	{
		m_vertices[i].r = ur[i]; m_vertices[i].g = ug[i]; m_vertices[i].b = ub[i];
	}


	delete ur; delete ug; delete ub;
	delete urold; delete ugold; delete ubold;
	delete b;
	delete pr; delete pg; delete pb;
	delete lambdar; delete lambdag; delete lambdab;
}

void TriangularMesh::ALM_TVL2Denoising_usub_nr_linbcg(unsigned long n,double b[],double x[],double fidParam, double penParam,int itol,double tol,int itmax,int* iter,double* err)
{
/* For PBCG.

   Solves the u-sub problem of ALM_TVL2Denoising by PBCG.
   fidParam: the fidelity parameter of the TVL2 model.
   penParam: the penalty parameter of the ALM.
   These two parameters will be involved into the coefficient matrix A.
   Specifically, A = fidParam*Identity - penParam*Laplacian.

   Solves Ax=b for x[1..n],given b[1..n],by the iterative preconditional biconjugate gradient method.
   on input x[1..n] should be set to an initial guess of the solution or all zeros;
   itol is 1,2,3,4, specialfying which convergence test is applied;
   itmax is allowed maximum number of iterations; tol is the desired convergence tolerance;
   on output x[1..n] is set to the improved solution; iter is the number of iterations actually taken;
   err is the estimated error.
   the matrix A is referenced only through the user-supplied routines, which computes
   the product of A or its transpose on a vector; or, which solves for precondition, 
   whose coefficient matrix is usually the diagonal part of A.
 */
    unsigned long j;
    double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zminrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p = new double[n];
	pp = new double[n];
	r = new double[n];
	rr = new double[n];
	z = new double[n];
	zz = new double[n];

	//calculate initial residual.
	*iter = 0;
	ALM_TVL2Denoising_usub_nr_atimes(n,x,r,fidParam,penParam,0);
	for(j=0;j<n;j++)
	{
		r[j] = b[j]-r[j];
		rr[j] = r[j];
	}
	/* ALM_TVL2Denoising_usub_nr_atimes(n,r,rr,fidParam,penParam,0) */
	if(itol==1)
	{
		bnrm = nr_snrm(n,b,itol);
		ALM_TVL2Denoising_usub_nr_asolve(n,r,z,fidParam,penParam,0);
	}
	else if(itol==2)
	{
		ALM_TVL2Denoising_usub_nr_asolve(n,b,z,fidParam,penParam,0);
		bnrm = nr_snrm(n,z,itol);
		ALM_TVL2Denoising_usub_nr_asolve(n,r,z,fidParam,penParam,0);
	}
	else if(itol==3 || itol==4)
	{
		ALM_TVL2Denoising_usub_nr_asolve(n,b,z,fidParam,penParam,0);
		bnrm = nr_snrm(n,z,itol);
		ALM_TVL2Denoising_usub_nr_asolve(n,r,z,fidParam,penParam,0);
		znrm = nr_snrm(n,z,itol);
	}
	else return;//
    while(*iter<=itmax)
	{// main loop.
		++(*iter);
		ALM_TVL2Denoising_usub_nr_asolve(n,rr,zz,fidParam,penParam,1);
		for(bknum=0.0,j=0;j<n;j++)
			bknum+= z[j]*rr[j];
		/*  calculate coefficient bk and direction vectors p and pp  */
		if(*iter==1)
		{
			for(j=0;j<n;j++)
			{
				p[j] = z[j];
				pp[j] = zz[j];
			}
		}
		else
		{
			bk = bknum/bkden;
			for(j=0;j<n;j++)
			{
				p[j] = bk*p[j]+z[j];
				pp[j] = bk*pp[j]+zz[j];
			}
		}
		bkden = bknum;//calculate new coefficients ak, new x, .
		ALM_TVL2Denoising_usub_nr_atimes(n,p,z,fidParam,penParam,0);
		for(akden=0.0,j=0;j<n;j++)
			akden+= z[j]*pp[j];
		ak = bknum/akden;
		ALM_TVL2Denoising_usub_nr_atimes(n,pp,zz,fidParam,penParam,1);
		for(j=0;j<n;j++)
		{
			x[j]+= ak*p[j];
			r[j]-= ak*z[j];
			rr[j]-= ak*zz[j];
		}
		ALM_TVL2Denoising_usub_nr_asolve(n,r,z,fidParam,penParam,0);//
		if(itol==1)
			*err = nr_snrm(n,r,itol)/bnrm;
		else if(itol==2)
			*err = nr_snrm(n,z,itol)/bnrm;
		else if(itol==3 || itol==4)
		{
			zminrm = znrm;
			znrm = nr_snrm(n,z,itol);
			if(fabs(zminrm-znrm)>EPS*znrm)
			{
				dxnrm = fabs(ak)*nr_snrm(n,p,itol);
				*err = znrm/fabs(zminrm-znrm)*dxnrm;
			}
			else
			{
				*err = znrm/bnrm;
				continue;
			}
			xnrm = nr_snrm(n,x,itol);
			if(*err<=0.5*xnrm) *err/= xnrm;
			else
			{
				*err = znrm/bnrm;
				continue;
			}
		}
		if(*err<=tol)
			break;
	}
	delete r;
	delete rr;
	delete p;
	delete pp;
	delete z;
	delete zz;
}

void TriangularMesh::ALM_TVgL1_usub_nr_linbcg(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itol,double tol,int itmax,int* iter,double* err)
{
/* For PBCG.

   Solves the u-sub problem of ALM_TVgL1 by PBCG.
   penParam_p: the penalty parameter for p of the ALM.
   penParam_z: the penalty parameter for z of the ALM.
   These two parameters will be involved into the coefficient matrix A.
   Specifically, A = penParam_z*Identity - penParam_p*Laplacian.

   Solves Ax=b for x[1..n],given b[1..n],by the iterative preconditional biconjugate gradient method.
   on input x[1..n] should be set to an initial guess of the solution or all zeros;
   itol is 1,2,3,4, specialfying which convergence test is applied;
   itmax is allowed maximum number of iterations; tol is the desired convergence tolerance;
   on output x[1..n] is set to the improved solution; iter is the number of iterations actually taken;
   err is the estimated error.
   the matrix A is referenced only through the user-supplied routines, which computes
   the product of A or its transpose on a vector; or, which solves for precondition, 
   whose coefficient matrix is usually the diagonal part of A.
 */
    unsigned long j;
    double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zminrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p = new double[n];
	pp = new double[n];
	r = new double[n];
	rr = new double[n];
	z = new double[n];
	zz = new double[n];

	//calculate initial residual.
	*iter = 0;
	ALM_TVgL1_usub_nr_atimes(n,x,r,penParam_p,penParam_z,0);
	for(j=0;j<n;j++)
	{
		r[j] = b[j]-r[j];
		rr[j] = r[j];
	}
 
	if(itol==1)
	{
		bnrm = nr_snrm(n,b,itol);
		ALM_TVgL1_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);
	}
	else if(itol==2)
	{
		ALM_TVgL1_usub_nr_asolve(n,b,z,penParam_p,penParam_z,0);
		bnrm = nr_snrm(n,z,itol);
		ALM_TVgL1_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);
	}
	else if(itol==3 || itol==4)
	{
		ALM_TVgL1_usub_nr_asolve(n,b,z,penParam_p,penParam_z,0);
		bnrm = nr_snrm(n,z,itol);
		ALM_TVgL1_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);
		znrm = nr_snrm(n,z,itol);
	}
	else return;//
    while(*iter<=itmax)
	{// main loop.
		++(*iter);
		ALM_TVgL1_usub_nr_asolve(n,rr,zz,penParam_p,penParam_z,1);
		for(bknum=0.0,j=0;j<n;j++)
			bknum+= z[j]*rr[j];
		/*  calculate coefficient bk and direction vectors p and pp  */
		if(*iter==1)
		{
			for(j=0;j<n;j++)
			{
				p[j] = z[j];
				pp[j] = zz[j];
			}
		}
		else
		{
			bk = bknum/bkden;
			for(j=0;j<n;j++)
			{
				p[j] = bk*p[j]+z[j];
				pp[j] = bk*pp[j]+zz[j];
			}
		}
		bkden = bknum;//calculate new coefficients ak, new x, .
		ALM_TVgL1_usub_nr_atimes(n,p,z,penParam_p,penParam_z,0);
		for(akden=0.0,j=0;j<n;j++)
			akden+= z[j]*pp[j];
		ak = bknum/akden;
		ALM_TVgL1_usub_nr_atimes(n,pp,zz,penParam_p,penParam_z,1);
		for(j=0;j<n;j++)
		{
			x[j]+= ak*p[j];
			r[j]-= ak*z[j];
			rr[j]-= ak*zz[j];
		}
		ALM_TVgL1_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);//
		if(itol==1)
			*err = nr_snrm(n,r,itol)/bnrm;
		else if(itol==2)
			*err = nr_snrm(n,z,itol)/bnrm;
		else if(itol==3 || itol==4)
		{
			zminrm = znrm;
			znrm = nr_snrm(n,z,itol);
			if(fabs(zminrm-znrm)>EPS*znrm)
			{
				dxnrm = fabs(ak)*nr_snrm(n,p,itol);
				*err = znrm/fabs(zminrm-znrm)*dxnrm;
			}
			else
			{
				*err = znrm/bnrm;
				continue;
			}
			xnrm = nr_snrm(n,x,itol);
			if(*err<=0.5*xnrm) *err/= xnrm;
			else
			{
				*err = znrm/bnrm;
				continue;
			}
		}
		if(*err<=tol)
			break;
	}
	delete r;
	delete rr;
	delete p;
	delete pp;
	delete z;
	delete zz;
}

void TriangularMesh::ALM_MulRegionLabelling_usub_nr_linbcg(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itol,double tol,int itmax,int* iter,double* err)
{
/* For PBCG.

   Solves the u-sub problem of ALM_MulRegionLabelling by PBCG.
   penParam_p: the penalty parameter for p of the ALM.
   penParam_z: the penalty parameter for z of the ALM.
   These two parameters will be involved into the coefficient matrix A.
   Specifically, A = penParam_z*Identity - penParam_p*Laplacian.

   Solves Ax=b for x[1..n],given b[1..n],by the iterative preconditional biconjugate gradient method.
   on input x[1..n] should be set to an initial guess of the solution or all zeros;
   itol is 1,2,3,4, specialfying which convergence test is applied;
   itmax is allowed maximum number of iterations; tol is the desired convergence tolerance;
   on output x[1..n] is set to the improved solution; iter is the number of iterations actually taken;
   err is the estimated error.
   the matrix A is referenced only through the user-supplied routines, which computes
   the product of A or its transpose on a vector; or, which solves for precondition, 
   whose coefficient matrix is usually the diagonal part of A.
 */
    unsigned long j;
    double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zminrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p = new double[n];
	pp = new double[n];
	r = new double[n];
	rr = new double[n];
	z = new double[n];
	zz = new double[n];

	//calculate initial residual.
	*iter = 0;
	ALM_MulRegionLabelling_usub_nr_atimes(n,x,r,penParam_p,penParam_z,0);
	for(j=0;j<n;j++)
	{
		r[j] = b[j]-r[j];
		rr[j] = r[j];
	}
 
	if(itol==1)
	{
		bnrm = nr_snrm(n,b,itol);
		ALM_MulRegionLabelling_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);
	}
	else if(itol==2)
	{
		ALM_MulRegionLabelling_usub_nr_asolve(n,b,z,penParam_p,penParam_z,0);
		bnrm = nr_snrm(n,z,itol);
		ALM_MulRegionLabelling_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);
	}
	else if(itol==3 || itol==4)
	{
		ALM_MulRegionLabelling_usub_nr_asolve(n,b,z,penParam_p,penParam_z,0);
		bnrm = nr_snrm(n,z,itol);
		ALM_MulRegionLabelling_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);
		znrm = nr_snrm(n,z,itol);
	}
	else return;//
    while(*iter<=itmax)
	{// main loop.
		++(*iter);
		ALM_MulRegionLabelling_usub_nr_asolve(n,rr,zz,penParam_p,penParam_z,1);
		for(bknum=0.0,j=0;j<n;j++)
			bknum+= z[j]*rr[j];
		/*  calculate coefficient bk and direction vectors p and pp  */
		if(*iter==1)
		{
			for(j=0;j<n;j++)
			{
				p[j] = z[j];
				pp[j] = zz[j];
			}
		}
		else
		{
			bk = bknum/bkden;
			for(j=0;j<n;j++)
			{
				p[j] = bk*p[j]+z[j];
				pp[j] = bk*pp[j]+zz[j];
			}
		}
		bkden = bknum;//calculate new coefficients ak, new x, .
		ALM_MulRegionLabelling_usub_nr_atimes(n,p,z,penParam_p,penParam_z,0);
		for(akden=0.0,j=0;j<n;j++)
			akden+= z[j]*pp[j];
		ak = bknum/akden;
		ALM_MulRegionLabelling_usub_nr_atimes(n,pp,zz,penParam_p,penParam_z,1);
		for(j=0;j<n;j++)
		{
			x[j]+= ak*p[j];
			r[j]-= ak*z[j];
			rr[j]-= ak*zz[j];
		}
		ALM_MulRegionLabelling_usub_nr_asolve(n,r,z,penParam_p,penParam_z,0);//
		if(itol==1)
			*err = nr_snrm(n,r,itol)/bnrm;
		else if(itol==2)
			*err = nr_snrm(n,z,itol)/bnrm;
		else if(itol==3 || itol==4)
		{
			zminrm = znrm;
			znrm = nr_snrm(n,z,itol);
			if(fabs(zminrm-znrm)>EPS*znrm)
			{
				dxnrm = fabs(ak)*nr_snrm(n,p,itol);
				*err = znrm/fabs(zminrm-znrm)*dxnrm;
			}
			else
			{
				*err = znrm/bnrm;
				continue;
			}
			xnrm = nr_snrm(n,x,itol);
			if(*err<=0.5*xnrm) *err/= xnrm;
			else
			{
				*err = znrm/bnrm;
				continue;
			}
		}
		if(*err<=tol)
			break;
	}
	delete r;
	delete rr;
	delete p;
	delete pp;
	delete z;
	delete zz;
}

void TriangularMesh::ALM_TVL2Denoising_usub_nr_atimes(unsigned long n,double x[],double r[],double fidParam, double penParam,int itrnsp)
{
/* For PBCG.
 */
    if(itrnsp)
		ALM_TVL2Denoising_usub_nr_dsprsax(n,x,r,fidParam,penParam);//转置乘。
	else
	    ALM_TVL2Denoising_usub_nr_dsprsax(n,x,r,fidParam,penParam);
}

void TriangularMesh::ALM_TVgL1_usub_nr_atimes(unsigned long n,double x[],double r[],double penParam_p, double penParam_z,int itrnsp)
{
/* For PBCG.
 */
    if(itrnsp)
		ALM_TVgL1_usub_nr_dsprsax(n,x,r,penParam_p,penParam_z);//转置乘。
	else
	    ALM_TVgL1_usub_nr_dsprsax(n,x,r,penParam_p,penParam_z);
}

void TriangularMesh::ALM_MulRegionLabelling_usub_nr_atimes(unsigned long n,double x[],double r[],double penParam_p, double penParam_z,int itrnsp)
{
/* For PBCG.
 */
    if(itrnsp)
		ALM_MulRegionLabelling_usub_nr_dsprsax(n,x,r,penParam_p,penParam_z);//转置乘。
	else
	    ALM_MulRegionLabelling_usub_nr_dsprsax(n,x,r,penParam_p,penParam_z);
}

void TriangularMesh::ALM_TVL2Denoising_usub_nr_dsprsax(unsigned long n,double x[],double b[],double fidParam, double penParam)
{
/* For PBCG.
 */
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<n;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]-= penParam*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= (fidParam*m_vertices[i].BCDArea-penParam*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::ALM_TVgL1_usub_nr_dsprsax(unsigned long n,double x[],double b[],double penParam_p, double penParam_z)
{
/* For PBCG.
 */
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<n;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]-= penParam_p*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= (penParam_z*m_vertices[i].BCDArea-penParam_p*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::ALM_MulRegionLabelling_usub_nr_dsprsax(unsigned long n,double x[],double b[],double penParam_p, double penParam_z)
{
/* For PBCG.
 */
	unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<n;i++)
	{
		b[i] = 0;
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{//计算非对角元乘积之和。
			b[i]-= penParam_p*or[j].co*x[or[j].iver];
			diag-= or[j].co;
		}
		b[i]+= (penParam_z*m_vertices[i].BCDArea-penParam_p*diag)*x[i]; //对角项。
	}
}

void TriangularMesh::ALM_TVL2Denoising_usub_nr_asolve(unsigned long n,double b[],double x[],double fidParam, double penParam,int itrnsp)
{
/* For PBCG.
 */
    unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<n;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = fidParam*m_vertices[i].BCDArea-penParam*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::ALM_TVgL1_usub_nr_asolve(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itrnsp)
{
/* For PBCG.
 */
    unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<n;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = penParam_z*m_vertices[i].BCDArea-penParam_p*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::ALM_MulRegionLabelling_usub_nr_asolve(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itrnsp)
{
/* For PBCG.
 */
    unsigned long i,j;
	ONERING or;
	double diag;
	for(i=0;i<n;i++)
	{
		diag = 0;
		or = m_ONERINGS[i];
		for(j=0;j<or.size();j++)
		{
			diag-= or[j].co;
		}
		diag = penParam_z*m_vertices[i].BCDArea-penParam_p*diag;
		x[i] = (diag!=0.0 ? b[i]/diag : b[i]);
	}
}

void TriangularMesh::GetPPIGradients(double q[],VECTOR3D *g)
{
// be sure the index of q is the same as the index of m_vertices.
// and the index of g is the same as the index of m_triangles.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		g[i].x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		g[i].y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		g[i].z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

void TriangularMesh::GetDivergence(long ntri,VECTOR3D vf[])
{
// Given a vector field vf on the mesh, to compute its divergence, which reaches a value at each vertex of the mesh.
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double div;
	for(i=0;i<m_vnum;i++)
	{
		div = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				div-= m_trianglesArea[dt.itri]*DotProduct(vf[dt.itri],m_TPPIBG[dt.itri].v0);
			else if(i==m_triangles[dt.itri].ver1)
				div-= m_trianglesArea[dt.itri]*DotProduct(vf[dt.itri],m_TPPIBG[dt.itri].v1);
			else
				div-= m_trianglesArea[dt.itri]*DotProduct(vf[dt.itri],m_TPPIBG[dt.itri].v2);
		}
		div/= m_vertices[i].BCDArea;
		m_vertices[i].divergence = div;
	}
}

void TriangularMesh::ALM_vTVL2Denoising_usub_nr_linbcg(unsigned long n,double b[],double x[],double fidParam, double penParam,int itol,double tol,int itmax,int* iter,double* err)
{
	ALM_TVL2Denoising_usub_nr_linbcg(n,b,x,fidParam,penParam,itol,tol,itmax,iter,err);
}

void TriangularMesh::ALM_vTVL2Denoising_usub_nr_atimes(unsigned long n,double x[],double r[],double fidParam, double penParam,int itrnsp)
{
	ALM_TVL2Denoising_usub_nr_atimes(n,x,r,fidParam,penParam,itrnsp);
}

void TriangularMesh::ALM_vTVL2Denoising_usub_nr_dsprsax(unsigned long n,double x[],double b[],double fidParam, double penParam)
{
	ALM_TVL2Denoising_usub_nr_dsprsax(n,x,b,fidParam,penParam);
}

void TriangularMesh::ALM_vTVL2Denoising_usub_nr_asolve(unsigned long n,double b[],double x[],double fidParam, double penParam,int itrnsp)
{
	ALM_TVL2Denoising_usub_nr_asolve(n,b,x,fidParam,penParam,itrnsp);
}

double TriangularMesh::GaussianRandom()
{
	double x1, x2, w, y1, y2;

    do
	{
        x1 = 2.0 * double(rand())/RAND_MAX - 1.0;
        x2 = 2.0 * double(rand())/RAND_MAX - 1.0;
        w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    
	return y1;
}

void TriangularMesh::AddGaussianNoise(double variance)
{
	unsigned long i;
	double rd;
	for(i=0;i<m_vnum;i++)
	{
		rd = variance*GaussianRandom();
		m_vertices[i].r+= rd;
		rd = variance*GaussianRandom();
		m_vertices[i].g+= rd;
		rd = variance*GaussianRandom();
		m_vertices[i].b+= rd;
	}
}

void TriangularMesh::AddGaussianNoise2VertexGeometry(double variance)
{
	unsigned long i;
	double rd;
	for(i=0;i<m_vnum;i++)
	{
		rd = variance*GaussianRandom();
		m_vertices[i].x+= rd;
		rd = variance*GaussianRandom();
		m_vertices[i].y+= rd;
		rd = variance*GaussianRandom();
		m_vertices[i].z+= rd;
	}
}

void TriangularMesh::ALM_vTVL2NormalFiltering()
{
/*      
	Normal filtering: Solve vectorial TV L2 denoising problem using ALM, and then normalization to get unit length normal vector.

    The problem:
	    min_u \sum|\nabla u|s_tau + fidParam/2 \|u-f\|_2.

    Here u=(unx,uny,unz) is the normal vector of the mesh.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   outTole  : the tolerance for outer iteration;
 */
	
    unsigned long nver,ntri,i;
	double *unx,*uny,*unz,*b;
	VECTOR3D *pnx,*pny,*pnz,*lambdanx,*lambdany,*lambdanz,wnx,wny,wnz,w;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 300;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	unx = new double[nver];
	uny = new double[nver];
	unz = new double[nver];
	b = new double[nver];
	pnx = new VECTOR3D[ntri];
	pny = new VECTOR3D[ntri];
	pnz = new VECTOR3D[ntri];
	lambdanx = new VECTOR3D[ntri];
	lambdany = new VECTOR3D[ntri];
	lambdanz = new VECTOR3D[ntri];

    unsigned int innerL,l;
	double fidParam, penParam, outTole, stoppingCond;
/*   Parameters setting  */
	innerL = 1;
	fidParam = m_Lambda*3;
	penParam = 0.5;
	outTole = 10;

/*   Initialization   */
	for(i=0;i<nver;i++)
	{	
		unx[i] = 0; uny[i] = 0; unz[i] = 0;
	}
	for(i=0;i<ntri;i++)
	{
		pnx[i].x=0; lambdanx[i].x=0;
		pnx[i].y=0; lambdanx[i].y=0;
		pnx[i].z=0; lambdanx[i].z=0;

		pny[i].x=0; lambdany[i].x=0;
		pny[i].y=0; lambdany[i].y=0;
		pny[i].z=0; lambdany[i].z=0;

		pnz[i].x=0; lambdanz[i].x=0;
		pnz[i].y=0; lambdanz[i].y=0;
		pnz[i].z=0; lambdanz[i].z=0;
	}
/*   Iteration        */
	do
	{
		/*  Inner iteration             */
		for(l=0;l<innerL;l++)
		{
		/*  Solve the u-sub problem  (solve for unx,uny,unz).   */

		   for(i=0;i<ntri;i++)
		   {
		    	pnx[i].x = lambdanx[i].x+penParam*pnx[i].x;
			    pnx[i].y = lambdanx[i].y+penParam*pnx[i].y;
			    pnx[i].z = lambdanx[i].z+penParam*pnx[i].z;
		   }
		   GetDivergence(ntri,pnx);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_vertices[i].normal_x-m_vertices[i].divergence);
           ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,unx,fidParam,penParam,itol,tol,itmax,&iter,&err);

		   for(i=0;i<ntri;i++)
		   {
		    	pny[i].x = lambdany[i].x+penParam*pny[i].x;
			    pny[i].y = lambdany[i].y+penParam*pny[i].y;
			    pny[i].z = lambdany[i].z+penParam*pny[i].z;
		   }
		   GetDivergence(ntri,pny);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_vertices[i].normal_y-m_vertices[i].divergence);
           ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,uny,fidParam,penParam,itol,tol,itmax,&iter,&err);

		   for(i=0;i<ntri;i++)
		   {
		    	pnz[i].x = lambdanz[i].x+penParam*pnz[i].x;
			    pnz[i].y = lambdanz[i].y+penParam*pnz[i].y;
			    pnz[i].z = lambdanz[i].z+penParam*pnz[i].z;
		   }
		   GetDivergence(ntri,pnz);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_vertices[i].normal_z-m_vertices[i].divergence);
           ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,unz,fidParam,penParam,itol,tol,itmax,&iter,&err);

        /*  Solve the p-sub problem  (solve for pnx,pny,pnz).   */
           
		   GetPPIGradients(nver,unx);
		   GetPPIGradients1(nver,uny);
		   GetPPIGradients2(nver,unz);

		   for(i=0;i<ntri;i++)
		   {
			    wnx.x = m_triangles[i].grad.x - lambdanx[i].x/penParam;
		    	wnx.y = m_triangles[i].grad.y - lambdanx[i].y/penParam;
			    wnx.z = m_triangles[i].grad.z - lambdanx[i].z/penParam;

			    wny.x = m_triangles[i].grad1.x - lambdany[i].x/penParam;
		    	wny.y = m_triangles[i].grad1.y - lambdany[i].y/penParam;
			    wny.z = m_triangles[i].grad1.z - lambdany[i].z/penParam;

			    wnz.x = m_triangles[i].grad2.x - lambdanz[i].x/penParam;
		    	wnz.y = m_triangles[i].grad2.y - lambdanz[i].y/penParam;
			    wnz.z = m_triangles[i].grad2.z - lambdanz[i].z/penParam;

				if(DotProduct(wnx,wnx)+DotProduct(wny,wny)+DotProduct(wnz,wnz) <= POWER(1/penParam))
				{
				     pnx[i].x = 0; pny[i].x = 0; pnz[i].x = 0;
				     pnx[i].y = 0; pny[i].y = 0; pnz[i].y = 0;
				     pnx[i].z = 0; pny[i].z = 0; pnz[i].z = 0;
				}
			    else
				{
					 double tempNorm;
					 tempNorm = sqrt(DotProduct(wnx,wnx)+DotProduct(wny,wny)+DotProduct(wnz,wnz));

					 pnx[i].x = (1-1/penParam/tempNorm)*wnx.x; pny[i].x = (1-1/penParam/tempNorm)*wny.x; pnz[i].x = (1-1/penParam/tempNorm)*wnz.x;
				     pnx[i].y = (1-1/penParam/tempNorm)*wnx.y; pny[i].y = (1-1/penParam/tempNorm)*wny.y; pnz[i].y = (1-1/penParam/tempNorm)*wnz.y;
				     pnx[i].z = (1-1/penParam/tempNorm)*wnx.z; pny[i].z = (1-1/penParam/tempNorm)*wny.z; pnz[i].z = (1-1/penParam/tempNorm)*wnz.z;
				}
		   } 
		}
        /*   Update Lagrange multipliers  (lambdanx,lambdany,lambdanz).      */

		for(i=0;i<ntri;i++)
		{
            lambdanx[i].x+= penParam*(pnx[i].x-m_triangles[i].grad.x);
			lambdanx[i].y+= penParam*(pnx[i].y-m_triangles[i].grad.y);
			lambdanx[i].z+= penParam*(pnx[i].z-m_triangles[i].grad.z);

            lambdany[i].x+= penParam*(pny[i].x-m_triangles[i].grad1.x);
			lambdany[i].y+= penParam*(pny[i].y-m_triangles[i].grad1.y);
			lambdany[i].z+= penParam*(pny[i].z-m_triangles[i].grad1.z);

            lambdanz[i].x+= penParam*(pnz[i].x-m_triangles[i].grad2.x);
			lambdanz[i].y+= penParam*(pnz[i].y-m_triangles[i].grad2.y);
			lambdanz[i].z+= penParam*(pnz[i].z-m_triangles[i].grad2.z);
		}
		/*   Compute the stopping condition     */
        stoppingCond = 0;
		VECTOR3D tempnx,tempny,tempnz;
		for(i=0;i<ntri;i++)
		{
			tempnx.x = pnx[i].x-m_triangles[i].grad.x;
			tempnx.y = pnx[i].y-m_triangles[i].grad.y;
			tempnx.z = pnx[i].z-m_triangles[i].grad.z;
			
			tempny.x = pny[i].x-m_triangles[i].grad1.x;
			tempny.y = pny[i].y-m_triangles[i].grad1.y;
			tempny.z = pny[i].z-m_triangles[i].grad1.z;
			
			tempnz.x = pnz[i].x-m_triangles[i].grad2.x;
			tempnz.y = pnz[i].y-m_triangles[i].grad2.y;
			tempnz.z = pnz[i].z-m_triangles[i].grad2.z;

			stoppingCond+= m_trianglesArea[i]*(POWER(tempnx.x)+POWER(tempnx.y)+POWER(tempnx.z)
				+POWER(tempny.x)+POWER(tempny.y)+POWER(tempny.z)
				+POWER(tempnz.x)+POWER(tempnz.y)+POWER(tempnz.z));
		}
	}
	while(stoppingCond>outTole);

    VECTOR3D normal;
	for(i=0;i<nver;i++)
	{
		normal.x = unx[i]; normal.y = uny[i]; normal.z = unz[i];
		NormalizeVector3D(&normal);
		m_vertices[i].normal_x = normal.x; 
		m_vertices[i].normal_y = normal.y; 
		m_vertices[i].normal_z = normal.z;
	}


	delete unx; delete uny; delete unz;
	delete b;
	delete pnx; delete pny; delete pnz;
	delete lambdanx; delete lambdany; delete lambdanz;
}

void TriangularMesh::AddGaussianNoise2VertexNormal(double variance)
{
	unsigned long i;
	double rd;
	for(i=0;i<m_vnum;i++)
	{
		rd = variance*GaussianRandom();
		m_vertices[i].normal_x+= rd;
		rd = variance*GaussianRandom();
		m_vertices[i].normal_y+= rd;
		rd = variance*GaussianRandom();
		m_vertices[i].normal_z+= rd;
	}
}

bool TriangularMesh::HasObtuseAngle()
{
/*  This function is to test whether the loaded triangle mesh has obtuse triangles. */
	int i;
	VECTOR3D v1,v2;
	for(i=0;i<m_trinum;i++)
	{
		v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver0]),GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver1]));
		v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver0]),GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver2]));
		if (DotProduct(v1,v2)<0)
			return true;
        v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver1]),GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver0]));
		v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver1]),GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver2]));
		if (DotProduct(v1,v2)<0)
			return true;
		v1 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver2]),GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver1]));
		v2 = GetVector3Dfrom2Vertices(GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver2]),GetVertex3DFromPoint3D(m_vertices[m_triangles[i].ver0]));
		if (DotProduct(v1,v2)<0)
			return true;
	}
	return false;
}

void TriangularMesh::AddSaltPepperNoise(float level)
{
	int num;
	unsigned int i,rd;
	num = floor(m_vnum*level);

	for(i=0;i<num;i++)
	{
		rd = floor(double(rand())/RAND_MAX*m_vnum);
		if(i%2==0)
		{
		    m_vertices[rd].r = 1;
		    m_vertices[rd].g = 1;
		    m_vertices[rd].b = 1;
		}
		else
		{
			m_vertices[rd].r = 0;
		    m_vertices[rd].g = 0;
		    m_vertices[rd].b = 0;
		}
	}
}

void TriangularMesh::ALM_MulRegionLabelling_sDev_grayimage(double *u, double *s, unsigned long nRegion, unsigned long nver)
{
	int iver,iRegion;
	double mean,denom;

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
		mean = 0;
		denom = 0;
		for(iver=0;iver<nver;iver++)
		{
			mean+= u[iRegion*nver+iver]*m_vertices[iver].r*m_vertices[iver].BCDArea;
			denom+= u[iRegion*nver+iver]*m_vertices[iver].BCDArea;
		}
		mean/= denom;
		for(iver=0;iver<nver;iver++)
		{
			*(s+iRegion*nver+iver) = POWER(m_vertices[iver].r-mean);
		}
	}
}

void TriangularMesh::ALM_MulRegionLabelling_sDev_grayimage_givenMean(double *mean, double *s, unsigned long nRegion, unsigned long nver)
{
	int iver,iRegion;

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
		for(iver=0;iver<nver;iver++)
		{
			*(s+iRegion*nver+iver) = POWER(m_vertices[iver].r-mean[iRegion]);
		}
	}
}

void TriangularMesh::ALM_MulRegionLabelling_sDev_colorimage(double *u, double *s, unsigned long nRegion, unsigned long nver)
{
	int iver,iRegion;
	double meanr,meang,meanb,denom;

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
		meanr = 0;
		meang = 0;
		meanb = 0;
		denom = 0;
		for(iver=0;iver<nver;iver++)
		{
			meanr+= u[iRegion*nver+iver]*m_vertices[iver].r*m_vertices[iver].BCDArea;
			meang+= u[iRegion*nver+iver]*m_vertices[iver].g*m_vertices[iver].BCDArea;
			meanb+= u[iRegion*nver+iver]*m_vertices[iver].b*m_vertices[iver].BCDArea;
			denom+= u[iRegion*nver+iver]*m_vertices[iver].BCDArea;
		}
		meanr/= denom;
		meang/= denom;
		meanb/= denom;
		for(iver=0;iver<nver;iver++)
		{
			*(s+iRegion*nver+iver) = POWER(m_vertices[iver].r-meanr)+POWER(m_vertices[iver].g-meang)
				+POWER(m_vertices[iver].b-meanb);
		}
	}
}

void TriangularMesh::ALM_MulRegionLabelling_sDev_colorimage_givenMean(double *meanr, double *meang, double *meanb, double *s, unsigned long nRegion, unsigned long nver)
{
	int iver,iRegion;

	for(iRegion=0;iRegion<nRegion;iRegion++)
	{
		for(iver=0;iver<nver;iver++)
		{
			*(s+iRegion*nver+iver) = pow(m_vertices[iver].r-meanr[iRegion],2)+pow(m_vertices[iver].g-meang[iRegion],2)
				+pow(m_vertices[iver].b-meanb[iRegion],2);
		}
	}
}

void TriangularMesh::ALM_MulRegionLabelling_Project2K(double *u, unsigned int nRegion, unsigned long nver)
{
// to be implemented.
	int iRegion,iver;
	bool * ISet;
	bool bContinue;
	double cardinality;
	ISet = new bool[nRegion];

	for(iver=0;iver<nver;iver++)
	{
        memset(ISet, false, sizeof(bool)*nRegion);
        bContinue = true;
		cardinality = nRegion;
        while(bContinue)
        {
            bContinue = false;
            double sum = 0;
            for(iRegion = 0; iRegion < nRegion; iRegion++)
            {
                sum += u[iRegion*nver+iver];
            }
            for(iRegion = 0; iRegion < nRegion; iRegion++)
            {
                if(!(ISet[iRegion]))
                {
                    u[iRegion*nver+iver] -= (sum - 1.0)/cardinality;
                    if(u[iRegion*nver+iver] < 0)
                    {
                        bContinue = true;
                        ISet[iRegion] = true;
                        u[iRegion*nver+iver] = 0;
                        cardinality -= 1.0;
                    }
                }
                else
                {
                    u[iRegion*nver+iver] = 0;
                }
            }
        }
	}
}

void TriangularMesh::ALM_MulRegionLabelling_Binarization(double *u, unsigned int nRegion, unsigned long nver)
{
	int iRegion;
	long iver;
	double normTemp;
	for(iver=0;iver<nver;iver++)
	{
		normTemp = 0;
        for(iRegion = 0; iRegion < nRegion; iRegion++)
        {
            normTemp += POWER(u[iRegion*nver+iver]);
        }

        int id;
        double MValueMin = 1.0e8;
        for(iRegion = 0; iRegion < nRegion; iRegion++)
        {
            double Ni = normTemp - POWER(u[iRegion*nver+iver]) + POWER(1.0 - u[iRegion*nver+iver]);
            if (Ni < MValueMin)
            {
                MValueMin = Ni;
                id = iRegion;
            }
        }
        m_vertices[iver].RegionId = id;
	}
}


// add by duan qi for alm_mesh_refinement
bool TriangularMesh::LoadMeshFile(const char* filename)
{
	OpenMesh::IO::Options read_options; 
	bool ok = OpenMesh::IO::read_mesh(this->m_ObjTriMesh, string(filename), read_options);
	if (!ok) {
		cout << "Error in load the off model " << filename << endl;
		return false;
	}
	cout << "#Vertex: " << this->m_ObjTriMesh.n_vertices() << ",  #Edges: " << this->m_ObjTriMesh.n_edges() << ",  #Faces: " << this->m_ObjTriMesh.n_faces() << endl;
	POINT3d tmp;
	TRIANGLE tmp2;
	COLORRGB crgb;
	int feature_id,ss,tt;//
	double s0,t0,s1,t1,s2,t2;//variables for loading data models from Turk.
	//ifp.open(filename);
	m_vertices.clear();
	m_triangles.clear();
	m_tempRGB.clear();
	m_tempContour.clear();

	m_WeightedCurveLength.clear();//
	m_tSteps.clear();//
	m_WeightedCurveLength_Forobservation.clear();//
	m_tSteps_Forobservation.clear();//
	m_CPUtime.clear();//

	m_vnum = this->m_ObjTriMesh.n_vertices(); m_trinum = this->m_ObjTriMesh.n_faces();
	for (MyMesh::VertexIter v_it = this->m_ObjTriMesh.vertices_begin();v_it != this->m_ObjTriMesh.vertices_end(); ++v_it) {
		tmp.x = this->m_ObjTriMesh.point(v_it.handle())[0];
		tmp.y = this->m_ObjTriMesh.point(v_it.handle())[1];
		tmp.z = this->m_ObjTriMesh.point(v_it.handle())[2];
		tmp.RegionId = 0;
		m_vertices.push_back(tmp);
	}
	for (MyMesh::FaceIter f_it = this->m_ObjTriMesh.faces_begin(); f_it != this->m_ObjTriMesh.faces_end(); ++f_it) {
		MyMesh::ConstFaceVertexIter cfv_it = this->m_ObjTriMesh.cfv_iter(f_it.handle());
		tmp2.ver0 = cfv_it.handle().idx(); ++cfv_it;
		tmp2.ver1 = cfv_it.handle().idx(); ++cfv_it;
		tmp2.ver2 = cfv_it.handle().idx();
		m_triangles.push_back(tmp2);
	}
	//UnifyData();
	//for (MyMesh::VertexIter v_it = this->m_ObjTriMesh.vertices_begin();v_it != this->m_ObjTriMesh.vertices_end(); ++v_it) {
	//	this->m_ObjTriMesh.set_point(v_it, OpenMesh::Vec3f(m_vertices[v_it.handle().idx()].x, m_vertices[v_it.handle().idx()].y, m_vertices[v_it.handle().idx()].z));
	//}
	this->m_ObjTriMesh.request_face_normals();		this->m_ObjTriMesh.update_face_normals();
	this->m_ObjTriMesh.request_vertex_normals();	this->m_ObjTriMesh.update_vertex_normals();
	for(int i=0;i<m_vnum;i++)
	{
		crgb.r = m_vertices[i].x;	m_vertices[i].ref_x = m_vertices[i].x;
		crgb.g = m_vertices[i].y;	m_vertices[i].ref_y = m_vertices[i].y;
		crgb.b = m_vertices[i].z;	m_vertices[i].ref_z = m_vertices[i].z;
		m_tempRGB.push_back(crgb);
	}
	for(int i=0;i<m_vnum;i++)
	{//initialize the normals of vertices.
		m_vertices[i].normal_x = 0;
		m_vertices[i].normal_y = 0;
		m_vertices[i].normal_z = 0;
	}
	FindNormals();

	BuildOneDisks();
	//  BuildCCDualVerticesForTriangles();//this is to construct cc dual.
	BuildBaryCentersForTriangles();//this is to construct bc dual.
	BuildTPPIBG();//to compute the gradient of basis functions of primal-primal interpolation.
	//  WriteTPPIBG();
	BuildONEDISKS();
	BuildONEDISKSANISOTROPIC();
	ComputeBCDualArea();
	BuildOnerings();
	BuildTrianglesArea();
	if (AnisotropicLaplace) {
		CalculateVertexVoronoiArea();
	}

	BuildEdgeIndicator(2);

	m_WeightedCurveLength.push_back(ComputeWeightedCurveLength());
	m_tSteps.push_back(m_tStep);
	m_WeightedCurveLength_Forobservation.push_back(ComputeWeightedCurveLength());
	m_tSteps_Forobservation.push_back(m_tStep);

	bool hasOA = HasObtuseAngle();

	return true;
}

bool TriangularMesh::LoadPSNormalFile(const char* filename)
{
	fstream fin(filename, ios::in);
	if (!fin) {
		cout << "Can not open " << filename<< " to read photometric normal data..." << endl;
		return false;
	}
	int c_size; fin>>c_size;
	if (c_size != m_vnum&&m_vnum>0) {
		cout << "The size of PS normal data "<< c_size <<" is not match with mesh vertex number " << m_vnum << endl;
		return false;
	}
	this->m_ObjTriMesh.update_face_normals();
	this->m_ObjTriMesh.update_vertex_normals();
	OpenMesh::Vec3f v_normal;
	for (int i = 0; i < c_size; i ++) {
		fin >> v_normal[0] >> v_normal[1] >> v_normal[2]; v_normal.normalize();
		m_vertices[i].ps_normal_x = v_normal[0];
		m_vertices[i].ps_normal_y = v_normal[1];
		m_vertices[i].ps_normal_z = v_normal[2];
	}
	fin.close();
	cout << "Load PS normal file: " << filename << " Done..." << endl;
	fstream	fout("CompareNorm.txt", ios::out);
	for (int i = 0; i < m_vnum; ++ i) {
		fout << "ps norm: " << m_vertices[i].ps_normal_x <<" " << m_vertices[i].ps_normal_y <<" " << m_vertices[i].ps_normal_z 
			 <<";  openmesh norm: " << m_ObjTriMesh.normal(MyMesh::VertexHandle(i)).data()[0] <<" " << m_ObjTriMesh.normal(MyMesh::VertexHandle(i)).data()[1] <<" " 
			 << m_ObjTriMesh.normal(MyMesh::VertexHandle(i)).data()[2] << endl;
	}
	fout.close();
	return true;
}

bool TriangularMesh::LoadVertexColorFile(const char* filename, double m_sigma)
{
	fstream fin(filename, ios::in);
	if (!fin) {
		cout << "Can not open " << filename<< " to read color data..." << endl;
		return false;
	}
	int c_size; fin>>c_size;
	if (c_size != m_vnum&&m_vnum>0) {
		cout << "The size of color data " << c_size << " is not match with mesh vertex number " << m_vnum << endl;
		return false;
	}
	// read all visible vertex intensity and corresponding camera center
	this->m_ObjTriMesh.request_vertex_colors();
	MyMesh org_mesh = this->m_ObjTriMesh; 	org_mesh.request_vertex_colors();
	double tintensity, tcx, tcy, tcz; VECTOR3D camra_center;
	for (int i = 0; i < m_vnum; i ++) {
		fin >> c_size; m_vertices[i].intensity_list.clear(); m_vertices[i].Camera_center_list.clear(); m_vertices[i].intensity = -1;
		for (int k = 0; k < c_size; ++ k) {
			fin >> tintensity >> camra_center.x >> camra_center.y >> camra_center.z;
			m_vertices[i].intensity_list.push_back(tintensity);
			m_vertices[i].intensity += tintensity;
			m_vertices[i].Camera_center_list.push_back(camra_center);
		}
		if (c_size > 0) {
			m_vertices[i].intensity = (m_vertices[i].intensity+1)/c_size;
		}
		org_mesh.set_color(MyMesh::VertexHandle(i), MyMesh::Color(m_vertices[i].intensity,m_vertices[i].intensity,m_vertices[i].intensity));
	}
	fin.close();
	cout << "Load color file: " << filename << " Done..." << endl;
	this->UpdateVertexWeight(m_sigma);

	string Obj_prefix = string(filename); 
	Obj_prefix.erase(Obj_prefix.find_last_of("\\"), Obj_prefix.length());
	Obj_prefix = Obj_prefix.substr(Obj_prefix.find_last_of("\\")+1, Obj_prefix.length());

	OpenMesh::IO::Options write_options; 
	write_options.set(OpenMesh::IO::Options::VertexColor); 
	OpenMesh::IO::write_mesh(org_mesh, Obj_prefix+"MeshColorOrg.off", write_options);
	return true;
}

bool TriangularMesh::LoadTargetMeshFile(const char* filename)
{
	OpenMesh::IO::Options read_options; MyMesh ObjTriMesh;
	bool ok = OpenMesh::IO::read_mesh(ObjTriMesh, string(filename), read_options);
	if (!ok) {
		cout << "Error in load the off model " << filename << endl;
		return false;
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		m_vertices[v_it.handle().idx()].tar_x = ObjTriMesh.point(v_it.handle())[0];
		m_vertices[v_it.handle().idx()].tar_y = ObjTriMesh.point(v_it.handle())[1];
		m_vertices[v_it.handle().idx()].tar_z = ObjTriMesh.point(v_it.handle())[2];
	}
	return true;
}

bool TriangularMesh::UpdateVertexWeight(double m_sigma)
{
	//calculate the vertex weight for tv term
	double max_weight = -100;
	for (int i = 0; i < m_vnum; ++ i) {
		int count = 0; int avg_intensity = 0.0;
		for (int idx = 0; idx < m_ONERINGS[i].size(); ++ idx) {
			int temp_nei_id = m_ONERINGS[i][idx].iver;
			count ++; avg_intensity += m_vertices[m_ONERINGS[i][idx].iver].intensity;
		}
		avg_intensity = avg_intensity/count;
		double cur_weight = std::exp(-std::pow((m_vertices[i].intensity - avg_intensity),2.0)/(m_sigma) );
		m_vertices[i].ps_weight = cur_weight;
		if (cur_weight > max_weight) {
			max_weight = cur_weight;
		}
	}
	for (int i = 0; i < m_vnum; ++ i) {
		m_vertices[i].ps_weight = m_vertices[i].ps_weight/max_weight;
	}
	cout << "Update mesh vertex weight for PS normal term Done..." << endl;
#ifdef _DEBUG
	fstream of("VertexWeight.txt",std::ios::out);
	if (!of) {
		cout << "Can not open VertexWeight.txt to save data..." << endl;
		return false;
	}
	for (int i = 0; i < m_vnum; ++ i) {
		m_vertices[i].ps_weight = m_vertices[i].ps_weight/max_weight;
		of << i << ": " << m_vertices[i].ps_weight << endl;
	}
	of.close();
#endif
	return true;
}

void TriangularMesh::GetDivergence(long ntri, vector<VECTOR3D> vf)
{
	// Given a vector field vf on the mesh, to compute its divergence, which reaches a value at each vertex of the mesh.
	DISK d;
	DISKTRIANGLE dt;
	unsigned long i,j;
	double div;
	for(i=0;i<m_vnum;i++)
	{
		div = 0;
		d = m_ONEDISKS[i];//每个顶点对应的ONEDISK.
		for(j=0;j<d.size();j++)
		{//对该顶点的ONEDISK中的所有三角形循环。
			dt = d[j];
			if(i==m_triangles[dt.itri].ver0)
				div-= m_trianglesArea[dt.itri]*DotProduct(vf[dt.itri],m_TPPIBG[dt.itri].v0);
			else if(i==m_triangles[dt.itri].ver1)
				div-= m_trianglesArea[dt.itri]*DotProduct(vf[dt.itri],m_TPPIBG[dt.itri].v1);
			else
				div-= m_trianglesArea[dt.itri]*DotProduct(vf[dt.itri],m_TPPIBG[dt.itri].v2);
		}
		div/= m_vertices[i].BCDArea;
		m_vertices[i].divergence = div;
	}
}

void TriangularMesh::GetPPIGradients(long n, vector<double> q)
{
	//be sure that n is the number of the vertices of the whole mesh.
	//and the index of q is the same as the index of m_vertices.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		m_triangles[i].grad.x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		m_triangles[i].grad.y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		m_triangles[i].grad.z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

void TriangularMesh::GetPPIGradients1(long n, vector<double> q)
{
	//be sure that n is the number of the vertices of the whole mesh.
	//and the index of q is the same as the index of m_vertices.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		m_triangles[i].grad1.x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		m_triangles[i].grad1.y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		m_triangles[i].grad1.z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

void TriangularMesh::GetPPIGradients2(long n, vector<double> q)
{
	//be sure that n is the number of the vertices of the whole mesh.
	//and the index of q is the same as the index of m_vertices.
	long i;
	for(i=0;i<m_trinum;i++)
	{
		m_triangles[i].grad2.x = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.x
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.x+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.x;
		m_triangles[i].grad2.y = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.y
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.y+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.y;
		m_triangles[i].grad2.z = q[m_triangles[i].ver0]*m_TPPIBG[i].v0.z
			+q[m_triangles[i].ver1]*m_TPPIBG[i].v1.z+q[m_triangles[i].ver2]*m_TPPIBG[i].v2.z;
	}
}

double TriangularMesh::CalculateVEnergy(MyMesh &T_Mesh, bool UseFaceArea = false)
{
	if (T_Mesh.n_vertices() != m_vnum) {
		cout << "The size of CurV and InputV is not equal..." << endl;
		return 0.0;
	}
	double result = 0.00; OpenMesh::Vec3f CurV;
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		CurV[0] = m_vertices[v_it.handle().idx()].ref_x; CurV[1] = m_vertices[v_it.handle().idx()].ref_y; CurV[2] = m_vertices[v_it.handle().idx()].ref_z;
		result += (T_Mesh.point(v_it.handle()) - CurV).sqrnorm()*(UseFaceArea?m_vertices[v_it.handle().idx()].BCDArea:1.0);
	}
	return result;
}

double TriangularMesh::CalculateNEnergy(MyMesh &T_Mesh, vector<double> &vec_pervertex_energy, bool UseFaceArea = false)
{
	if (T_Mesh.n_vertices() != m_vnum) {
		cout << "The size of CurN and InputN is not equal..." << endl;
		return 0.0;
	}
	if (vec_pervertex_energy.size() != m_vnum) {
		vec_pervertex_energy.resize(m_vnum);
	}
	T_Mesh.update_face_normals();T_Mesh.update_vertex_normals();
	double result = 0.00; OpenMesh::Vec3f CurN; double Cur_diff;
	int num_increase = 0, num_decrease = 0;
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		CurN[0] = m_vertices[v_it.handle().idx()].ps_normal_x; CurN[1] = m_vertices[v_it.handle().idx()].ps_normal_y; CurN[2] = m_vertices[v_it.handle().idx()].ps_normal_z;
		Cur_diff = (T_Mesh.normal(v_it.handle()) - CurN).sqrnorm()*(UseFaceArea?m_vertices[v_it.handle().idx()].BCDArea:1.0);	
		result += Cur_diff;
		if (vec_pervertex_energy[v_it.handle().idx()] <= Cur_diff) {
			num_increase ++;
			//T_Mesh.set_color(v_it, MyMesh::Color(255,0,0));
		} else {
			num_decrease ++;
			//T_Mesh.set_color(v_it, MyMesh::Color(0,255,0));
		}
		vec_pervertex_energy[v_it.handle().idx()] = Cur_diff;
	}
	//cout << "For per vertex normal energy: " << num_increase << " increase, " << num_decrease << " decrease." << endl;
	return result;
}

double TriangularMesh::CalculateNDPEnergy(MyMesh &T_Mesh, double varsigma = 200, bool UseFaceArea = false)
{
	double result = 0.0;
	T_Mesh.update_face_normals(); T_Mesh.update_vertex_normals();
	double max_omega = 1, min_omega = 0;
	max_omega = -1; min_omega = 10000;
	for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
		int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
		int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
		double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);

#ifdef TEST_MESHREFINE
		OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
		OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
		omegaij = std::exp(-(psni-psnj).norm());
#endif
		if (omegaij < min_omega) {
			min_omega = omegaij;
		}
		if (omegaij > max_omega) {
			max_omega = omegaij;
		}
	}
	for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
		int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
		int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
		double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
		omegaij = (omegaij-min_omega)/(max_omega-min_omega);

#ifdef TEST_MESHREFINE
		OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
		OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
		omegaij = (std::exp(-(psni-psnj).norm()) - min_omega)/(max_omega - min_omega);
#endif
		result += omegaij*(T_Mesh.normal(MyMesh::VertexHandle(iid)) - T_Mesh.normal(MyMesh::VertexHandle(jid))).sqrnorm() * (UseFaceArea?m_vertices[iid].BCDArea:1.0); 
	}
	return result;
}

double TriangularMesh::CalculateNDFEnergy(MyMesh &T_Mesh, double varsigma, bool UseFaceArea = false)
{
	double max_omega = 1, min_omega = 0;
	max_omega = -1; min_omega = 10000;
	for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
		MyMesh::FaceHandle fh1 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,0)); 
		MyMesh::FaceHandle fh2 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,1));
		int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
		int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();

		MyMesh::ConstFaceVertexIter cfv_it = T_Mesh.cfv_iter(fh1);
		while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
			++cfv_it;
		}  int k1id = cfv_it.handle().idx();

		cfv_it = T_Mesh.cfv_iter(fh2);
		while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
			++cfv_it;
		}  int k2id = cfv_it.handle().idx();
		double omegaij = std::exp(-pow((m_vertices[k1id].intensity - m_vertices[k2id].intensity), 2.0)/varsigma);

#ifdef TEST_MESHREFINE
		OpenMesh::Vec3f psni(m_vertices[k1id].ps_normal_x,m_vertices[k1id].ps_normal_y,m_vertices[k1id].ps_normal_z);
		OpenMesh::Vec3f psnj(m_vertices[k2id].ps_normal_x,m_vertices[k2id].ps_normal_y,m_vertices[k2id].ps_normal_z);
		omegaij = std::exp(-(psni-psnj).norm());
#endif
		if (omegaij < min_omega) {
			min_omega = omegaij;
		}
		if (omegaij > max_omega) {
			max_omega = omegaij;
		}
	}
	double result = 0.0;
	T_Mesh.update_face_normals(); T_Mesh.update_vertex_normals();
	for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
		int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
		int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
		MyMesh::FaceHandle fh1 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,0)); 
		MyMesh::FaceHandle fh2 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,1));

		MyMesh::ConstFaceVertexIter cfv_it = T_Mesh.cfv_iter(fh1);
		while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
			++cfv_it;
		} int k1id = cfv_it.handle().idx();

		cfv_it = T_Mesh.cfv_iter(fh2);
		while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
			++cfv_it;
		} int k2id = cfv_it.handle().idx();
		double omegaij = std::exp(-pow((m_vertices[k1id].intensity - m_vertices[k2id].intensity), 2.0)/varsigma);
		omegaij = (omegaij-min_omega)/(max_omega-min_omega);

#ifdef TEST_MESHREFINE
		OpenMesh::Vec3f psni(m_vertices[k1id].ps_normal_x,m_vertices[k1id].ps_normal_y,m_vertices[k1id].ps_normal_z);
		OpenMesh::Vec3f psnj(m_vertices[k2id].ps_normal_x,m_vertices[k2id].ps_normal_y,m_vertices[k2id].ps_normal_z);
		omegaij = (std::exp(-(psni-psnj).norm()) - min_omega)/(max_omega - min_omega);
#endif
		result += omegaij*(T_Mesh.normal(fh1) - T_Mesh.normal(fh2)).sqrnorm() * (UseFaceArea?m_trianglesArea[fh1.idx()]:1.0); // need to times some var, e.g. edge length, tri area
	}
	return result;
}

double TriangularMesh::CalculateVTVEnergy(MyMesh &T_Mesh, vector<double> &ux, vector<double> &uy, vector<double> &uz, bool UseFaceArea = false)
{
	GetPPIGradients (m_vnum,ux);
	GetPPIGradients1(m_vnum,uy);
	GetPPIGradients2(m_vnum,uz);
	double result = 0.0;
	for (int i = 0; i < m_trinum; ++ i) {
		result += (NormSquare(m_triangles[i].grad) + NormSquare(m_triangles[i].grad1) + NormSquare(m_triangles[i].grad2))*(UseFaceArea?m_trianglesArea[i]:1.0);
	}
	return result;
}

double TriangularMesh::CalculateVTVNormEnergy(MyMesh &T_Mesh, vector<double> &nux, vector<double> &nuy, vector<double> &nuz, bool UseFaceArea = false)
{
	return CalculateVTVEnergy(T_Mesh, nux, nuy, nuz, UseFaceArea);
}

double TriangularMesh::CalculateLaplaceEnergy(MyMesh &T_Mesh, bool UseFaceArea = false)
{
	double result = 0.0;
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		int vertex_id = v_it.handle().idx();

		double degree = 0.0; OpenMesh::Vec3f avg_nei_point; avg_nei_point.vectorize(0.0);
		for (MyMesh::VertexVertexIter vv_it = T_Mesh.vv_iter(v_it); vv_it; ++vv_it) {
			degree++;
			avg_nei_point += T_Mesh.point(vv_it.handle());
		} avg_nei_point = avg_nei_point/degree;

		result += (T_Mesh.point(v_it.handle()) - avg_nei_point).sqrnorm()*(UseFaceArea?m_vertices[v_it.handle().idx()].BCDArea:1.0);
	}
	return result;
}

void TriangularMesh::ALM_TVU_MeshRefinement(string meshname, double fidParam, double pld_eta, double pcd_eta, double fcd_eta, double pnd_eta, double fnd_eta, double varsigma, double pc_eta, double penParam, double regParam = 1.0, double lapParam = 100, bool UseTVU = true, bool UseTVNorm = false, int iter_step = 0, bool UseFaceArea = false, bool UseMatlabSolver = true)
{
	unsigned long nver = m_vnum, ntri = m_trinum; int itmax = 300, itol = 1;//1,2,3, or 4.
	double tol = 1.0e-15;  char buffer[255]; 
	vector<double>ux, uy, uz; ux.resize(nver); uy.resize(nver); uz.resize(nver);
	vector<double>nux, nuy, nuz; nux.resize(nver); nuy.resize(nver); nuz.resize(nver);
	vector<double>uxold, uyold, uzold; uxold.resize(nver); uyold.resize(nver); uzold.resize(nver);
	vector<VECTOR3D> px, py, pz; px.resize(ntri); py.resize(ntri); pz.resize(ntri);
	vector<VECTOR3D> lambda_x, lambda_y, lambda_z; lambda_x.resize(ntri); lambda_y.resize(ntri); lambda_z.resize(ntri);

	/*   Parameters setting  */
	unsigned int innerL = 1; int outL = -1;
	double outTole = 1.0e-12, stoppingCond;
	MyMesh T_Mesh = this->m_ObjTriMesh; OpenMesh::Vec3f CurV; 
	OpenMesh::IO::Options write_options;
	if (RecordColor) {
		write_options.set(OpenMesh::IO::Options::VertexColor); 
	}
	T_Mesh.request_face_normals();		T_Mesh.update_face_normals();
	T_Mesh.request_vertex_normals();	T_Mesh.update_vertex_normals();
	bool UsePositionFidelity = fidParam>0?true:false;	bool UsePointLightDiff = (pld_eta>0)?true:false;	
	bool UsePointColor = pc_eta>0?true:false;			bool UsePointColorDiff = pcd_eta>0?true:false;		bool UseFaceColorDiff = fcd_eta>0?true:false;
	bool UsePointNormalDiff = pnd_eta>0?true:false;		bool UseFaceNormalDiff = fnd_eta>0?true:false;		bool UseLaplace = (lapParam>0)?true:false; 
	cout << endl << "Start the mesh refinement process: ";		UsePositionFidelity?cout<<"Pfid "<<fidParam<<" ":cout<<" ";
	UsePointLightDiff?cout<<"Pld "<<pld_eta<<" ":cout<<" ";		UsePointColor?cout<<"PCol "<<pc_eta<<" ":cout<<" ";	
	UsePointColorDiff?cout<<"PCdiff "<<pcd_eta<<" ":cout<<" ";	UseFaceColorDiff?cout<<"FCdiff "<<fcd_eta<<" ":cout<<" "; 
	UsePointNormalDiff?cout<<"PND "<<pnd_eta<<" ":cout<<" ";	UseFaceNormalDiff?cout<<"FND "<<fnd_eta<<" ":cout<<" ";	
	UseLaplace?cout<<"Lap "<<lapParam<<" ":cout<<" ";
	UseTVU?cout<<"TVU "<<penParam<<" ":cout<<" ";				UseTVNorm?cout<<"TVNorm "<<penParam<<" ":cout<<" ";cout << endl;

	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[0]; 
		uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[1];
		uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[2];
		m_vertices[v_it.handle().idx()].light_x = m_vertices[v_it.handle().idx()].light_y = m_vertices[v_it.handle().idx()].light_z = 0.0;
	}
	if (UseTVNorm) {
		T_Mesh.update_face_normals();  T_Mesh.update_vertex_normals();
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
			nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
			nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
			nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
		}
		GetPPIGradients (m_vnum,nux);
		GetPPIGradients1(m_vnum,nuy);
		GetPPIGradients2(m_vnum,nuz);
	} 
	if (UseTVU) {
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
			uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
			uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
		}
		GetPPIGradients (nver,ux);
		GetPPIGradients1(nver,uy);
		GetPPIGradients2(nver,uz);
	}
	for(int i=0;i<ntri;i++)
	{
		px[i].x=m_triangles[i].grad.x;	lambda_x[i].x=0;	px[i].y=m_triangles[i].grad.y;	lambda_x[i].y=0;	px[i].z=m_triangles[i].grad.z;	lambda_x[i].z=0;
		py[i].x=m_triangles[i].grad1.x; lambda_y[i].x=0;	py[i].y=m_triangles[i].grad1.y; lambda_y[i].y=0;	py[i].z=m_triangles[i].grad1.z; lambda_y[i].z=0;
		pz[i].x=m_triangles[i].grad2.x; lambda_z[i].x=0;	pz[i].y=m_triangles[i].grad2.y; lambda_z[i].y=0;	pz[i].z=m_triangles[i].grad2.z; lambda_z[i].z=0;
	}
	vector<vector<double>> energy_result; energy_result.clear();
	vector<double> vec_n_energy; vec_n_energy.resize(m_vnum); 
	for (int i = 0; i < vec_n_energy.size(); ++ i) {
		vec_n_energy[i] = 0.0;
	}
	double v_energy, n_energy, ndp_energy, ndf_energy, tv_energy, lap_energy, total_energy;   
	do { // start the iteration
		cout << "Main iteration step: " << ++outL << endl;
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			uxold[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[0]; 
			uyold[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[1];
			uzold[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[2];
		}
		// inner iteration
		for (int iter = 0; iter < innerL; ++ iter) {
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				CurV[0] = ux[v_it.handle().idx()]; CurV[1] = uy[v_it.handle().idx()]; CurV[2] = uz[v_it.handle().idx()];
				T_Mesh.set_point(v_it.handle(), CurV);
			}

			// solve the u-sub problem, using Levenberg–Marquardt Algorithm
			double sigma, eOld, eNew;
			double var = 2.0, epsilon1 = 1.0e-12, epsilon2 = 1.0e-12, avgEdgeLength = 0.0, minEdgeLength = 10000;
			for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
				int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
				int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
				double edge_length = (T_Mesh.point(MyMesh::VertexHandle(iid)) - T_Mesh.point(MyMesh::VertexHandle(jid))).length();
				avgEdgeLength += edge_length;
				if (minEdgeLength > edge_length) {
					minEdgeLength = edge_length;
				}
			}
			avgEdgeLength = avgEdgeLength/T_Mesh.n_edges();    //avgEdgeLength = minEdgeLength;

			int dimension = m_vnum*3 + m_vnum*3;
			DenseMatrix m_b(dimension, 1);
			DenseMatrix m_x(dimension, 1);
			std::vector<double> tau(dimension, 0);
			std::vector<OpenMesh::Vec3f> originalPosition(m_vnum);
			std::vector<OpenMesh::Vec3f> originalLight(m_vnum);
			RowSparseMatrix mat_J;	DenseMatrix m_f;
			gmm::csc_matrix<double> JtJ;

			TV_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, fidParam, pld_eta, pcd_eta, fcd_eta, pnd_eta, fnd_eta, varsigma, pc_eta, penParam, lapParam,
				px, py, pz, lambda_x, lambda_y, lambda_z, UseTVU, UseTVNorm, UseFaceArea);

			eOld = gmm::mat_euclidean_norm_sqr(m_f);
			vector<double> t_energy; t_energy.clear(); cout <<"00:";
			CalculateEnergyTerm(m_f, t_energy, UsePositionFidelity, UsePointLightDiff, UsePointColor, UsePointColorDiff, UseFaceColorDiff, 
				UsePointNormalDiff, UseFaceNormalDiff, UseLaplace, UseTVNorm, UseTVU);		energy_result.push_back(t_energy); 
			sprintf(buffer, "Results\\LMResults\\%sLMResult%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s_%d_%d.off", 
				meshname.c_str(), iter_step, UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
				UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
				UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
				UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_", outL, 0);
			if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
				std::cerr << "Cannot write mesh to file " << buffer << std::endl;
			}

			RowSparseMatrix mat_JTJ(dimension, dimension);  RowSparseMatrix mat_A(dimension, dimension);
			gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
			for (int i = 0; i < dimension; i++) {
				tau[i] = mat_JTJ(i, i);
			}
			gmm::scale(tau, 1.e-5);
			gmm::mult(gmm::transposed(mat_J), m_f, m_b);
			gmm::scale(m_b, -1.0);
			double normG = gmm::mat_norminf(m_b); 

			numc::SparseSolver solver;	double *b = new double[dimension];	double *vx = new double[dimension];
			int iteration = 0; int maxStep = 40;
			if (UseTVU || UseTVNorm) {
				maxStep = ceil(maxStep/2.0);
			}
			while (iteration++ < maxStep) {
				sprintf(buffer, "%.2d:", iteration);
				cout << endl << buffer;
				if (normG < epsilon1) {
					break;
				}
				gmm::copy(mat_JTJ, mat_A);
				for (int i = 0; i < dimension; i++) {
					mat_A(i, i) += tau[i];
				}

				if (UseMatlabSolver) {
					gmm::copy(gmm::transposed(mat_A), JtJ);
					engEvalString( m_ep, "clear all;" );
					int nnz= gmm::nnz(JtJ);
					mxArray* arrayA=mxCreateSparse(JtJ.nrows(),JtJ.ncols(),nnz,mxREAL);
					double*  pr = mxGetPr( arrayA );
					mwIndex* ir = mxGetIr( arrayA );
					mwIndex* jc = mxGetJc( arrayA );
					int i=0;
					for(std::vector<double>::const_iterator it=JtJ.pr.begin(); it!=JtJ.pr.end();++it,i++)
						pr[i]=(*it);
					i=0;
					for(std::vector<unsigned int>::const_iterator it=JtJ.ir.begin(); it!=JtJ.ir.end();++it,i++)
						ir[i]=(*it);
					i=0;
					for(std::vector<unsigned int>::const_iterator it=JtJ.jc.begin(); it!=JtJ.jc.end();++it,i++)
						jc[i]=(*it);

					engPutVariable( m_ep, "A", arrayA );
					mxDestroyArray(arrayA);

					mxArray* arrayb=mxCreateDoubleMatrix(m_b.size(),1,mxREAL);
					double*  prb = mxGetPr( arrayb );
					for( int r=0;r<m_b.size();++r)
						prb[r]= m_b[r];
					engPutVariable( m_ep, "b", arrayb );
					mxDestroyArray(arrayb);

					engEvalString( m_ep, "x = A\\b;" );

					mxArray* arrayx = engGetVariable( m_ep, "x");
					int m = mxGetDimensions( arrayx )[0];
					int n = mxGetDimensions( arrayx )[1];
					if( !arrayx || mxIsSparse( arrayx ) || mxIsComplex( arrayx ) || n!=1 )
					{
						mxDestroyArray( arrayx );
						cout << "Get x error"<<endl;
					}
					int mM= m_x.size();
					if( m != mM ){
						m_x.resize(m,1);
					}
					double*  prx = mxGetPr( arrayx );
					for( int r=0;r<m;++r)
						m_x(r,0)=prx[r];
				} 
				else {
					// solve mat_A*m_x = m_b
					numc::RowMat<double> RM_A(dimension, dimension); int r = 0;
					for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(mat_A);
						rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(mat_A); ++ rIt, ++ r) {
							gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
							gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
							for (; it != ite; ++ it) {
								RM_A(r, it.index()) = mat_A(r, it.index());
							}
							b[r] = m_b[r];
					}
					solver.getMatA() = RM_A;
					solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
					solver.init();
					// solve the nonlinear equation to calculate the new vertex position in vx;
					solver.solve(b, vx);
					solver.clear();
					for (int r = 0; r < dimension; ++ r) {
						m_x(r,0) = vx[r];
					}
				}

				if (ScaleDelta) {
					//modify m_x, do not change too much according to the value of MinEdgeLength
					double maxChange = 0.0;
					for (int i = 0; i < m_vnum; ++ i) {
						OpenMesh::Vec3f deltaP(m_x[i+m_vnum*0], m_x[i+m_vnum*1], m_x[i+m_vnum*2]);
						maxChange = (deltaP.length() > maxChange)?deltaP.length():maxChange;
					}
					//cout << endl << "The max delta p is " << maxChange/avgEdgeLength << " of avg edge length" << endl;
					double stepSize = (avgEdgeLength/maxChange > 1.0)?1.0:(avgEdgeLength/maxChange);
					gmm::scale(m_x, stepSize);
					//end of modify m_x
				}

				double normH = 0.0, normX = 0.0, l0h = 0;
				normH = gmm::mat_euclidean_norm(m_x); 
				for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
					int idx = v_it.handle().idx();
					originalLight[v_it.handle().idx()] = 
						OpenMesh::Vec3f(m_vertices[v_it.handle().idx()].light_x, m_vertices[v_it.handle().idx()].light_y, m_vertices[v_it.handle().idx()].light_z);
					originalPosition[v_it.handle().idx()] = T_Mesh.point(v_it.handle());
					normX += originalPosition[v_it.handle().idx()].sqrnorm() + originalLight[v_it.handle().idx()].sqrnorm();
					for (int i = 0; i < 6; ++i) {
						l0h += 0.5*m_x[i*m_vnum + idx]*(tau[i*m_vnum + idx]*m_x[i*m_vnum + idx] + m_b[i*m_vnum + idx]);
					}
				}
				if (normH < epsilon2*(sqrt(normX) + epsilon2)) {
					break;
				} else {
					//change the coordinate value of each vertex...
					for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
						int vIndex = v_it.handle().idx();
						T_Mesh.point(v_it.handle()) += OpenMesh::Vec3f(m_x[vIndex+m_vnum*0], m_x[vIndex+m_vnum*1], m_x[vIndex+m_vnum*2]);
						m_vertices[vIndex].light_x += m_x[vIndex+m_vnum*0+m_vnum*3];
						m_vertices[vIndex].light_y += m_x[vIndex+m_vnum*1+m_vnum*3];
						m_vertices[vIndex].light_z += m_x[vIndex+m_vnum*2+m_vnum*3];
					}
					UpdateMeshInfo(T_Mesh); 

					TV_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, fidParam, pld_eta, pcd_eta, fcd_eta, pnd_eta, fnd_eta, varsigma, pc_eta, penParam, lapParam,
						px, py, pz, lambda_x, lambda_y, lambda_z, UseTVU, UseTVNorm, UseFaceArea);
					eNew = gmm::mat_euclidean_norm_sqr(m_f); 

					vector<double> t_energy; t_energy.clear(); 
					CalculateEnergyTerm(m_f, t_energy, UsePositionFidelity, UsePointLightDiff, UsePointColor, UsePointColorDiff, UseFaceColorDiff, 
						UsePointNormalDiff, UseFaceNormalDiff, UseLaplace, UseTVNorm, UseTVU);		energy_result.push_back(t_energy); 
					//v_energy = UseFidelity?CalculateVEnergy(T_Mesh,UseFaceArea):0.0;				t_energy.push_back(0.5*fidParam*v_energy);
					//n_energy = UseGDN?CalculateNEnergy(T_Mesh, vec_n_energy,UseFaceArea):0.0;		t_energy.push_back(0.5*pld_eta*n_energy); 
					//ndp_energy = UsePointND?CalculateNDPEnergy(T_Mesh, varsigma,UseFaceArea):0.0;	t_energy.push_back(0.5*pcd_eta*ndp_energy);
					//ndf_energy = UseFaceND?CalculateNDFEnergy(T_Mesh, varsigma,UseFaceArea):0.0;	t_energy.push_back(0.5*fcd_eta*ndf_energy);
					//lap_energy = UseLaplace?CalculateLaplaceEnergy(T_Mesh,UseFaceArea):0.0;		t_energy.push_back(0.5*lapParam*lap_energy);
					//tv_energy = 0.0;
					//if (UseTVNorm) {
					//	T_Mesh.update_face_normals();  T_Mesh.update_vertex_normals();
					//	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
					//		nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
					//		nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
					//		nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
					//	}
					//	tv_energy = CalculateVTVNormEnergy(T_Mesh, nux, nuy, nuz,UseFaceArea);		t_energy.push_back(0.5*penParam*tv_energy);
					//}
					//if (UseTVU) {
					//	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
					//		ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[0];
					//		uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[1];
					//		uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[2];
					//	}
					//	tv_energy = CalculateVTVEnergy(T_Mesh, ux, uy, uz,UseFaceArea);				t_energy.push_back(0.5*penParam*tv_energy);
					//}
					if (iteration%1 == 0) {
						sprintf(buffer, "Results\\LMResults\\%sLMResult%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s_%d_%d.off", 
							meshname.c_str(), iter_step, UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
							UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
							UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
							UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_", outL, iteration);
						if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
							std::cerr << "Cannot write mesh to file " << buffer << std::endl;
						}
						sprintf(buffer, "Results\\LMResults\\%sLMResult%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s_%d_%d.txt", 
							meshname.c_str(), iter_step, UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
							UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
							UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
							UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_", outL, iteration);
						fstream fout(buffer,ios::out);
						for (int i=0;i<m_vnum;++i) {
							fout<<m_vertices[i].light_x<<", "<<m_vertices[i].light_y<<", "<<m_vertices[i].light_z<<endl;
						}
						fout.close();
					}

					sigma = (eOld - eNew)/(l0h);
					if (sigma >= 0) {
						gmm::scale(tau, std::max(1./3., 1. - pow(2.*sigma - 1., 3.)));
						var = 2.0;
						eOld = eNew;
						gmm::mult(gmm::transposed(mat_J), m_f, m_b);
						gmm::scale(m_b, -1.0);
						normG = gmm::mat_norminf(m_b);
						gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
					} else {
						cout << " Sigma = " << sigma;// << "; eOld = " << eOld << "; eNew = " << eNew;
						gmm::scale(tau, var);
						var *= 2.0;  
						for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
							T_Mesh.point(v_it) = originalPosition[v_it.handle().idx()];
							m_vertices[v_it.handle().idx()].light_x = originalLight[v_it.handle().idx()].data()[0];
							m_vertices[v_it.handle().idx()].light_y = originalLight[v_it.handle().idx()].data()[1];
							m_vertices[v_it.handle().idx()].light_z = originalLight[v_it.handle().idx()].data()[2];
						}
					}
				}
			}
			if (!(UseTVU||UseTVNorm)) {
				sprintf(buffer, "Results\\ALM_%sResult_Final%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s.off", 
					meshname.c_str(), iter_step,  UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
					UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
					UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
					UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_");
				if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
					std::cerr << "Cannot write mesh to file " << buffer << std::endl;
				}
				sprintf(buffer, "Results\\%sEnergyResult%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s.off", 
					meshname.c_str(), iter_step, UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
					UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
					UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
					UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_");
				fstream nof(buffer,std::ios::out);
				if (!nof) {
					cout << "Can not open "<<buffer<<" to save data..." << endl;
					return;
				}
				for (int i = 0; i < energy_result.size(); ++ i) {
					for (int j = 0; j < energy_result[i].size(); ++ j) {
						nof << setiosflags(ios::fixed)<<setprecision(6)<<setw(10)<<setiosflags(ios::left) << energy_result[i][j] << "  ";
					}
					nof << endl;
				}
				nof.close();
				return;
			}
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
				uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
				uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
			}

			/*  Solve the p-sub problem  (solve for px,py,pz).   */

			UpdateMeshInfo(T_Mesh);
			if (UseTVNorm) {
				T_Mesh.update_face_normals();  T_Mesh.update_vertex_normals();
				for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
					nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
					nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
					nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
				}
				GetPPIGradients (m_vnum,nux);
				GetPPIGradients1(m_vnum,nuy);
				GetPPIGradients2(m_vnum,nuz);
			} 
			if (UseTVU) {
				for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
					ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
					uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
					uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
				}
				GetPPIGradients (nver,ux);
				GetPPIGradients1(nver,uy);
				GetPPIGradients2(nver,uz);
			}
			VECTOR3D wx, wy, wz;
			double pdiff = 0.0, avg_wxy = 0.0; int i1 = 0, i2 = 0;
			for(int i=0;i<ntri;i++)
			{
				wx.x = m_triangles[i].grad.x - lambda_x[i].x/penParam;
				wx.y = m_triangles[i].grad.y - lambda_x[i].y/penParam;
				wx.z = m_triangles[i].grad.z - lambda_x[i].z/penParam;

				wy.x = m_triangles[i].grad1.x - lambda_y[i].x/penParam;
				wy.y = m_triangles[i].grad1.y - lambda_y[i].y/penParam;
				wy.z = m_triangles[i].grad1.z - lambda_y[i].z/penParam;

				wz.x = m_triangles[i].grad2.x - lambda_z[i].x/penParam;
				wz.y = m_triangles[i].grad2.y - lambda_z[i].y/penParam;
				wz.z = m_triangles[i].grad2.z - lambda_z[i].z/penParam;
				avg_wxy += sqrt(NormSquare(wx)+NormSquare(wy)+NormSquare(wz));

				if(DotProduct(wx,wx)+DotProduct(wy,wy)+DotProduct(wz,wz) <= POWER(regParam/penParam)) {
					px[i].x = 0; py[i].x = 0; pz[i].x = 0;
					px[i].y = 0; py[i].y = 0; pz[i].y = 0;
					px[i].z = 0; py[i].z = 0; pz[i].z = 0;
					i1 ++;
				} else {
					double tempNorm; i2++;
					tempNorm = sqrt(DotProduct(wx,wx)+DotProduct(wy,wy)+DotProduct(wz,wz));
					px[i].x = (1-regParam/penParam/tempNorm)*wx.x; py[i].x = (1-regParam/penParam/tempNorm)*wy.x; pz[i].x = (1-regParam/penParam/tempNorm)*wz.x;
					px[i].y = (1-regParam/penParam/tempNorm)*wx.y; py[i].y = (1-regParam/penParam/tempNorm)*wy.y; pz[i].y = (1-regParam/penParam/tempNorm)*wz.y;
					px[i].z = (1-regParam/penParam/tempNorm)*wx.z; py[i].z = (1-regParam/penParam/tempNorm)*wy.z; pz[i].z = (1-regParam/penParam/tempNorm)*wz.z;
				}
				pdiff += Norm(px[i]) + Norm(py[i]) + Norm(pz[i]);
			}
			cout << endl << "Solve p sub-problem: Pdiff is: " << pdiff << "  " << i1 << "/" << i2 << " under " << avg_wxy/ntri << endl;
		}
		/*   Update Lagrange multipliers  (lambda_x,lambda_y,lambda_z).      */
		cout << "Update lambda: ";
		if (UseTVNorm) {
			T_Mesh.update_face_normals();  T_Mesh.update_vertex_normals();
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
				nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
				nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
				nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
			}
			GetPPIGradients (m_vnum,nux);
			GetPPIGradients1(m_vnum,nuy);
			GetPPIGradients2(m_vnum,nuz);
		} 
		if (UseTVU) {
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
				uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
				uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
			}
			GetPPIGradients (nver,ux);
			GetPPIGradients1(nver,uy);
			GetPPIGradients2(nver,uz);
		}
		for(int i=0;i<ntri;i++)
		{
			lambda_x[i].x+= penParam*(px[i].x-m_triangles[i].grad.x);
			lambda_x[i].y+= penParam*(px[i].y-m_triangles[i].grad.y);
			lambda_x[i].z+= penParam*(px[i].z-m_triangles[i].grad.z);

			lambda_y[i].x+= penParam*(py[i].x-m_triangles[i].grad1.x);
			lambda_y[i].y+= penParam*(py[i].y-m_triangles[i].grad1.y);
			lambda_y[i].z+= penParam*(py[i].z-m_triangles[i].grad1.z);

			lambda_z[i].x+= penParam*(pz[i].x-m_triangles[i].grad2.x);
			lambda_z[i].y+= penParam*(pz[i].y-m_triangles[i].grad2.y);
			lambda_z[i].z+= penParam*(pz[i].z-m_triangles[i].grad2.z);
		}
		/*   Compute the stopping condition     */
		stoppingCond = 0.0;
		for(int i=0;i<nver;i++){
			stoppingCond+= (POWER(ux[i]-uxold[i])+POWER(uy[i]-uyold[i])+POWER(uz[i]-uzold[i])) * m_vertices[i].BCDArea;
		}
		cout << "The error is: " << setprecision(15) <<setw(10) << stoppingCond << "; the threshold is: " << setprecision(15) <<setw(10) << outTole << endl;
	} while (stoppingCond>outTole);
	this->m_ObjTriMesh = T_Mesh;
	for(int i=0;i<nver;i++){
		m_vertices[i].x = ux[i]; m_vertices[i].y = uy[i]; m_vertices[i].z = uz[i];
		m_vertices[i].ref_x = ux[i]; m_vertices[i].ref_y = uy[i]; m_vertices[i].ref_z = uz[i];
	}
	sprintf(buffer, "Results\\ALM_%sResult_Final%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s.off", 
		meshname.c_str(), iter_step, UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
		UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
		UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
		UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_");
	if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
		std::cerr << "Cannot write mesh to file " << buffer << std::endl;
	}
	sprintf(buffer, "Results\\%sEnergyResult%d_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%.0f_%s%s%.5f_%.0f_%.0f_%s.off", 
		meshname.c_str(), iter_step, UsePositionFidelity?"Fid":"_", fidParam, UsePointLightDiff?"PLD":"_", pld_eta, 
		UsePointColor?"PC":"_", pc_eta, UsePointColorDiff?"PCD":"_", pcd_eta, UseFaceColorDiff?"FCD":"_", fcd_eta, 
		UsePointNormalDiff?"PND":"_", pnd_eta, UseFaceNormalDiff?"FND":"_", fnd_eta, UseLaplace?"Lap":"_", lapParam, 
		UseTVU?"TVU":"_",UseTVNorm?"TVNorm":"_", penParam, regParam, varsigma, UseFaceArea?"FA":"_");
	fstream nof(buffer,std::ios::out);
	if (!nof) {
		cout << "Can not open "<<buffer<<" to save data..." << endl;
		return;
	}
	for (int i = 0; i < energy_result.size(); ++ i) {
		for (int j = 0; j < energy_result[i].size(); ++ j) {
			nof << setiosflags(ios::fixed)<<setprecision(6)<<setw(10)<<setiosflags(ios::left) << energy_result[i][j] << "  ";
		}
		nof << endl;
	}
	nof.close();
	//delete b; delete vx;
	return;
}

void TriangularMesh::CalculateEnergyTerm(DenseMatrix mat_f, vector<double>& energy, bool UsePositionFidelity, bool UsePointLightDiff, bool UsePointColor, bool UsePointColorDiff, bool UseFaceColorDiff, bool UsePointNormalDiff, bool UseFaceNormalDiff, bool UseLaplace, bool UseTVNorm, bool UseTVU)
{
	energy.clear(); int Start_id = 0; DenseMatrix temp;
	if (UsePositionFidelity) {
		temp.resize(m_vnum*3,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < m_vnum*3; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " UFID: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp); 
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += m_vnum*3;
	}
	if (UsePointLightDiff) {
		temp.resize(this->m_ObjTriMesh.n_edges()*3,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < this->m_ObjTriMesh.n_edges()*3; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " PLD: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += this->m_ObjTriMesh.n_edges()*3;
	}
	if (UsePointColor) {
		temp.resize(m_vnum,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < m_vnum; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " PCOL: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += m_vnum;
	}
	if (UsePointColorDiff) {
		temp.resize(this->m_ObjTriMesh.n_edges(),1); gmm::scale(temp, 0.0);
		for (int i = 0; i < this->m_ObjTriMesh.n_edges(); ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " PCD: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += this->m_ObjTriMesh.n_edges();
	}
	if (UseFaceColorDiff) {
		temp.resize(this->m_ObjTriMesh.n_edges(),1); gmm::scale(temp, 0.0);
		for (int i = 0; i < this->m_ObjTriMesh.n_edges(); ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " FCD: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += this->m_ObjTriMesh.n_edges();
	}
	if (UsePointNormalDiff) {
		temp.resize(this->m_ObjTriMesh.n_edges()*3,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < this->m_ObjTriMesh.n_edges()*3; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " PND: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += this->m_ObjTriMesh.n_edges()*3;
	}
	if (UseFaceNormalDiff) {
		temp.resize(this->m_ObjTriMesh.n_edges()*3,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < this->m_ObjTriMesh.n_edges()*3; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " FND: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += this->m_ObjTriMesh.n_edges()*3;
	}
	if (UseLaplace) {
		temp.resize(m_vnum*3,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < m_vnum*3; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " LAP: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += m_vnum*3;
	}
	if (UseTVNorm) {
		temp.resize(m_trinum*9,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < m_trinum*9; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " TVN: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += UseTVNorm*m_trinum*9;
	}
	if (UseTVU) {
		temp.resize(m_trinum*9,1); gmm::scale(temp, 0.0);
		for (int i = 0; i < m_trinum*9; ++i) {
			temp(i,0) = mat_f(Start_id+i,0);
		}
		cout << " TVU: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< gmm::mat_euclidean_norm_sqr(temp);
		energy.push_back(gmm::mat_euclidean_norm_sqr(temp));
		Start_id += UseTVNorm*m_trinum*9;
	}
	double total_energy = 0.0;
	for (int k = 0; k < energy.size(); ++ k) {
		total_energy += energy[k];
	}
	energy.push_back(total_energy);
	cout << " Sum: " <<setprecision(6)<<setw(6)<<setiosflags(ios::left)<< energy[energy.size()-1];
}

void TriangularMesh::TV_JacobianMatrix_Construction(MyMesh& T_Mesh, RowSparseMatrix& mat_J, DenseMatrix& mat_f, double fidParam, double pld_eta, 
	double pcd_eta, double fcd_eta, double pnd_eta, double fnd_eta, double varsigma, double pc_eta, double penParam, double lapParam, 
	vector<VECTOR3D> &px, vector<VECTOR3D> &py, vector<VECTOR3D> &pz, vector<VECTOR3D> &lambda_x, vector<VECTOR3D> &lambda_y, vector<VECTOR3D> &lambda_z, 
	bool UseTVU = false, bool UseTVNorm = false, bool UseFaceArea = false)
{
	double epsilon = 1.0e-5; UpdateMeshInfo(T_Mesh);
	bool UsePositionFidelity = fidParam>0?true:false;	bool UsePointLightDiff = pld_eta>0?true:false;			// fidelity terms
	bool UsePointColor = pc_eta>0?true:false;			bool UsePointColorDiff = pcd_eta>0?true:false;		bool UseFaceColorDiff = fcd_eta>0?true:false;
	bool UsePointNormalDiff = pnd_eta>0?true:false;		bool UseFaceNormalDiff = fnd_eta?true:false;		bool UseLaplace = lapParam>0?true:false;
	// construct the jacobian matrix according to the new model
	int dimension = UsePositionFidelity*m_vnum*3 + UsePointLightDiff*T_Mesh.n_edges()*3 + UsePointColor*m_vnum + 
		UsePointColorDiff*T_Mesh.n_edges() + UseFaceColorDiff*T_Mesh.n_edges() + UsePointNormalDiff*T_Mesh.n_edges()*3 + UseFaceNormalDiff*T_Mesh.n_edges()*3 + 
		UseLaplace*m_vnum*3 + UseTVNorm*m_trinum*9 + UseTVU*m_trinum*9;
	mat_J.resize(dimension, m_vnum*3+m_vnum*3); mat_f.resize(dimension, 1);
	gmm::scale(mat_J, 0.0);				gmm::scale(mat_f, 0.0);
	double fid_scale = sqrt(fidParam*0.5), pld_scale = sqrt(pld_eta*0.5),	pcd_scale = sqrt(pcd_eta*0.5),	fcd_scale = sqrt(fcd_eta*0.5), 
		   pen_scale = sqrt(penParam*0.5), lap_scale = sqrt(lapParam*0.5),	pnd_scale = sqrt(pnd_eta*0.5),	fnd_scale = sqrt(fnd_eta*0.5),		
		   pc_scale = sqrt(pc_eta*0.5);

	RowSparseMatrix gradJ; std::vector<double> grad_f;
	this->LM_JacobianMatrix_Construction(T_Mesh, gradJ, grad_f, 0);

	int Start_id = 0;
	if (UsePositionFidelity) {
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			int vertex_id = v_it.handle().idx();											double area_scale = UseFaceArea?sqrt(m_vertices[vertex_id].BCDArea):1.0;
			// the fidelity term, \alpha/2 \|u-f\|^2
			OpenMesh::Vec3f ref_p = OpenMesh::Vec3f(m_vertices[vertex_id].ref_x, m_vertices[vertex_id].ref_y,m_vertices[vertex_id].ref_z);
			for (int k = 0; k < 3; ++ k) {
				mat_J(vertex_id*3+k+Start_id, vertex_id+m_vnum*k) += area_scale*fid_scale*1.0;
				mat_f(vertex_id*3+k+Start_id, 0) += area_scale*fid_scale*(T_Mesh.point(v_it.handle())[k] - ref_p[k]);
			}
		} 
		Start_id += m_vnum*3;
	}

	if (UsePointLightDiff) {
		//gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite; // UseGDN
		//for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		//	const int vertex_id = v_it.handle().idx();									double area_scale = UseFaceArea?sqrt(m_vertices[vertex_id].BCDArea):1.0;
		//	for (int k = 0; k < 3; ++ k) {
		//		for (it = vect_const_begin(gmm::mat_const_row(gradJ, vertex_id*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, vertex_id*3+k)); ++ it) {
		//			mat_J(vertex_id*3+k+Start_id, it.index()) += area_scale*beta_scale*gradJ(vertex_id*3+k, it.index());
		//		}
		//	}
		//	mat_f(vertex_id*3+0+Start_id, 0) += area_scale*beta_scale*(T_Mesh.normal(v_it).data()[0] - m_vertices[vertex_id].ps_normal_x);
		//	mat_f(vertex_id*3+1+Start_id, 0) += area_scale*beta_scale*(T_Mesh.normal(v_it).data()[1] - m_vertices[vertex_id].ps_normal_y);
		//	mat_f(vertex_id*3+2+Start_id, 0) += area_scale*beta_scale*(T_Mesh.normal(v_it).data()[2] - m_vertices[vertex_id].ps_normal_z);
		//}
//		double max_omega = 1, min_omega = 0;
//		max_omega = -1; min_omega = 10000;
//		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
//			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
//			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
//			double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
//#ifdef TEST_MESHREFINE
//			OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
//			OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
//			omegaij = std::exp(-(psni-psnj).norm());
//#endif
//			if (omegaij < min_omega) {
//				min_omega = omegaij;
//			}
//			if (omegaij > max_omega) {
//				max_omega = omegaij;
//			}
//		}
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		//the light of neighbor points should be the same
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx(); 
			double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
			//omegaij = (omegaij-min_omega)/(max_omega-min_omega);
#ifdef TEST_MESHREFINE
			OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
			OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
			omegaij = (std::exp(-(psni-psnj).norm()) - min_omega)/(max_omega - min_omega);
#endif
			OpenMesh::Vec3f i_light = OpenMesh::Vec3f(m_vertices[iid].light_x,m_vertices[iid].light_y,m_vertices[iid].light_z);
			OpenMesh::Vec3f j_light = OpenMesh::Vec3f(m_vertices[jid].light_x,m_vertices[jid].light_y,m_vertices[jid].light_z);

			int edge_id = e_it.handle().idx();												double area_scale = UseFaceArea?sqrt(m_vertices[iid].BCDArea):1.0;
			for (int k = 0; k < 3; ++ k) {
				mat_J(edge_id*3+k+Start_id, iid+m_vnum*k+m_vnum*3) += area_scale*pld_scale*sqrt(omegaij)*1.0;
				mat_J(edge_id*3+k+Start_id, jid+m_vnum*k+m_vnum*3) -= area_scale*pld_scale*sqrt(omegaij)*1.0;
				mat_f(edge_id*3+k+Start_id, 0) += area_scale*pld_scale*sqrt(omegaij)*(i_light[k] - j_light[k]);
			}
		}
		Start_id += T_Mesh.n_edges()*3;
	}	

	if (UsePointColor) {
		// the point color constraint, color should equal to light*normal
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		T_Mesh.update_face_normals();T_Mesh.update_vertex_normals();
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			const int vertex_id = v_it.handle().idx();										double area_scale = UseFaceArea?sqrt(m_vertices[vertex_id].BCDArea):1.0;
			OpenMesh::Vec3f cur_light = OpenMesh::Vec3f(m_vertices[vertex_id].light_x,m_vertices[vertex_id].light_y,m_vertices[vertex_id].light_z);
			for (int k = 0; k < 3; ++ k) {
				for (it = vect_const_begin(gmm::mat_const_row(gradJ, vertex_id*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, vertex_id*3+k)); ++ it) {
					mat_J(vertex_id+Start_id, it.index()) += area_scale*pc_scale*cur_light[k]*gradJ(vertex_id*3+k, it.index());
				}
				mat_J(vertex_id+Start_id, vertex_id+m_vnum*k+m_vnum*3) += area_scale*pc_scale*T_Mesh.normal(v_it)[k];
			}
			mat_f(vertex_id+Start_id, 0) += area_scale*pc_scale*(OpenMesh::dot(cur_light, T_Mesh.normal(v_it))-m_vertices[vertex_id].intensity);
		}
		Start_id += m_vnum;
	}
	
	if (UsePointColorDiff) {
		// the normal difference term, eta*\omega_ij\|n_i - n_j\|^2
//		double max_omega = 1, min_omega = 0;
//		max_omega = -1; min_omega = 10000;
//		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
//			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
//			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
//			double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
//			omegaij = 1.0/(pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)+epsilon);
//#ifdef TEST_MESHREFINE
//			OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
//			OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
//			omegaij = std::exp(-(psni-psnj).norm());
//#endif
//			if (omegaij < min_omega) {
//				min_omega = omegaij;
//			}
//			if (omegaij > max_omega) {
//				max_omega = omegaij;
//			}
//		}
//		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
//			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();			double area_scale = UseFaceArea?sqrt(m_vertices[iid].BCDArea):1.0;
//			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
//
//			double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
//			omegaij = 1.0/(pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)+epsilon);
//			omegaij = (omegaij-min_omega)/(max_omega-min_omega);
//
//#ifdef TEST_MESHREFINE
//			OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
//			OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
//			omegaij = (std::exp(-(psni-psnj).norm()) - min_omega)/(max_omega - min_omega);
//#endif
//			for (int k = 0; k < 3; k++) {
//				for (it = vect_const_begin(gmm::mat_const_row(gradJ, iid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, iid*3+k)); ++ it) {
//					mat_J(e_it.handle().idx()*3+k+Start_id, it.index()) += area_scale*ndp_scale*sqrt(omegaij)*gradJ(iid*3+k, it.index());
//				}
//				for (it = vect_const_begin(gmm::mat_const_row(gradJ, jid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, jid*3+k)); ++ it) {
//					mat_J(e_it.handle().idx()*3+k+Start_id, it.index()) -= area_scale*ndp_scale*sqrt(omegaij)*gradJ(jid*3+k, it.index());
//				}
//				mat_f(e_it.handle().idx()*3+k+Start_id, 0) += area_scale*ndp_scale*sqrt(omegaij)*(T_Mesh.normal(MyMesh::VertexHandle(iid))[k] - T_Mesh.normal(MyMesh::VertexHandle(jid))[k]);
//			}
//		}


		// the new model of normal enhancement according to ci-cj
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		T_Mesh.update_face_normals();T_Mesh.update_vertex_normals();
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();		double area_scale = UseFaceArea?sqrt(m_vertices[iid].BCDArea):1.0;
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();

			OpenMesh::Vec3f i_light = OpenMesh::Vec3f(m_vertices[iid].light_x,m_vertices[iid].light_y,m_vertices[iid].light_z);
			OpenMesh::Vec3f j_light = OpenMesh::Vec3f(m_vertices[jid].light_x,m_vertices[jid].light_y,m_vertices[jid].light_z);
			for (int k = 0; k < 3; k++) {
				for (it = vect_const_begin(gmm::mat_const_row(gradJ, iid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, iid*3+k)); ++ it) {
					mat_J(e_it.handle().idx()+Start_id, it.index()) += area_scale*pcd_scale*i_light[k]*gradJ(iid*3+k, it.index());
				}	mat_J(e_it.handle().idx()+Start_id, iid+m_vnum*k+m_vnum*3) += area_scale*pcd_scale*T_Mesh.normal(MyMesh::VertexHandle(iid))[k];

				for (it = vect_const_begin(gmm::mat_const_row(gradJ, jid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, jid*3+k)); ++ it) {
					mat_J(e_it.handle().idx()+Start_id, it.index()) -= area_scale*pcd_scale*j_light[k]*gradJ(jid*3+k, it.index());
				}	mat_J(e_it.handle().idx()+Start_id, jid+m_vnum*k+m_vnum*3) -= area_scale*pcd_scale*T_Mesh.normal(MyMesh::VertexHandle(jid))[k];
			}
			mat_f(e_it.handle().idx()+Start_id, 0) += area_scale*pcd_scale*
				( (OpenMesh::dot(i_light, T_Mesh.normal(MyMesh::VertexHandle(iid))) - OpenMesh::dot(j_light, T_Mesh.normal(MyMesh::VertexHandle(jid)))) - 
					(m_vertices[iid].intensity - m_vertices[jid].intensity) );
		}
		Start_id += T_Mesh.n_edges();
	}

	if (UseFaceColorDiff) {
		T_Mesh.update_face_normals();T_Mesh.update_vertex_normals();
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
			MyMesh::FaceHandle fh1 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,0)); 
			MyMesh::FaceHandle fh2 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,1));	double area_scale = UseFaceArea?sqrt(m_trianglesArea[fh1.idx()]):1.0; 
			OpenMesh::Vec3f i_light = OpenMesh::Vec3f(m_vertices[iid].light_x,m_vertices[iid].light_y,m_vertices[iid].light_z);
			OpenMesh::Vec3f j_light = OpenMesh::Vec3f(m_vertices[jid].light_x,m_vertices[jid].light_y,m_vertices[jid].light_z);

			MyMesh::ConstFaceVertexIter cfv_it = T_Mesh.cfv_iter(fh1);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			} int k1id = cfv_it.handle().idx();
			OpenMesh::Vec3f k1_light = OpenMesh::Vec3f(m_vertices[k1id].light_x,m_vertices[k1id].light_y,m_vertices[k1id].light_z);
			double Cf1 = (m_vertices[iid].intensity + m_vertices[jid].intensity + m_vertices[k1id].intensity)/3.0;
			OpenMesh::Vec3f avg_light1 = (i_light+j_light+k1_light)/3.0;

			cfv_it = T_Mesh.cfv_iter(fh2);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			} int k2id = cfv_it.handle().idx();
			OpenMesh::Vec3f k2_light = OpenMesh::Vec3f(m_vertices[k2id].light_x,m_vertices[k2id].light_y,m_vertices[k2id].light_z);
			double Cf2 = (m_vertices[iid].intensity + m_vertices[jid].intensity + m_vertices[k2id].intensity)/3.0;
			OpenMesh::Vec3f avg_light2 = (i_light+j_light+k2_light)/3.0;

			cfv_it = T_Mesh.cfv_iter(fh1);
			OpenMesh::Vec3f pa = T_Mesh.point(cfv_it.handle()); int ida = cfv_it.handle().idx(); ++ cfv_it;
			OpenMesh::Vec3f pb = T_Mesh.point(cfv_it.handle()); int idb = cfv_it.handle().idx(); ++ cfv_it;
			OpenMesh::Vec3f pc = T_Mesh.point(cfv_it.handle()); int idc = cfv_it.handle().idx(); 
			double pax = pa[0]; double pay = pa[1]; double paz = pa[2];
			double pbx = pb[0]; double pby = pb[1]; double pbz = pb[2];
			double pcx = pc[0]; double pcy = pc[1]; double pcz = pc[2];
			OpenMesh::Vec3f n2bx, n2by, n2bz, n2ax, n2ay, n2az, n2cx, n2cy, n2cz;
			
			OpenMesh::Vec3f Crx(1,0,0), Cry(0,1,0), Crz(0,0,1), Nijk, Eik, Eki;
			Nijk = OpenMesh::cross(pb-pa, pc-pa)/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon); Eik = pc-pa; Eki = pa-pc;
			n2bx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2by = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2bz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);

			Nijk = OpenMesh::cross(pa-pc, pb-pc)/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon); Eik = pb-pc; Eki = pc-pb;
			n2ax = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2ay = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2az = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);

			Nijk = OpenMesh::cross(pc-pb, pa-pb)/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon); Eik = pa-pb; Eki = pb-pa;
			n2cx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cy = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);

			for (int k = 0; k < 3; ++ k) {
				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*0) += area_scale*fcd_scale*avg_light1[k]*n2ax[k];
				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*1) += area_scale*fcd_scale*avg_light1[k]*n2ay[k];
				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*2) += area_scale*fcd_scale*avg_light1[k]*n2az[k];

				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*0) += area_scale*fcd_scale*avg_light1[k]*n2bx[k];
				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*1) += area_scale*fcd_scale*avg_light1[k]*n2by[k];
				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*2) += area_scale*fcd_scale*avg_light1[k]*n2bz[k];

				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*0) += area_scale*fcd_scale*avg_light1[k]*n2cx[k];
				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*1) += area_scale*fcd_scale*avg_light1[k]*n2cy[k];
				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*2) += area_scale*fcd_scale*avg_light1[k]*n2cz[k];

				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*k+m_vnum*3) += area_scale*fcd_scale*Nijk[k]/3.0;
				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*k+m_vnum*3) += area_scale*fcd_scale*Nijk[k]/3.0;
				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*k+m_vnum*3) += area_scale*fcd_scale*Nijk[k]/3.0;
			}

			cfv_it = T_Mesh.cfv_iter(fh2);
			pa = T_Mesh.point(cfv_it.handle()); ida = cfv_it.handle().idx(); ++ cfv_it;
			pb = T_Mesh.point(cfv_it.handle()); idb = cfv_it.handle().idx(); ++ cfv_it;
			pc = T_Mesh.point(cfv_it.handle()); idc = cfv_it.handle().idx(); 
			pax = pa[0]; pay = pa[1]; paz = pa[2];
			pbx = pb[0]; pby = pb[1]; pbz = pb[2];
			pcx = pc[0]; pcy = pc[1]; pcz = pc[2];		

			Nijk = OpenMesh::cross(pb-pa, pc-pa)/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon); Eik = pc-pa; Eki = pa-pc;
			n2bx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2by = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2bz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);

			Nijk = OpenMesh::cross(pa-pc, pb-pc)/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon); Eik = pb-pc; Eki = pc-pb;
			n2ax = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2ay = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2az = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);

			Nijk = OpenMesh::cross(pc-pb, pa-pb)/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon); Eik = pa-pb; Eki = pb-pa;
			n2cx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cy = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);

			for (int k = 0; k < 3; ++ k) {
				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*0) -= area_scale*fcd_scale*avg_light2[k]*n2ax[k];
				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*1) -= area_scale*fcd_scale*avg_light2[k]*n2ay[k];
				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*2) -= area_scale*fcd_scale*avg_light2[k]*n2az[k];

				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*0) -= area_scale*fcd_scale*avg_light2[k]*n2bx[k];
				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*1) -= area_scale*fcd_scale*avg_light2[k]*n2by[k];
				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*2) -= area_scale*fcd_scale*avg_light2[k]*n2bz[k];

				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*0) -= area_scale*fcd_scale*avg_light2[k]*n2cx[k];
				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*1) -= area_scale*fcd_scale*avg_light2[k]*n2cy[k];
				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*2) -= area_scale*fcd_scale*avg_light2[k]*n2cz[k];

				mat_J(e_it.handle().idx()+Start_id, ida+m_vnum*k+m_vnum*3) -= area_scale*fcd_scale*Nijk[k]/3.0;
				mat_J(e_it.handle().idx()+Start_id, idb+m_vnum*k+m_vnum*3) -= area_scale*fcd_scale*Nijk[k]/3.0;
				mat_J(e_it.handle().idx()+Start_id, idc+m_vnum*k+m_vnum*3) -= area_scale*fcd_scale*Nijk[k]/3.0;
			}
			mat_f(e_it.handle().idx()+Start_id, 0) += area_scale*fcd_scale*
				( OpenMesh::dot(avg_light1, T_Mesh.normal(fh1)) - OpenMesh::dot(avg_light2, T_Mesh.normal(fh2)) - (Cf1 - Cf2) );
		}
		Start_id += T_Mesh.n_edges();
	} 

	if (UsePointNormalDiff) {
		// the normal difference term, eta*\omega_ij\|n_i - n_j\|^2
		double max_omega = 1, min_omega = 0;
		max_omega = -1; min_omega = 10000;
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
			double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
			//omegaij = 1.0/(pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)+epsilon);
#ifdef TEST_MESHREFINE
			OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
			OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
			omegaij = std::exp(-(psni-psnj).norm());
#endif
			if (omegaij < min_omega) {
				min_omega = omegaij;
			}
			if (omegaij > max_omega) {
				max_omega = omegaij;
			}
		}
		T_Mesh.update_face_normals();T_Mesh.update_vertex_normals();
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();			double area_scale = UseFaceArea?sqrt(m_vertices[iid].BCDArea):1.0;
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();

			double omegaij = std::exp(-pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);
			//omegaij = 1.0/(pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)+epsilon);
			omegaij = (omegaij-min_omega)/(max_omega-min_omega);

#ifdef TEST_MESHREFINE
			OpenMesh::Vec3f psni(m_vertices[iid].ps_normal_x,m_vertices[iid].ps_normal_y,m_vertices[iid].ps_normal_z);
			OpenMesh::Vec3f psnj(m_vertices[jid].ps_normal_x,m_vertices[jid].ps_normal_y,m_vertices[jid].ps_normal_z);
			omegaij = (std::exp(-(psni-psnj).norm()) - min_omega)/(max_omega - min_omega);
#endif
			for (int k = 0; k < 3; k++) {
				for (it = vect_const_begin(gmm::mat_const_row(gradJ, iid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, iid*3+k)); ++ it) {
					mat_J(e_it.handle().idx()*3+k+Start_id, it.index()) += area_scale*pnd_scale*sqrt(omegaij)*gradJ(iid*3+k, it.index());
				}
				for (it = vect_const_begin(gmm::mat_const_row(gradJ, jid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, jid*3+k)); ++ it) {
					mat_J(e_it.handle().idx()*3+k+Start_id, it.index()) -= area_scale*pnd_scale*sqrt(omegaij)*gradJ(jid*3+k, it.index());
				}
				mat_f(e_it.handle().idx()*3+k+Start_id, 0) += area_scale*pnd_scale*sqrt(omegaij)*(T_Mesh.normal(MyMesh::VertexHandle(iid))[k] - T_Mesh.normal(MyMesh::VertexHandle(jid))[k]);
			}
		}
		Start_id += T_Mesh.n_edges()*3;
	}

	if (UseFaceNormalDiff) {
		double max_omega = 1, min_omega = 0;
		max_omega = -1; min_omega = 10000;
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			MyMesh::FaceHandle fh1 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,0)); 
			MyMesh::FaceHandle fh2 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,1));
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();

			MyMesh::ConstFaceVertexIter cfv_it = T_Mesh.cfv_iter(fh1);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			}  int k1id = cfv_it.handle().idx();

			cfv_it = T_Mesh.cfv_iter(fh2);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			}  int k2id = cfv_it.handle().idx();

			double omegaij = std::exp(-pow((m_vertices[k1id].intensity - m_vertices[k2id].intensity), 2.0)/varsigma);
#ifdef TEST_MESHREFINE
			OpenMesh::Vec3f psni(m_vertices[k1id].ps_normal_x,m_vertices[k1id].ps_normal_y,m_vertices[k1id].ps_normal_z);
			OpenMesh::Vec3f psnj(m_vertices[k2id].ps_normal_x,m_vertices[k2id].ps_normal_y,m_vertices[k2id].ps_normal_z);
			omegaij = std::exp(-(psni-psnj).norm());
#endif
			if (omegaij < min_omega) {
				min_omega = omegaij;
			}
			if (omegaij > max_omega) {
				max_omega = omegaij;
			}
		}
		T_Mesh.update_face_normals();T_Mesh.update_vertex_normals();
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,0)).idx();
			int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(e_it,1)).idx();
			MyMesh::FaceHandle fh1 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,0)); 
			MyMesh::FaceHandle fh2 = T_Mesh.face_handle(T_Mesh.halfedge_handle(e_it,1));	double area_scale = UseFaceArea?sqrt(m_trianglesArea[fh1.idx()]):1.0; \

			MyMesh::ConstFaceVertexIter cfv_it = T_Mesh.cfv_iter(fh1);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			} int k1id = cfv_it.handle().idx();

			cfv_it = T_Mesh.cfv_iter(fh2);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			} int k2id = cfv_it.handle().idx();

			double omegaij = std::exp(-pow((m_vertices[k1id].intensity - m_vertices[k2id].intensity), 2.0)/varsigma);
			omegaij = (omegaij-min_omega)/(max_omega-min_omega);
#ifdef TEST_MESHREFINE
			OpenMesh::Vec3f psni(m_vertices[k1id].ps_normal_x,m_vertices[k1id].ps_normal_y,m_vertices[k1id].ps_normal_z);
			OpenMesh::Vec3f psnj(m_vertices[k2id].ps_normal_x,m_vertices[k2id].ps_normal_y,m_vertices[k2id].ps_normal_z);
			omegaij = (std::exp(-(psni-psnj).norm()) - min_omega)/(max_omega - min_omega);
#endif
			cfv_it = T_Mesh.cfv_iter(fh1);
			OpenMesh::Vec3f pa = T_Mesh.point(cfv_it.handle()); int ida = cfv_it.handle().idx(); ++ cfv_it;
			OpenMesh::Vec3f pb = T_Mesh.point(cfv_it.handle()); int idb = cfv_it.handle().idx(); ++ cfv_it;
			OpenMesh::Vec3f pc = T_Mesh.point(cfv_it.handle()); int idc = cfv_it.handle().idx(); 
			double pax = pa[0]; double pay = pa[1]; double paz = pa[2];
			double pbx = pb[0]; double pby = pb[1]; double pbz = pb[2];
			double pcx = pc[0]; double pcy = pc[1]; double pcz = pc[2];
			OpenMesh::Vec3f n2bx, n2by, n2bz, n2ax, n2ay, n2az, n2cx, n2cy, n2cz;

			OpenMesh::Vec3f Crx(1,0,0), Cry(0,1,0), Crz(0,0,1), Nijk, Eik, Eki;
			Nijk = OpenMesh::cross(pb-pa, pc-pa)/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon); Eik = pc-pa; Eki = pa-pc;
			n2bx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2by = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2bz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);

			Nijk = OpenMesh::cross(pa-pc, pb-pc)/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon); Eik = pb-pc; Eki = pc-pb;
			n2ax = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2ay = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2az = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);

			Nijk = OpenMesh::cross(pc-pb, pa-pb)/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon); Eik = pa-pb; Eki = pb-pa;
			n2cx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cy = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);

			for (int k = 0; k < 3; k++) {
				mat_J(e_it.handle().idx()*3+k+Start_id, ida+m_vnum*0) += area_scale*fnd_scale*sqrt(omegaij)*n2ax[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, ida+m_vnum*1) += area_scale*fnd_scale*sqrt(omegaij)*n2ay[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, ida+m_vnum*2) += area_scale*fnd_scale*sqrt(omegaij)*n2az[k];

				mat_J(e_it.handle().idx()*3+k+Start_id, idb+m_vnum*0) += area_scale*fnd_scale*sqrt(omegaij)*n2bx[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idb+m_vnum*1) += area_scale*fnd_scale*sqrt(omegaij)*n2by[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idb+m_vnum*2) += area_scale*fnd_scale*sqrt(omegaij)*n2bz[k];

				mat_J(e_it.handle().idx()*3+k+Start_id, idc+m_vnum*0) += area_scale*fnd_scale*sqrt(omegaij)*n2cx[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idc+m_vnum*1) += area_scale*fnd_scale*sqrt(omegaij)*n2cy[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idc+m_vnum*2) += area_scale*fnd_scale*sqrt(omegaij)*n2cz[k];
				//mat_f[e_it.handle().idx()*3+k+Start_id, 0] += area_scale*ndf_scale*sqrt(omegaij)*(Nijk.data()[k]);
			}

			cfv_it = T_Mesh.cfv_iter(fh2);
			pa = T_Mesh.point(cfv_it.handle()); ida = cfv_it.handle().idx(); ++ cfv_it;
			pb = T_Mesh.point(cfv_it.handle()); idb = cfv_it.handle().idx(); ++ cfv_it;
			pc = T_Mesh.point(cfv_it.handle()); idc = cfv_it.handle().idx(); 
			pax = pa[0]; pay = pa[1]; paz = pa[2];
			pbx = pb[0]; pby = pb[1]; pbz = pb[2];
			pcx = pc[0]; pcy = pc[1]; pcz = pc[2];		

			Nijk = OpenMesh::cross(pb-pa, pc-pa)/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon); Eik = pc-pa; Eki = pa-pc;
			n2bx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2by = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);
			n2bz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pb-pa, pc-pa).norm()+epsilon);

			Nijk = OpenMesh::cross(pa-pc, pb-pc)/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon); Eik = pb-pc; Eki = pc-pb;
			n2ax = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2ay = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);
			n2az = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pa-pc, pb-pc).norm()+epsilon);

			Nijk = OpenMesh::cross(pc-pb, pa-pb)/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon); Eik = pa-pb; Eki = pb-pa;
			n2cx = ( OpenMesh::cross(Crx, Eik) - Nijk*OpenMesh::dot(Crx, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cy = ( OpenMesh::cross(Cry, Eik) - Nijk*OpenMesh::dot(Cry, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);
			n2cz = ( OpenMesh::cross(Crz, Eik) - Nijk*OpenMesh::dot(Crz, OpenMesh::cross(Nijk, Eki)) )/(OpenMesh::cross(pc-pb, pa-pb).norm()+epsilon);

			for (int k = 0; k < 3; k++) {
				mat_J(e_it.handle().idx()*3+k+Start_id, ida+m_vnum*0) -= area_scale*fnd_scale*sqrt(omegaij)*n2ax[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, ida+m_vnum*1) -= area_scale*fnd_scale*sqrt(omegaij)*n2ay[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, ida+m_vnum*2) -= area_scale*fnd_scale*sqrt(omegaij)*n2az[k];

				mat_J(e_it.handle().idx()*3+k+Start_id, idb+m_vnum*0) -= area_scale*fnd_scale*sqrt(omegaij)*n2bx[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idb+m_vnum*1) -= area_scale*fnd_scale*sqrt(omegaij)*n2by[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idb+m_vnum*2) -= area_scale*fnd_scale*sqrt(omegaij)*n2bz[k];

				mat_J(e_it.handle().idx()*3+k+Start_id, idc+m_vnum*0) -= area_scale*fnd_scale*sqrt(omegaij)*n2cx[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idc+m_vnum*1) -= area_scale*fnd_scale*sqrt(omegaij)*n2cy[k];
				mat_J(e_it.handle().idx()*3+k+Start_id, idc+m_vnum*2) -= area_scale*fnd_scale*sqrt(omegaij)*n2cz[k];
				//mat_f[e_it.handle().idx()*3+k+Start_id, 0] -= area_scale*ndf_scale*sqrt(omegaij)*(Nijk.data()[k]);
			}

			for (int k = 0; k < 3; ++ k) {
				mat_f(e_it.handle().idx()*3+k+Start_id, 0) += area_scale*fnd_scale*sqrt(omegaij)*(T_Mesh.normal(fh1).data()[k] - T_Mesh.normal(fh2).data()[k]);
			}
		}
		Start_id += T_Mesh.n_edges()*3;
	} 

	if (UseLaplace) {
		if (AnisotropicLaplace) {
			double CValue = 50; // according to cvpr paper, to set the weight for anisotropic laplacian term
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				int vertex_id = v_it.handle().idx();									double area_scale = UseFaceArea?sqrt(m_vertices[vertex_id].BCDArea):1.0;
				//add the laplacian term for singularity mesh
				double degree = 0.0, avg_color = 0.0, sum_weight = 0.0; 
				for (MyMesh::ConstVertexVertexIter vv_it = T_Mesh.cvv_iter(v_it); vv_it; ++vv_it) {
					degree++;	avg_color += m_vertices[vv_it.handle().idx()].intensity; 
				}  avg_color = avg_color/degree;
				double out_weight = 1.0;//sqrt(exp(-pow(m_vertices[v_it.handle().idx()].intensity - avg_color,2.0)/varsigma))+epsilon;

				OpenMesh::Vec3f avg_nei_point; avg_nei_point.vectorize(0.0); vector<double> nei_weight; nei_weight.clear();
				for (MyMesh::VertexEdgeIter ve_it = T_Mesh.ve_iter(v_it); ve_it; ++ ve_it) {
					MyMesh::FaceHandle fh1 = T_Mesh.face_handle(T_Mesh.halfedge_handle(ve_it,0)); 
					MyMesh::FaceHandle fh2 = T_Mesh.face_handle(T_Mesh.halfedge_handle(ve_it,1));
					int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(ve_it,0)).idx();
					int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(ve_it,1)).idx();
					MyMesh::ConstFaceVertexIter cfv_it = T_Mesh.cfv_iter(fh1);
					while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
						++cfv_it;
					}  int k1id = cfv_it.handle().idx();
					cfv_it = T_Mesh.cfv_iter(fh2);
					while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
						++cfv_it;
					}  int k2id = cfv_it.handle().idx();
					if (iid != vertex_id) { // jid = vertex_id
						jid = iid; iid = vertex_id;
					}
					OpenMesh::Vec3f pPi  = T_Mesh.point(MyMesh::VertexHandle(iid)),		pPj = T_Mesh.point(MyMesh::VertexHandle(jid));
					OpenMesh::Vec3f pPk1 = T_Mesh.point(MyMesh::VertexHandle(k1id)),	pPk2 = T_Mesh.point(MyMesh::VertexHandle(k2id));
					double cotaij = OpenMesh::dot(pPi-pPk1, pPj-pPk1)/OpenMesh::cross(pPi-pPk1, pPj-pPk1).norm();
					double cotbij = OpenMesh::dot(pPi-pPk2, pPj-pPk2)/OpenMesh::cross(pPi-pPk2, pPj-pPk2).norm();

					double col_edge_weight = 1 - min(abs(m_vertices[iid].intensity - m_vertices[jid].intensity), CValue)/CValue;
					double lap_edge_weight = (cotaij+cotbij)/(2.0*m_vertices[vertex_id].Voronoi_Area);

					avg_nei_point += pPj*(col_edge_weight*lap_edge_weight);
					nei_weight.push_back(col_edge_weight*lap_edge_weight);	sum_weight += col_edge_weight*lap_edge_weight;
				} 

				for (int k = 0; k < 3; ++ k) {
					int i = 0;		mat_J(vertex_id*3+k+Start_id, vertex_id+m_vnum*k) += area_scale*lap_scale*out_weight*sum_weight; 
					for (MyMesh::VertexEdgeIter ve_it = T_Mesh.ve_iter(v_it); ve_it; ++ ve_it, ++ i) {
						int iid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(ve_it,0)).idx();
						int jid = T_Mesh.to_vertex_handle(T_Mesh.halfedge_handle(ve_it,1)).idx();
						if (iid != vertex_id) { // jid = vertex_id
							jid = iid; iid = vertex_id;
						}
						mat_J(vertex_id*3+k+Start_id, jid+m_vnum*k) -= area_scale*lap_scale*out_weight*nei_weight[i];
					}
					mat_f(vertex_id*3+k+Start_id, 0) += area_scale*lap_scale*out_weight*(T_Mesh.point(v_it.handle()).data()[k]*sum_weight - avg_nei_point[k]);
				}
			}
		} else {
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				int vertex_id = v_it.handle().idx();										double area_scale = UseFaceArea?sqrt(m_vertices[vertex_id].BCDArea):1.0;
				//add the laplacian term for singularity mesh
				double degree = 0.0, avg_color = 0.0, sum_weight = 0.0; OpenMesh::Vec3f avg_nei_point; avg_nei_point.vectorize(0.0); vector<double> nei_weight; nei_weight.clear();
				for (MyMesh::ConstVertexVertexIter vv_it = T_Mesh.cvv_iter(v_it); vv_it; ++vv_it) {
					degree++;	avg_color += m_vertices[vv_it.handle().idx()].intensity; 

					double col_edge_weight = exp(-pow(m_vertices[v_it.handle().idx()].intensity - m_vertices[vv_it.handle().idx()].intensity,2.0)/varsigma);
					avg_nei_point += T_Mesh.point(vv_it.handle())*col_edge_weight;
					nei_weight.push_back(col_edge_weight);	sum_weight += col_edge_weight;
				} avg_color = avg_color/degree;

				double out_weight = exp(-pow(m_vertices[v_it.handle().idx()].intensity - avg_color,2.0)/varsigma);
				for (int k = 0; k < 3; ++ k) {
					int i = 0;		mat_J(vertex_id*3+k+Start_id, vertex_id+m_vnum*k) += area_scale*lap_scale*out_weight*sum_weight; 
					for (MyMesh::ConstVertexVertexIter vv_it = T_Mesh.cvv_iter(v_it); vv_it; ++vv_it, ++ i) {
						mat_J(vertex_id*3+k+Start_id, vv_it.handle().idx()+m_vnum*k) -= area_scale*lap_scale*out_weight*nei_weight[i];
					}
					mat_f(vertex_id*3+k+Start_id, 0) += area_scale*lap_scale*out_weight*(T_Mesh.point(v_it.handle()).data()[k]*sum_weight - avg_nei_point[k]);
				}
			}
		}
		Start_id += m_vnum*3;
	}

	if (UseTVNorm) {
		// the tvn term, \gamma/2 \|(p+\lambda^k/\gamma)-\nabla n(u)\|^2, this is defined on each triangle
		vector<double>nux, nuy, nuz; nux.resize(m_vnum); nuy.resize(m_vnum); nuz.resize(m_vnum);
		T_Mesh.update_face_normals();  T_Mesh.update_vertex_normals();
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
			nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
			nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
			nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
		}
		GetPPIGradients (m_vnum,nux);
		GetPPIGradients1(m_vnum,nuy);
		GetPPIGradients2(m_vnum,nuz);
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		for (int i = 0; i < m_trinum; ++ i) {
			double area_scale = UseFaceArea?sqrt(m_trianglesArea[i]):1.0;
			// need to calculate the derivative of m_TPPIBG[i] to each vertex;
			gmm::dense_matrix<double> h0_ver0(3,3), h0_ver1(3,3), h0_ver2(3,3), h1_ver0(3,3), h1_ver1(3,3), h1_ver2(3,3), h2_ver0(3,3), h2_ver1(3,3), h2_ver2(3,3), Mat_I(3,3);
			gmm::copy(gmm::identity_matrix(), Mat_I); 

			gmm::dense_matrix<double> pPi(3,1), pPj(3,1), pPk(3,1), Eij(3,1), eij(3,1), Eik(3,1), eijeijT(3,3), m_eijeijT(3,3), temp133(3,3), temp233(3,3), temp333(3,3), temp11(1,1);
			gmm::dense_matrix<double> Hk(3,1), hk(3,1), NHk(3,3);
			// calculate the derivative of h2
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk); //cout << NHk;

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h2_ver2); // calculate h22ver2
			gmm::mult(NHk, h2_ver2, temp333); gmm::copy(temp333, h2_ver2); //cout << h2_ver2;

			gmm::copy(Mat_I, h2_ver0); gmm::scale(h2_ver0, -1.0); gmm::add(eijeijT, h2_ver0); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h2_ver0, temp133, temp333);	gmm::copy(temp333, h2_ver0); //cout << h1_ver2;// calculate h12ver2
			gmm::mult(NHk, h2_ver0, temp333);		gmm::copy(temp333, h2_ver0); //cout << h2_ver0;

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij2vj
			gmm::mult(temp333, temp233, h2_ver1); //cout << h1_ver0;
			gmm::mult(NHk, h2_ver1, temp333);		gmm::copy(temp333, h2_ver1); //cout << h2_ver1;

			// calculate the derivative of h1
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk); //cout << NHk;

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h1_ver1); //cout<< h1_ver1;// calculate h12ver1
			gmm::mult(NHk, h1_ver1, temp333); gmm::copy(temp333, h1_ver1); //cout << h1_ver1;

			gmm::copy(Mat_I, h1_ver2); gmm::scale(h1_ver2, -1.0); gmm::add(eijeijT, h1_ver2); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h1_ver2, temp133, temp333);	gmm::copy(temp333, h1_ver2); //cout << h1_ver2;// calculate h12ver2
			gmm::mult(NHk, h1_ver2, temp333);		gmm::copy(temp333, h1_ver2); //cout << h1_ver2;

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij2vj
			gmm::mult(temp333, temp233, h1_ver0); //cout << h1_ver0;
			gmm::mult(NHk, h1_ver0, temp333);		gmm::copy(temp333, h1_ver0); //cout << h1_ver0;

			// calculate the derivative of h0
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk); //cout << NHk;

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h0_ver0); //cout << h0_ver0;// calculate h12ver1
			gmm::mult(NHk, h0_ver0, temp333); gmm::copy(temp333, h0_ver0); //cout << h0_ver0;

			gmm::copy(Mat_I, h0_ver1); gmm::scale(h0_ver1, -1.0); gmm::add(eijeijT, h0_ver1); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h0_ver1, temp133, temp333);	gmm::copy(temp333, h0_ver1); //cout << h0_ver1;// calculate h12ver2
			gmm::mult(NHk, h0_ver1, temp333);		gmm::copy(temp333, h0_ver1); //cout << h0_ver1;

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij2vj
			gmm::mult(temp333, temp233, h0_ver2); //cout << h0_ver2;
			gmm::mult(NHk, h0_ver2, temp333);		gmm::copy(temp333, h0_ver2);  //cout << h0_ver2;

			// \partial grad.x/\partial nx \times \partial nx/\partial u
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(0,0) + nux[m_triangles[i].ver1]*h1_ver0(0,0) + nux[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(0,1) + nux[m_triangles[i].ver1]*h1_ver0(0,1) + nux[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(0,2) + nux[m_triangles[i].ver1]*h1_ver0(0,2) + nux[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(0,0) + nux[m_triangles[i].ver1]*h1_ver1(0,0) + nux[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(0,1) + nux[m_triangles[i].ver1]*h1_ver1(0,1) + nux[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(0,2) + nux[m_triangles[i].ver1]*h1_ver1(0,2) + nux[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(0,0) + nux[m_triangles[i].ver1]*h1_ver2(0,0) + nux[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(0,1) + nux[m_triangles[i].ver1]*h1_ver2(0,1) + nux[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(0,2) + nux[m_triangles[i].ver1]*h1_ver2(0,2) + nux[m_triangles[i].ver2]*h2_ver2(0,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); ++ it) {
				mat_J(i*9+0+Start_id , it.index()) += m_TPPIBG[i].v0.x*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+0, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); ++ it) {
				mat_J(i*9+0+Start_id , it.index()) += m_TPPIBG[i].v1.x*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+0, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); ++ it) {
				mat_J(i*9+0+Start_id , it.index()) += m_TPPIBG[i].v2.x*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+0, it.index());
			}
			mat_f(i*9+0+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad.x - (lambda_x[i].x/penParam + px[i].x));

			// \partial grad.y/\partial nx \times \partial nx/\partial u
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(1,0) + nux[m_triangles[i].ver1]*h1_ver0(1,0) + nux[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(1,1) + nux[m_triangles[i].ver1]*h1_ver0(1,1) + nux[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(1,2) + nux[m_triangles[i].ver1]*h1_ver0(1,2) + nux[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(1,0) + nux[m_triangles[i].ver1]*h1_ver1(1,0) + nux[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(1,1) + nux[m_triangles[i].ver1]*h1_ver1(1,1) + nux[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(1,2) + nux[m_triangles[i].ver1]*h1_ver1(1,2) + nux[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(1,0) + nux[m_triangles[i].ver1]*h1_ver2(1,0) + nux[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(1,1) + nux[m_triangles[i].ver1]*h1_ver2(1,1) + nux[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(1,2) + nux[m_triangles[i].ver1]*h1_ver2(1,2) + nux[m_triangles[i].ver2]*h2_ver2(1,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); ++ it) {
				mat_J(i*9+1+Start_id , it.index()) += m_TPPIBG[i].v0.y*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+0, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); ++ it) {
				mat_J(i*9+1+Start_id , it.index()) += m_TPPIBG[i].v1.y*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+0, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); ++ it) {
				mat_J(i*9+1+Start_id , it.index()) += m_TPPIBG[i].v2.y*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+0, it.index());
			}
			mat_f(i*9+1+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad.y - (lambda_x[i].y/penParam + px[i].y));

			// \partial grad.z/\partial nx \times \partial nx/\partial u
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(2,0) + nux[m_triangles[i].ver1]*h1_ver0(2,0) + nux[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(2,1) + nux[m_triangles[i].ver1]*h1_ver0(2,1) + nux[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(2,2) + nux[m_triangles[i].ver1]*h1_ver0(2,2) + nux[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(2,0) + nux[m_triangles[i].ver1]*h1_ver1(2,0) + nux[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(2,1) + nux[m_triangles[i].ver1]*h1_ver1(2,1) + nux[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(2,2) + nux[m_triangles[i].ver1]*h1_ver1(2,2) + nux[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(2,0) + nux[m_triangles[i].ver1]*h1_ver2(2,0) + nux[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(2,1) + nux[m_triangles[i].ver1]*h1_ver2(2,1) + nux[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(2,2) + nux[m_triangles[i].ver1]*h1_ver2(2,2) + nux[m_triangles[i].ver2]*h2_ver2(2,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); ++ it) {
				mat_J(i*9+2+Start_id , it.index()) += m_TPPIBG[i].v0.z*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+0, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); ++ it) {
				mat_J(i*9+2+Start_id , it.index()) += m_TPPIBG[i].v1.z*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+0, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); ++ it) {
				mat_J(i*9+2+Start_id , it.index()) += m_TPPIBG[i].v2.z*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+0, it.index());
			}
			mat_f(i*9+2+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad.z - (lambda_x[i].z/penParam + px[i].z));
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// \partial grad1.x/\partial ny \times \partial ny/\partial u
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(0,0) + nuy[m_triangles[i].ver1]*h1_ver0(0,0) + nuy[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(0,1) + nuy[m_triangles[i].ver1]*h1_ver0(0,1) + nuy[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(0,2) + nuy[m_triangles[i].ver1]*h1_ver0(0,2) + nuy[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(0,0) + nuy[m_triangles[i].ver1]*h1_ver1(0,0) + nuy[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*1) +=  
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(0,1) + nuy[m_triangles[i].ver1]*h1_ver1(0,1) + nuy[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(0,2) + nuy[m_triangles[i].ver1]*h1_ver1(0,2) + nuy[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(0,0) + nuy[m_triangles[i].ver1]*h1_ver2(0,0) + nuy[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(0,1) + nuy[m_triangles[i].ver1]*h1_ver2(0,1) + nuy[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(0,2) + nuy[m_triangles[i].ver1]*h1_ver2(0,2) + nuy[m_triangles[i].ver2]*h2_ver2(0,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); ++ it) {
				mat_J(i*9+3+Start_id , it.index()) += m_TPPIBG[i].v0.x*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+1, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); ++ it) {
				mat_J(i*9+3+Start_id , it.index()) += m_TPPIBG[i].v1.x*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+1, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); ++ it) {
				mat_J(i*9+3+Start_id , it.index()) += m_TPPIBG[i].v2.x*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+1, it.index());
			}
			mat_f(i*9+3+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad1.x - (lambda_y[i].x/penParam + py[i].x));

			// \partial grad1.y/\partial ny \times \partial ny/\partial u
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(1,0) + nuy[m_triangles[i].ver1]*h1_ver0(1,0) + nuy[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(1,1) + nuy[m_triangles[i].ver1]*h1_ver0(1,1) + nuy[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(1,2) + nuy[m_triangles[i].ver1]*h1_ver0(1,2) + nuy[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(1,0) + nuy[m_triangles[i].ver1]*h1_ver1(1,0) + nuy[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(1,1) + nuy[m_triangles[i].ver1]*h1_ver1(1,1) + nuy[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(1,2) + nuy[m_triangles[i].ver1]*h1_ver1(1,2) + nuy[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(1,0) + nuy[m_triangles[i].ver1]*h1_ver2(1,0) + nuy[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(1,1) + nuy[m_triangles[i].ver1]*h1_ver2(1,1) + nuy[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(1,2) + nuy[m_triangles[i].ver1]*h1_ver2(1,2) + nuy[m_triangles[i].ver2]*h2_ver2(1,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); ++ it) {
				mat_J(i*9+4+Start_id , it.index()) += m_TPPIBG[i].v0.y*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+1, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); ++ it) {
				mat_J(i*9+4+Start_id , it.index()) += m_TPPIBG[i].v1.y*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+1, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); ++ it) {
				mat_J(i*9+4+Start_id , it.index()) += m_TPPIBG[i].v2.y*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+1, it.index());
			}
			mat_f(i*9+4+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad1.y - (lambda_y[i].y/penParam + py[i].y));

			// \partial grad1.z/\partial ny \times \partial ny/\partial u
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(2,0) + nuy[m_triangles[i].ver1]*h1_ver0(2,0) + nuy[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(2,1) + nuy[m_triangles[i].ver1]*h1_ver0(2,1) + nuy[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(2,2) + nuy[m_triangles[i].ver1]*h1_ver0(2,2) + nuy[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(2,0) + nuy[m_triangles[i].ver1]*h1_ver1(2,0) + nuy[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(2,1) + nuy[m_triangles[i].ver1]*h1_ver1(2,1) + nuy[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(2,2) + nuy[m_triangles[i].ver1]*h1_ver1(2,2) + nuy[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(2,0) + nuy[m_triangles[i].ver1]*h1_ver2(2,0) + nuy[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(2,1) + nuy[m_triangles[i].ver1]*h1_ver2(2,1) + nuy[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(2,2) + nuy[m_triangles[i].ver1]*h1_ver2(2,2) + nuy[m_triangles[i].ver2]*h2_ver2(2,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); ++ it) {
				mat_J(i*9+5+Start_id , it.index()) += m_TPPIBG[i].v0.z*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+1, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); ++ it) {
				mat_J(i*9+5+Start_id , it.index()) += m_TPPIBG[i].v1.z*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+1, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); ++ it) {
				mat_J(i*9+5+Start_id , it.index()) += m_TPPIBG[i].v2.z*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+1, it.index());
			}
			mat_f(i*9+5+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad1.z - (lambda_y[i].z/penParam + py[i].z));
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// \partial grad2.x/\partial nz \times \partial nz/\partial u
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(0,0) + nuz[m_triangles[i].ver1]*h1_ver0(0,0) + nuz[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(0,1) + nuz[m_triangles[i].ver1]*h1_ver0(0,1) + nuz[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(0,2) + nuz[m_triangles[i].ver1]*h1_ver0(0,2) + nuz[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(0,0) + nuz[m_triangles[i].ver1]*h1_ver1(0,0) + nuz[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(0,1) + nuz[m_triangles[i].ver1]*h1_ver1(0,1) + nuz[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(0,2) + nuz[m_triangles[i].ver1]*h1_ver1(0,2) + nuz[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(0,0) + nuz[m_triangles[i].ver1]*h1_ver2(0,0) + nuz[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(0,1) + nuz[m_triangles[i].ver1]*h1_ver2(0,1) + nuz[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(0,2) + nuz[m_triangles[i].ver1]*h1_ver2(0,2) + nuz[m_triangles[i].ver2]*h2_ver2(0,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); ++ it) {
				mat_J(i*9+6+Start_id , it.index()) += m_TPPIBG[i].v0.x*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+2, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); ++ it) {
				mat_J(i*9+6+Start_id , it.index()) += m_TPPIBG[i].v1.x*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+2, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); ++ it) {
				mat_J(i*9+6+Start_id , it.index()) += m_TPPIBG[i].v2.x*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+2, it.index());
			}
			mat_f(i*9+6+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad2.x - (lambda_z[i].x/penParam + pz[i].x));

			// \partial grad2.y/\partial nz \times \partial nz/\partial u
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(1,0) + nuz[m_triangles[i].ver1]*h1_ver0(1,0) + nuz[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(1,1) + nuz[m_triangles[i].ver1]*h1_ver0(1,1) + nuz[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(1,2) + nuz[m_triangles[i].ver1]*h1_ver0(1,2) + nuz[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(1,0) + nuz[m_triangles[i].ver1]*h1_ver1(1,0) + nuz[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(1,1) + nuz[m_triangles[i].ver1]*h1_ver1(1,1) + nuz[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(1,2) + nuz[m_triangles[i].ver1]*h1_ver1(1,2) + nuz[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(1,0) + nuz[m_triangles[i].ver1]*h1_ver2(1,0) + nuz[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(1,1) + nuz[m_triangles[i].ver1]*h1_ver2(1,1) + nuz[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(1,2) + nuz[m_triangles[i].ver1]*h1_ver2(1,2) + nuz[m_triangles[i].ver2]*h2_ver2(1,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); ++ it) {
				mat_J(i*9+7+Start_id , it.index()) += m_TPPIBG[i].v0.y*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+2, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); ++ it) {
				mat_J(i*9+7+Start_id , it.index()) += m_TPPIBG[i].v1.y*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+2, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); ++ it) {
				mat_J(i*9+7+Start_id , it.index()) += m_TPPIBG[i].v2.y*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+2, it.index());
			}
			mat_f(i*9+7+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad2.y - (lambda_z[i].y/penParam + pz[i].y));

			// \partial grad2.z/\partial nz \times \partial nz/\partial u
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(2,0) + nuz[m_triangles[i].ver1]*h1_ver0(2,0) + nuz[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(2,1) + nuz[m_triangles[i].ver1]*h1_ver0(2,1) + nuz[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(2,2) + nuz[m_triangles[i].ver1]*h1_ver0(2,2) + nuz[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(2,0) + nuz[m_triangles[i].ver1]*h1_ver1(2,0) + nuz[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(2,1) + nuz[m_triangles[i].ver1]*h1_ver1(2,1) + nuz[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(2,2) + nuz[m_triangles[i].ver1]*h1_ver1(2,2) + nuz[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(2,0) + nuz[m_triangles[i].ver1]*h1_ver2(2,0) + nuz[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(2,1) + nuz[m_triangles[i].ver1]*h1_ver2(2,1) + nuz[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				area_scale*pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(2,2) + nuz[m_triangles[i].ver1]*h1_ver2(2,2) + nuz[m_triangles[i].ver2]*h2_ver2(2,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); ++ it) {
				mat_J(i*9+8+Start_id , it.index()) += m_TPPIBG[i].v0.z*area_scale*pen_scale*gradJ(m_triangles[i].ver0*3+2, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); ++ it) {
				mat_J(i*9+8+Start_id , it.index()) += m_TPPIBG[i].v1.z*area_scale*pen_scale*gradJ(m_triangles[i].ver1*3+2, it.index());
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); ++ it) {
				mat_J(i*9+8+Start_id , it.index()) += m_TPPIBG[i].v2.z*area_scale*pen_scale*gradJ(m_triangles[i].ver2*3+2, it.index());
			}
			mat_f(i*9+8+Start_id, 0) += area_scale*pen_scale*(m_triangles[i].grad2.z - (lambda_z[i].z/penParam + pz[i].z));
		}
		Start_id += m_trinum*9;
	}

	if (UseTVU) {
		// the tvu term, \gamma/2 \|(p+\lambda^k/\gamma)-\nabla u\|^2, this is defined on each triangle
		vector<double>ux, uy, uz; ux.resize(m_vnum); uy.resize(m_vnum); uz.resize(m_vnum);
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
			ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
			uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
			uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
		}
		GetPPIGradients (m_vnum,ux);
		GetPPIGradients1(m_vnum,uy);
		GetPPIGradients2(m_vnum,uz);
		for (int i = 0; i < m_trinum; ++ i) {
			double area_scale = UseFaceArea?sqrt(m_trianglesArea[i]):1.0;
			double cur_scale = pen_scale*area_scale;
			// need to calculate the derivative of m_TPPIBG[i] to each vertex;
			gmm::dense_matrix<double> h0_ver0(3,3), h0_ver1(3,3), h0_ver2(3,3), h1_ver0(3,3), h1_ver1(3,3), h1_ver2(3,3), h2_ver0(3,3), h2_ver1(3,3), h2_ver2(3,3), Mat_I(3,3);
			gmm::copy(gmm::identity_matrix(), Mat_I); 

			gmm::dense_matrix<double> pPi(3,1), pPj(3,1), pPk(3,1), Eij(3,1), eij(3,1), Eik(3,1), eijeijT(3,3), m_eijeijT(3,3), temp133(3,3), temp233(3,3), temp333(3,3), temp11(1,1);
			gmm::dense_matrix<double> Hk(3,1), hk(3,1), NHk(3,3);
			// calculate the derivative of h2
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk); //cout << NHk;

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h2_ver2); // calculate h22ver2
			gmm::mult(NHk, h2_ver2, temp333); gmm::copy(temp333, h2_ver2); //cout << h2_ver2;

			gmm::copy(Mat_I, h2_ver0); gmm::scale(h2_ver0, -1.0); gmm::add(eijeijT, h2_ver0); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h2_ver0, temp133, temp333);	gmm::copy(temp333, h2_ver0); //cout << h1_ver2;// calculate h12ver2
			gmm::mult(NHk, h2_ver0, temp333);		gmm::copy(temp333, h2_ver0); //cout << h2_ver0;

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij2vj
			gmm::mult(temp333, temp233, h2_ver1); //cout << h1_ver0;
			gmm::mult(NHk, h2_ver1, temp333);		gmm::copy(temp333, h2_ver1); //cout << h2_ver1;

			// calculate the derivative of h1
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk); //cout << NHk;

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h1_ver1); //cout<< h1_ver1;// calculate h12ver1
			gmm::mult(NHk, h1_ver1, temp333); gmm::copy(temp333, h1_ver1); //cout << h1_ver1;

			gmm::copy(Mat_I, h1_ver2); gmm::scale(h1_ver2, -1.0); gmm::add(eijeijT, h1_ver2); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h1_ver2, temp133, temp333);	gmm::copy(temp333, h1_ver2); //cout << h1_ver2;// calculate h12ver2
			gmm::mult(NHk, h1_ver2, temp333);		gmm::copy(temp333, h1_ver2); //cout << h1_ver2;

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij2vj
			gmm::mult(temp333, temp233, h1_ver0); //cout << h1_ver0;
			gmm::mult(NHk, h1_ver0, temp333);		gmm::copy(temp333, h1_ver0); //cout << h1_ver0;

			// calculate the derivative of h0
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/(gmm::mat_euclidean_norm_sqr(Hk)+epsilon));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk); //cout << NHk;

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h0_ver0); //cout << h0_ver0;// calculate h12ver1
			gmm::mult(NHk, h0_ver0, temp333); gmm::copy(temp333, h0_ver0); //cout << h0_ver0;

			gmm::copy(Mat_I, h0_ver1); gmm::scale(h0_ver1, -1.0); gmm::add(eijeijT, h0_ver1); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h0_ver1, temp133, temp333);	gmm::copy(temp333, h0_ver1); //cout << h0_ver1;// calculate h12ver2
			gmm::mult(NHk, h0_ver1, temp333);		gmm::copy(temp333, h0_ver1); //cout << h0_ver1;

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/(gmm::mat_euclidean_norm(Eij)+epsilon)); // calculate eij2vj
			gmm::mult(temp333, temp233, h0_ver2); //cout << h0_ver2;
			gmm::mult(NHk, h0_ver2, temp333);		gmm::copy(temp333, h0_ver2);  //cout << h0_ver2;
			

			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*0) += cur_scale*m_TPPIBG[i].v0.x + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(0,0) + ux[m_triangles[i].ver1]*h1_ver0(0,0) + ux[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(0,1) + ux[m_triangles[i].ver1]*h1_ver0(0,1) + ux[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(0,2) + ux[m_triangles[i].ver1]*h1_ver0(0,2) + ux[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*0) += cur_scale*m_TPPIBG[i].v1.x + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(0,0) + ux[m_triangles[i].ver1]*h1_ver1(0,0) + ux[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(0,1) + ux[m_triangles[i].ver1]*h1_ver1(0,1) + ux[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(0,2) + ux[m_triangles[i].ver1]*h1_ver1(0,2) + ux[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*0) += cur_scale*m_TPPIBG[i].v2.x + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(0,0) + ux[m_triangles[i].ver1]*h1_ver2(0,0) + ux[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(0,1) + ux[m_triangles[i].ver1]*h1_ver2(0,1) + ux[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(0,2) + ux[m_triangles[i].ver1]*h1_ver2(0,2) + ux[m_triangles[i].ver2]*h2_ver2(0,2));
			mat_f(i*9+0+Start_id, 0) += cur_scale*(m_triangles[i].grad.x - (lambda_x[i].x/penParam + px[i].x));
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*0) += cur_scale*m_TPPIBG[i].v0.y + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(1,0) + ux[m_triangles[i].ver1]*h1_ver0(1,0) + ux[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(1,1) + ux[m_triangles[i].ver1]*h1_ver0(1,1) + ux[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(1,2) + ux[m_triangles[i].ver1]*h1_ver0(1,2) + ux[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*0) += cur_scale*m_TPPIBG[i].v1.y +
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(1,0) + ux[m_triangles[i].ver1]*h1_ver1(1,0) + ux[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(1,1) + ux[m_triangles[i].ver1]*h1_ver1(1,1) + ux[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(1,2) + ux[m_triangles[i].ver1]*h1_ver1(1,2) + ux[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*0) += cur_scale*m_TPPIBG[i].v2.y + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(1,0) + ux[m_triangles[i].ver1]*h1_ver2(1,0) + ux[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(1,1) + ux[m_triangles[i].ver1]*h1_ver2(1,1) + ux[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(1,2) + ux[m_triangles[i].ver1]*h1_ver2(1,2) + ux[m_triangles[i].ver2]*h2_ver2(1,2));
			mat_f(i*9+1+Start_id, 0) += cur_scale*(m_triangles[i].grad.y - (lambda_x[i].y/penParam + px[i].y));
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*0) += cur_scale*m_TPPIBG[i].v0.z + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(2,0) + ux[m_triangles[i].ver1]*h1_ver0(2,0) + ux[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(2,1) + ux[m_triangles[i].ver1]*h1_ver0(2,1) + ux[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver0(2,2) + ux[m_triangles[i].ver1]*h1_ver0(2,2) + ux[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*0) += cur_scale*m_TPPIBG[i].v1.z +
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(2,0) + ux[m_triangles[i].ver1]*h1_ver1(2,0) + ux[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(2,1) + ux[m_triangles[i].ver1]*h1_ver1(2,1) + ux[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver1(2,2) + ux[m_triangles[i].ver1]*h1_ver1(2,2) + ux[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*0) += cur_scale*m_TPPIBG[i].v2.z + 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(2,0) + ux[m_triangles[i].ver1]*h1_ver2(2,0) + ux[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(2,1) + ux[m_triangles[i].ver1]*h1_ver2(2,1) + ux[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				cur_scale*(ux[m_triangles[i].ver0]*h0_ver2(2,2) + ux[m_triangles[i].ver1]*h1_ver2(2,2) + ux[m_triangles[i].ver2]*h2_ver2(2,2));
			mat_f(i*9+2+Start_id, 0) += cur_scale*(m_triangles[i].grad.z - (lambda_x[i].z/penParam + px[i].z));
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(0,0) + uy[m_triangles[i].ver1]*h1_ver0(0,0) + uy[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*1) += cur_scale*m_TPPIBG[i].v0.x + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(0,1) + uy[m_triangles[i].ver1]*h1_ver0(0,1) + uy[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(0,2) + uy[m_triangles[i].ver1]*h1_ver0(0,2) + uy[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(0,0) + uy[m_triangles[i].ver1]*h1_ver1(0,0) + uy[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*1) += cur_scale*m_TPPIBG[i].v1.x + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(0,1) + uy[m_triangles[i].ver1]*h1_ver1(0,1) + uy[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(0,2) + uy[m_triangles[i].ver1]*h1_ver1(0,2) + uy[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(0,0) + uy[m_triangles[i].ver1]*h1_ver2(0,0) + uy[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*1) += cur_scale*m_TPPIBG[i].v2.x + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(0,1) + uy[m_triangles[i].ver1]*h1_ver2(0,1) + uy[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(0,2) + uy[m_triangles[i].ver1]*h1_ver2(0,2) + uy[m_triangles[i].ver2]*h2_ver2(0,2));
			mat_f(i*9+3+Start_id, 0) += cur_scale*(m_triangles[i].grad1.x - (lambda_y[i].x/penParam + py[i].x));
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(1,0) + uy[m_triangles[i].ver1]*h1_ver0(1,0) + uy[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*1) += cur_scale*m_TPPIBG[i].v0.y + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(1,1) + uy[m_triangles[i].ver1]*h1_ver0(1,1) + uy[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(1,2) + uy[m_triangles[i].ver1]*h1_ver0(1,2) + uy[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(1,0) + uy[m_triangles[i].ver1]*h1_ver1(1,0) + uy[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*1) += cur_scale*m_TPPIBG[i].v1.y +
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(1,1) + uy[m_triangles[i].ver1]*h1_ver1(1,1) + uy[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(1,2) + uy[m_triangles[i].ver1]*h1_ver1(1,2) + uy[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(1,0) + uy[m_triangles[i].ver1]*h1_ver2(1,0) + uy[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*1) += cur_scale*m_TPPIBG[i].v2.y + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(1,1) + uy[m_triangles[i].ver1]*h1_ver2(1,1) + uy[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(1,2) + uy[m_triangles[i].ver1]*h1_ver2(1,2) + uy[m_triangles[i].ver2]*h2_ver2(1,2));
			mat_f(i*9+4+Start_id, 0) += cur_scale*(m_triangles[i].grad1.y - (lambda_y[i].y/penParam + py[i].y));
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(2,0) + uy[m_triangles[i].ver1]*h1_ver0(2,0) + uy[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*1) += cur_scale*m_TPPIBG[i].v0.z + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(2,1) + uy[m_triangles[i].ver1]*h1_ver0(2,1) + uy[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver0(2,2) + uy[m_triangles[i].ver1]*h1_ver0(2,2) + uy[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(2,0) + uy[m_triangles[i].ver1]*h1_ver1(2,0) + uy[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*1) += cur_scale*m_TPPIBG[i].v1.z +
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(2,1) + uy[m_triangles[i].ver1]*h1_ver1(2,1) + uy[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver1(2,2) + uy[m_triangles[i].ver1]*h1_ver1(2,2) + uy[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(2,0) + uy[m_triangles[i].ver1]*h1_ver2(2,0) + uy[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*1) += cur_scale*m_TPPIBG[i].v2.z + 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(2,1) + uy[m_triangles[i].ver1]*h1_ver2(2,1) + uy[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*2) += 
				cur_scale*(uy[m_triangles[i].ver0]*h0_ver2(2,2) + uy[m_triangles[i].ver1]*h1_ver2(2,2) + uy[m_triangles[i].ver2]*h2_ver2(2,2));
			mat_f(i*9+5+Start_id, 0) += cur_scale*(m_triangles[i].grad1.z - (lambda_y[i].z/penParam + py[i].z));
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(0,0) + uz[m_triangles[i].ver1]*h1_ver0(0,0) + uz[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(0,1) + uz[m_triangles[i].ver1]*h1_ver0(0,1) + uz[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*2) += cur_scale*m_TPPIBG[i].v0.x + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(0,2) + uz[m_triangles[i].ver1]*h1_ver0(0,2) + uz[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(0,0) + uz[m_triangles[i].ver1]*h1_ver1(0,0) + uz[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(0,1) + uz[m_triangles[i].ver1]*h1_ver1(0,1) + uz[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*2) += cur_scale*m_TPPIBG[i].v1.x + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(0,2) + uz[m_triangles[i].ver1]*h1_ver1(0,2) + uz[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(0,0) + uz[m_triangles[i].ver1]*h1_ver2(0,0) + uz[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(0,1) + uz[m_triangles[i].ver1]*h1_ver2(0,1) + uz[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*2) += cur_scale*m_TPPIBG[i].v2.x + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(0,2) + uz[m_triangles[i].ver1]*h1_ver2(0,2) + uz[m_triangles[i].ver2]*h2_ver2(0,2));
			mat_f(i*9+6+Start_id, 0) += cur_scale*(m_triangles[i].grad2.x - (lambda_z[i].x/penParam + pz[i].x));
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(1,0) + uz[m_triangles[i].ver1]*h1_ver0(1,0) + uz[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(1,1) + uz[m_triangles[i].ver1]*h1_ver0(1,1) + uz[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*2) += cur_scale*m_TPPIBG[i].v0.y + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(1,2) + uz[m_triangles[i].ver1]*h1_ver0(1,2) + uz[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(1,0) + uz[m_triangles[i].ver1]*h1_ver1(1,0) + uz[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(1,1) + uz[m_triangles[i].ver1]*h1_ver1(1,1) + uz[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*2) += cur_scale*m_TPPIBG[i].v1.y +
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(1,2) + uz[m_triangles[i].ver1]*h1_ver1(1,2) + uz[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(1,0) + uz[m_triangles[i].ver1]*h1_ver2(1,0) + uz[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(1,1) + uz[m_triangles[i].ver1]*h1_ver2(1,1) + uz[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*2) += cur_scale*m_TPPIBG[i].v2.y + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(1,2) + uz[m_triangles[i].ver1]*h1_ver2(1,2) + uz[m_triangles[i].ver2]*h2_ver2(1,2));
			mat_f(i*9+7+Start_id, 0) += cur_scale*(m_triangles[i].grad2.y - (lambda_z[i].y/penParam + pz[i].y));
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(2,0) + uz[m_triangles[i].ver1]*h1_ver0(2,0) + uz[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(2,1) + uz[m_triangles[i].ver1]*h1_ver0(2,1) + uz[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*2) += cur_scale*m_TPPIBG[i].v0.z + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver0(2,2) + uz[m_triangles[i].ver1]*h1_ver0(2,2) + uz[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(2,0) + uz[m_triangles[i].ver1]*h1_ver1(2,0) + uz[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(2,1) + uz[m_triangles[i].ver1]*h1_ver1(2,1) + uz[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*2) += cur_scale*m_TPPIBG[i].v1.z +
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver1(2,2) + uz[m_triangles[i].ver1]*h1_ver1(2,2) + uz[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*0) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(2,0) + uz[m_triangles[i].ver1]*h1_ver2(2,0) + uz[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*1) += 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(2,1) + uz[m_triangles[i].ver1]*h1_ver2(2,1) + uz[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*2) += cur_scale*m_TPPIBG[i].v2.z + 
				cur_scale*(uz[m_triangles[i].ver0]*h0_ver2(2,2) + uz[m_triangles[i].ver1]*h1_ver2(2,2) + uz[m_triangles[i].ver2]*h2_ver2(2,2));
			mat_f(i*9+8+Start_id, 0) += cur_scale*(m_triangles[i].grad2.z - (lambda_z[i].z/penParam + pz[i].z));
			////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////
		}
		Start_id += m_trinum*9;
	}
	return;
}

void TriangularMesh::LM_Testing(bool TestTV = false, int choice = 0)
{
	MyMesh T_Mesh = this->m_ObjTriMesh;
	//Test the partial normal equation using Levenberg–Marquardt Algorithm
	double sigma, eOld, eNew;
	double var = 2.0, epsilon1 = 1.0e-8, epsilon2 = 1.0e-8, minEdgeLength = 5.0e-3;
	int dimension = m_vnum * 3;
	std::vector<double> m_b(dimension, 0);
	std::vector<double> m_x(dimension, 0);
	std::vector<double> tau(dimension, 0);
	std::vector<OpenMesh::Vec3f> originalPhi(m_vnum);
	RowSparseMatrix mat_J;	std::vector<double> m_f;
	if (TestTV) {
		cout << "Testing the TV minimization process: " << choice << endl;
		this->LM_TV_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, choice);
	} else {
		cout << "Testing the normal minimization process: " << choice << endl;
		this->LM_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, choice);
	}

	cout <<"0: ";
	if (TestTV) {
		vector<double>nux, nuy, nuz; nux.resize(m_vnum); nuy.resize(m_vnum); nuz.resize(m_vnum);
		double total_energy = 0.0;
		if (choice > 0) {
			T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
			T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
				nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
				nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
				nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
			}
			double tvn_energy = CalculateVTVEnergy(T_Mesh, nux, nuy, nuz);
			cout << "The TV Norm energy is: " << tvn_energy; total_energy += tvn_energy;
		} 
		else {
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
				nux[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[0];
				nuy[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[1];
				nuz[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[2];
			}
			double tvn_energy = CalculateVTVNormEnergy(T_Mesh, nux, nuy, nuz);
			cout << "The TV energy is: " << tvn_energy; total_energy += tvn_energy;
		}
		cout << "; The V energy is: " << CalculateVEnergy(T_Mesh) << "; total is: " << total_energy+CalculateVEnergy(T_Mesh) << endl;
	} else {
		double total_energy = 0.0;
		if (choice < 4) {
			vector<double> vec_n_energy; vec_n_energy.resize(m_vnum); 
			cout << " N: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< CalculateNEnergy(T_Mesh,vec_n_energy);
			//<< CalculateVEnergy(T_Mesh) + CalculateNEnergy(T_Mesh, vec_n_energy);
			total_energy += CalculateNEnergy(T_Mesh,vec_n_energy);
			if (choice > 1 && choice < 4) {
				cout << " U: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< CalculateVEnergy(T_Mesh);
				total_energy += CalculateVEnergy(T_Mesh);
			}  
			// << ", The N energy is: " << CalculateNEnergy(T_Mesh,vec_n_energy) << "; total is: "
		}
		if (choice == 1 || choice == 3 || choice == 4) {
			cout << " ND: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< CalculateNDPEnergy(T_Mesh);
			total_energy += CalculateNDPEnergy(T_Mesh);
		}
		cout << "; Sum: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< total_energy;
	}

	eOld = gmm::vect_norm2_sqr(m_f);
	RowSparseMatrix mat_JTJ(dimension, dimension);
	RowSparseMatrix mat_A(dimension, dimension);
	gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
	for (int i = 0; i < dimension; i++) {
		tau[i] = mat_JTJ(i, i);
	}
	gmm::scale(tau, 1.e-5);
	gmm::mult(gmm::transposed(mat_J), m_f, m_b);
	gmm::scale(m_b, -1.0);
	double normG = gmm::vect_norminf(m_b);

	gmm::csc_matrix<double> JtJ;
	numc::SparseSolver solver;	double *b = new double[dimension];	double *vx = new double[dimension];
	OpenMesh::IO::Options write_options;	char buffer[255]; 
	int iteration = 0; int maxStep = 200;  
	while (iteration++ < maxStep) {
		cout << endl << iteration <<":";
		if (normG < epsilon1) {
			break;
		}
		for (int i = 0; i < dimension; i++) {
			mat_JTJ(i, i) += tau[i];
		}
		gmm::copy(mat_JTJ, JtJ);

		engEvalString( m_ep, "clear all;" );
		int nnz= gmm::nnz(JtJ);
		mxArray* arrayA=mxCreateSparse(JtJ.nrows(),JtJ.ncols(),nnz,mxREAL);
		double*  pr = mxGetPr( arrayA );
		mwIndex* ir = mxGetIr( arrayA );
		mwIndex* jc = mxGetJc( arrayA );
		int i=0;
		for(std::vector<double>::const_iterator it=JtJ.pr.begin(); it!=JtJ.pr.end();++it,i++)
			pr[i]=(*it);
		i=0;
		for(std::vector<unsigned int>::const_iterator it=JtJ.ir.begin(); it!=JtJ.ir.end();++it,i++)
			ir[i]=(*it);
		i=0;
		for(std::vector<unsigned int>::const_iterator it=JtJ.jc.begin(); it!=JtJ.jc.end();++it,i++)
			jc[i]=(*it);

		engPutVariable( m_ep, "A", arrayA );
		mxDestroyArray(arrayA);

		int m=m_b.size();
		mxArray* arrayb=mxCreateDoubleMatrix(m,1,mxREAL);
		double*  prb = mxGetPr( arrayb );
		for( int r=0;r<m;++r)
			prb[r]= m_b[r];
		engPutVariable( m_ep, "b", arrayb );
		mxDestroyArray(arrayb);

		engEvalString( m_ep, "x = A\\b;" );

		mxArray* arrayx = engGetVariable( m_ep, "x");
		m = mxGetDimensions( arrayx )[0];
		int n = mxGetDimensions( arrayx )[1];
		if( !arrayx || mxIsSparse( arrayx ) || mxIsComplex( arrayx ) || n!=1 )
		{
			mxDestroyArray( arrayx );
			cout << "Get x error"<<endl;
		}
		int mM= m_x.size();
		if( m != mM ){
			m_x.resize(m);
		}
		double*  prx = mxGetPr( arrayx );
		for( int r=0;r<m;++r)
			m_x[r]=prx[r];



		//// solve mat_A*m_x = m_b
		//numc::RowMat<double> RM_A(dimension, dimension); int r = 0;
		//for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(mat_A);
		//	rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(mat_A); ++ rIt, ++ r) {
		//	gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
		//	gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
		//	for (; it != ite; ++ it) {
		//		RM_A(r, it.index()) = mat_A(r, it.index());
		//	}
		//	b[r] = m_b[r];
		//}
		//solver.getMatA() = RM_A;
		//solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
		//solver.init();
		//// solve the nonlinear equation to calculate the new vertex position in vx;
		//solver.solve(b, vx);
		//solver.clear();
		//for (int r = 0; r < dimension; ++ r) {
		//	m_x[r] = vx[r];
		//}

		////modify m_x, do not change too much according to the value of MinEdgeLength
		//double maxChange = 0.0;
		//for (int i = 0; i < m_vnum; ++ i) {
		//	OpenMesh::Vec3f deltaP(m_x[i], m_x[i+m_vnum], m_x[i+m_vnum*2]);
		//	maxChange = (deltaP.length() > maxChange)?deltaP.length():maxChange;
		//}
		//double stepSize = (minEdgeLength/maxChange > 1.0)?1.0:(minEdgeLength/maxChange);
		//stepSize = 0.2;
		//gmm::scale(m_x, stepSize);
		////end of modify m_x

		double normH = 0.0, normX = 0.0, l0h = 0;
		normH = gmm::vect_norm2(m_x);
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			int idx = v_it.handle().idx();
			originalPhi[v_it.handle().idx()] = T_Mesh.point(v_it.handle());
			normX += originalPhi[v_it.handle().idx()].sqrnorm();
			for (int i = 0; i < 3; ++i) {
				l0h += 0.5*m_x[i*m_vnum + idx]*(tau[i*m_vnum + idx]*m_x[i*m_vnum + idx] + m_b[i*m_vnum + idx]);
			}
		}
		if (normH < epsilon2*(sqrt(normX) + epsilon2)) {
			break;
		} else {
			//change the coordinate value of each vertex...
			OpenMesh::Vec3f CurV;
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				int vIndex = v_it.handle().idx();
				T_Mesh.point(v_it.handle()) += OpenMesh::Vec3f(m_x[vIndex], m_x[m_vnum + vIndex], m_x[2*m_vnum + vIndex]);
			}
			UpdateMeshInfo(T_Mesh);
			if (iteration%1==0) {
				sprintf(buffer, "Results\\LMResults\\LMMeshResult_%d_%d.off", choice, iteration);
				if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
					std::cerr << "Cannot write mesh to file " << buffer << std::endl;
				}
			}
			if (TestTV) {
				vector<double>nux, nuy, nuz; nux.resize(m_vnum); nuy.resize(m_vnum); nuz.resize(m_vnum);
				double total_energy = 0.0;
				if (choice > 0) {
					T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
					T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
					for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
						nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
						nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
						nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
					}
					double tvn_energy = CalculateVTVEnergy(T_Mesh, nux, nuy, nuz);
					cout << "The TV Norm energy is: " << tvn_energy; total_energy += tvn_energy;
				} 
				else {
					for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
						nux[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[0];
						nuy[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[1];
						nuz[v_it.handle().idx()] = T_Mesh.point(v_it.handle()).data()[2];
					}
					double tvn_energy = CalculateVTVNormEnergy(T_Mesh, nux, nuy, nuz);
					cout << "The TV energy is: " << tvn_energy; total_energy += tvn_energy;
				}
				cout << "; The V energy is: " << CalculateVEnergy(T_Mesh) << "; total is: " << total_energy+CalculateVEnergy(T_Mesh) << endl;
			} else {
				double total_energy = 0.0;
				if (choice < 4) {
					vector<double> vec_n_energy; vec_n_energy.resize(m_vnum); 
					cout << " N: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< CalculateNEnergy(T_Mesh,vec_n_energy);
					//<< CalculateVEnergy(T_Mesh) + CalculateNEnergy(T_Mesh, vec_n_energy);
					total_energy += CalculateNEnergy(T_Mesh,vec_n_energy);
				if (choice > 1 && choice < 4) {
					cout << " U: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< CalculateVEnergy(T_Mesh);
					total_energy += CalculateVEnergy(T_Mesh);
				}  
					// << ", The N energy is: " << CalculateNEnergy(T_Mesh,vec_n_energy) << "; total is: "
				}
				if (choice == 1 || choice == 3 || choice == 4) {
					cout << " ND: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< CalculateNDPEnergy(T_Mesh);
					total_energy += CalculateNDPEnergy(T_Mesh);
				}
				cout << "; Sum: " <<setprecision(6)<<setw(10)<<setiosflags(ios::left)<< total_energy;
			}
			if (TestTV) {
				this->LM_TV_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, choice);
			} else {
				this->LM_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, choice);
			}
			eNew = gmm::vect_norm2_sqr(m_f);
			sigma = (eOld - eNew)/l0h;
			if (sigma >= 0) {
				gmm::scale(tau, std::max(1./3., 1. - pow(2.*sigma - 1., 3.)));
				var = 2.0;
				eOld = eNew;
				gmm::mult(gmm::transposed(mat_J), m_f, m_b);
				gmm::scale(m_b, -1.0);
				normG = gmm::vect_norminf(m_b);
				gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
				//if (iteration%1==0) {
				//	write_options.set(OpenMesh::IO::Options::VertexNormal); 
				//	sprintf(buffer, "Results\\LMResults\\LMMeshResult_%d_%d.off", choice, iteration);
				//	if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
				//		std::cerr << "Cannot write mesh to file " << buffer << std::endl;
				//	}
				//}
			} else {
				cout << "; Sigma = " << sigma ;
				gmm::scale(tau, var);
				var *= 2.0;
				for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
					T_Mesh.point(v_it) = originalPhi[v_it.handle().idx()];
				}
			}
		}
	}
	sprintf(buffer, "Results\\LMMeshResult_%d_Final.off", choice, iteration);
	if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
		std::cerr << "Cannot write mesh to file " << buffer << std::endl;
	}
	delete b, vx;
}

void TriangularMesh::LM_JacobianMatrix_Construction(MyMesh& T_Mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, int choice = 0)
{
	double epsilon = 1.0e-5;
	// construct the two matrix for Levenberg–Marquardt Algorithm
	if (choice == 0) {
		mat_J.resize(m_vnum*3, m_vnum*3);	mat_f.resize(m_vnum*3); // only norm term
	}
	if (choice == 1) {
		mat_J.resize(m_vnum*3+T_Mesh.n_edges()*3, m_vnum*3);	mat_f.resize(m_vnum*3+T_Mesh.n_edges()*3);
	}
	if (choice == 2) {
		mat_J.resize(m_vnum*6, m_vnum*3);	mat_f.resize(m_vnum*6); // norm term and u-f term
	}
	if (choice == 3) {
		mat_J.resize(m_vnum*6+T_Mesh.n_edges()*3, m_vnum*3);	mat_f.resize(m_vnum*6+T_Mesh.n_edges()*3);
	}
	if (choice == 4) {
		mat_J.resize(T_Mesh.n_edges()*3, m_vnum*3);	mat_f.resize(T_Mesh.n_edges()*3); // only n_i - n_j term
	}
	gmm::scale(mat_J, 0.0);		gmm::scale(mat_f, 0.0);

	//for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
	//	int vertex_id = v_it.handle().idx();
	//	//add the laplacian term for singularity mesh
	//	double degree = 0.0; OpenMesh::Vec3f avg_nei_point; avg_nei_point.vectorize(0.0);
	//	for (MyMesh::VertexVertexIter vv_it = T_Mesh.vv_iter(v_it); vv_it; ++vv_it) {
	//		degree++;
	//		avg_nei_point += T_Mesh.point(vv_it.handle());
	//	} avg_nei_point = avg_nei_point/degree;

	//	for (int k = 0; k < 3; ++ k) {
	//		mat_J(vertex_id*3+k+m_vnum*3, vertex_id+m_vnum*k) = 1.0;
	//		for (MyMesh::VertexVertexIter vv_it = T_Mesh.vv_iter(v_it); vv_it; ++vv_it) {
	//			mat_J(vertex_id*3+k+m_vnum*3, vv_it.handle().idx()+m_vnum*k) = -1.0/degree;
	//		}
	//		mat_f[vertex_id*3+k+m_vnum*3] = (T_Mesh.point(v_it.handle()).data()[k] - avg_nei_point[k]);
	//	}
	//}

	OpenMesh::Vec3f pointA , pointB , pointC; int ida, idb, idc;
	T_Mesh.update_face_normals(); T_Mesh.update_vertex_normals();
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
		// first for all neighbor triangle, store its vertex in counter clock wise, A(current vi), B, C
		vector< OpenMesh::Vec3f >	NeiVertexList;		NeiVertexList.clear();
		vector< int >				NeiVertexIDList;	NeiVertexIDList.clear();
		vector< double >			NeiFaceAreaList;	NeiFaceAreaList.clear();
		for (MyMesh::ConstVertexVertexIter vv_it = T_Mesh.cvv_iter(v_it.handle()); vv_it; ++ vv_it) { // the iterator is clockwise~~~~~~~
			NeiVertexIDList.push_back(vv_it.handle().idx());
			pointA = T_Mesh.point(vv_it.handle());
			NeiVertexList.push_back(pointA); 
		}
		std::reverse(NeiVertexIDList.begin(), NeiVertexIDList.end());
		std::reverse(NeiVertexList.begin(), NeiVertexList.end());

		//n(A) = pb*pc + pc*pd + pd*pe + pe*pf + pf*pg + pg*pb
		OpenMesh::Vec3f NA; NA.vectorize(0.0);
		for (int i = 0; i < NeiVertexList.size(); ++ i) { // need to add boundary condition judgment~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			NA += OpenMesh::cross(NeiVertexList[i], NeiVertexList[(i+1+NeiVertexList.size())%NeiVertexList.size()]);
		}
		double NAx = NA[0]; double NAy = NA[1]; double NAz = NA[2]; 
		double NAx_2 = NAx*NAx; double NAy_2 = NAy*NAy; double NAz_2 = NAz*NAz;
		double LNA = NA.norm()+epsilon; double LNA_2 = std::pow(LNA,2)+epsilon; double LNA_3 = std::pow(LNA,3)+epsilon; double LNA_4 = std::pow(LNA, 4)+epsilon;
		OpenMesh::Vec3f N_NA = NA; N_NA.normalize(); // normalized NA;
		OpenMesh::Vec3f T_NA(m_vertices[vertex_id].ps_normal_x,m_vertices[vertex_id].ps_normal_y,m_vertices[vertex_id].ps_normal_z); T_NA.normalize();

		for (int i = 0; i < NeiVertexList.size(); ++ i) {
			// to calculate the coefficient of point c based on the normal of point a, pb*pc+pc*pd
			OpenMesh::Vec3f pb = NeiVertexList[(i-1+NeiVertexList.size())%NeiVertexList.size()];
			OpenMesh::Vec3f pd = NeiVertexList[(i+1+NeiVertexList.size())%NeiVertexList.size()];
			OpenMesh::Vec3f pc = NeiVertexList[(i+0+NeiVertexList.size())%NeiVertexList.size()];
			double pbx = pb[0]; double pby = pb[1]; double pbz = pb[2];
			double pdx = pd[0]; double pdy = pd[1]; double pdz = pd[2];
			double pcx = pc[0]; double pcy = pc[1]; double pcz = pc[2];
			int cur_vid = NeiVertexIDList[i]; // for the center point pc

			OpenMesh::Vec3f Re, Rest_Norm; Re.vectorize(0.0);
			for (int j = 0; j < NeiVertexList.size(); ++ j) {
				if (j != i && (j+1+NeiVertexList.size())%NeiVertexList.size() != i) {
					//cout << j << "*" << (j+1+NeiVertexList.size())%NeiVertexList.size() << ", ";
					Re += OpenMesh::cross(NeiVertexList[j], NeiVertexList[(j+1+NeiVertexList.size())%NeiVertexList.size()]);
				}
			}
			double Rex = Re[0]; double Rey = Re[1]; double Rez = Re[2];

			if (choice < 4) {
				// JtN =(jacobian(N,C)./LNA - transpose(NA)*NA*jacobian(N,C)./LNA._3); // this is the correct normal
				mat_J(vertex_id*3+0, cur_vid+m_vnum*0) += (NAx*NAz*(pby - pdy) - NAx*NAy*(pbz - pdz))/(LNA_3);
				mat_J(vertex_id*3+0, cur_vid+m_vnum*1) += ((pbz - pdz)*NAx_2 - NAz*(pbx - pdx)*NAx)/LNA_3 - (pbz - pdz)/LNA;
				mat_J(vertex_id*3+0, cur_vid+m_vnum*2) += (pby - pdy)/LNA - ((pby - pdy)*NAx_2 - NAy*(pbx - pdx)*NAx)/LNA_3;

				mat_J(vertex_id*3+1, cur_vid+m_vnum*0) += (pbz - pdz)/LNA - ((pbz - pdz)*NAy_2 - NAz*(pby - pdy)*NAy)/LNA_3;
				mat_J(vertex_id*3+1, cur_vid+m_vnum*1) += -(NAy*NAz*(pbx - pdx) - NAx*NAy*(pbz - pdz))/LNA_3;
				mat_J(vertex_id*3+1, cur_vid+m_vnum*2) += ((pbx - pdx)*NAy_2 - NAx*(pby - pdy)*NAy)/LNA_3 - (pbx - pdx)/LNA;

				mat_J(vertex_id*3+2, cur_vid+m_vnum*0) += ((pby - pdy)*NAz_2 - NAy*(pbz - pdz)*NAz)/LNA_3 - (pby - pdy)/LNA;
				mat_J(vertex_id*3+2, cur_vid+m_vnum*1) += (pbx - pdx)/LNA - ((pbx - pdx)*NAz_2 - NAx*(pbz - pdz)*NAz)/LNA_3;
				mat_J(vertex_id*3+2, cur_vid+m_vnum*2) += (NAy*NAz*(pbx - pdx) - NAx*NAz*(pby - pdy))/LNA_3;
			}
		}
		if (choice < 4) {
			mat_f[vertex_id*3+0] += (T_Mesh.normal(v_it.handle()).data()[0] - m_vertices[vertex_id].ps_normal_x);
			mat_f[vertex_id*3+1] += (T_Mesh.normal(v_it.handle()).data()[1] - m_vertices[vertex_id].ps_normal_y);
			mat_f[vertex_id*3+2] += (T_Mesh.normal(v_it.handle()).data()[2] - m_vertices[vertex_id].ps_normal_z);
		}
		
		if (choice > 1 && choice < 4) {
			mat_J(vertex_id*3+0+m_vnum*3, vertex_id+m_vnum*0) += 1.0;
			mat_J(vertex_id*3+1+m_vnum*3, vertex_id+m_vnum*1) += 1.0;
			mat_J(vertex_id*3+2+m_vnum*3, vertex_id+m_vnum*2) += 1.0;

			mat_f[vertex_id*3+0+m_vnum*3] += (T_Mesh.point(v_it.handle()).data()[0] - m_vertices[vertex_id].ref_x);
			mat_f[vertex_id*3+1+m_vnum*3] += (T_Mesh.point(v_it.handle()).data()[1] - m_vertices[vertex_id].ref_y);
			mat_f[vertex_id*3+2+m_vnum*3] += (T_Mesh.point(v_it.handle()).data()[2] - m_vertices[vertex_id].ref_z);
		}
	}
	int Start_id = 0;
	if (choice == 1) {
		Start_id = m_vnum*3;
	}
	if (choice == 3) {
		Start_id = m_vnum*6;
	}
	if (choice == 1 || choice == 3 || choice == 4) {
		RowSparseMatrix gradJ; std::vector<double> grad_f;
		this->LM_JacobianMatrix_Construction(T_Mesh, gradJ, grad_f, 0);
		//fstream nof("Results\\grad_J.txt",std::ios::out);
		//if (!nof) {
		//	cout << "Can not open txt file to save grad_J data..." << endl;
		//	return;
		//}
		//int r = 0;
		//for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(gradJ);
		//	rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(gradJ); ++ rIt, ++ r) {
		//		nof << r << ": ";
		//		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
		//		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
		//		for (; it != ite; ++ it) {
		//			nof << it.index() << "," << gradJ(r, it.index()) << "; ";
		//		}
		//		nof << endl;
		//}
		//nof.close();
		// the normal difference term, eta*\omega_ij\|n_i - n_j\|^2
		T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
		T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
			MyMesh::HalfedgeHandle h1 = T_Mesh.halfedge_handle(e_it,0);
			MyMesh::HalfedgeHandle h2 = T_Mesh.halfedge_handle(e_it,1);
			int iid = T_Mesh.to_vertex_handle(h1).idx();
			int jid = T_Mesh.to_vertex_handle(h2).idx();

			for (int k = 0; k < 3; k++) {
				//cout << endl << iid*3+k << " row of grad_J: ";
				for (it = vect_const_begin(gmm::mat_const_row(gradJ, iid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, iid*3+k)); ++ it) {
				//	cout << it.index() << ": " << gradJ(iid*3+k, it.index()) << ", ";
					mat_J(e_it.handle().idx()*3+k+Start_id, it.index()) += gradJ(iid*3+k, it.index());
				}
				//cout << endl<<endl << e_it.handle().idx()*3+k << " row of mat_J: ";
				//for (it = vect_const_begin(gmm::mat_const_row(mat_J, e_it.handle().idx()*3+k)); it != vect_const_end(gmm::mat_const_row(mat_J, e_it.handle().idx()*3+k)); ++ it) {
				//	cout << it.index() << ": " << mat_J(e_it.handle().idx()*3+k, it.index()) << ", ";
				//}

				//cout << endl<<endl << jid*3+k << " row of grad_J: ";
				for (it = vect_const_begin(gmm::mat_const_row(gradJ, jid*3+k)); it != vect_const_end(gmm::mat_const_row(gradJ, jid*3+k)); ++ it) {
					//cout << it.index() << ": " << -1.0*gradJ(jid*3+k, it.index()) << ", ";
					mat_J(e_it.handle().idx()*3+k+Start_id, it.index()) -= gradJ(jid*3+k, it.index());
				}
				//cout << endl<<endl << e_it.handle().idx()*3+k << " row of mat_J: ";
				//for (it = vect_const_begin(gmm::mat_const_row(mat_J, e_it.handle().idx()*3+k)); it != vect_const_end(gmm::mat_const_row(mat_J, e_it.handle().idx()*3+k)); ++ it) {
				//	cout << it.index() << ": " << mat_J(e_it.handle().idx()*3+k, it.index()) << ", ";
				//}

				mat_f[e_it.handle().idx()*3+k+Start_id] = (T_Mesh.normal(MyMesh::VertexHandle(iid)).data()[k] - T_Mesh.normal(MyMesh::VertexHandle(jid)).data()[k]);
			}

		}
	}
	return;
}

void TriangularMesh::LM_TV_JacobianMatrix_Construction(MyMesh& T_Mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, bool UseTVNorm = false)
{
	// only test TVU and TVNorm
	RowSparseMatrix gradJ;  std::vector<double> fx;
	this->LM_JacobianMatrix_Construction(T_Mesh, gradJ, fx, 0);
	mat_J.resize(m_trinum*9+m_vnum*3, m_vnum*3); 	mat_f.resize(m_trinum*9+m_vnum*3);
	double pen_scale = 1.0, beta_scale = 1.0;
	int Start_id = 0;

	if (UseTVNorm) {
		// the tv term, \gamma/2 \|(p+\lambda^k/\gamma)-\nabla n(u)\|^2, this is defined on each triangle
		vector<double>nux, nuy, nuz; nux.resize(m_vnum); nuy.resize(m_vnum); nuz.resize(m_vnum);
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			m_vertices[v_it.handle().idx()].x = T_Mesh.point(v_it.handle()).data()[0];
			m_vertices[v_it.handle().idx()].y = T_Mesh.point(v_it.handle()).data()[1];
			m_vertices[v_it.handle().idx()].z = T_Mesh.point(v_it.handle()).data()[2];
		}
		BuildTPPIBG();
		T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
		T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
			nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
			nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
			nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
		}
		GetPPIGradients (m_vnum,nux);
		GetPPIGradients1(m_vnum,nuy);
		GetPPIGradients2(m_vnum,nuz);
		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
		for (int i = 0; i < m_trinum; ++ i) {
			// need to calculate the derivative of m_TPPIBG[i] to each vertex;
			gmm::dense_matrix<double> h0_ver0(3,3), h0_ver1(3,3), h0_ver2(3,3), h1_ver0(3,3), h1_ver1(3,3), h1_ver2(3,3), h2_ver0(3,3), h2_ver1(3,3), h2_ver2(3,3), Mat_I(3,3);
			gmm::copy(gmm::identity_matrix(), Mat_I); 

			gmm::dense_matrix<double> pPi(3,1), pPj(3,1), pPk(3,1), Eij(3,1), eij(3,1), Eik(3,1), eijeijT(3,3), m_eijeijT(3,3), temp133(3,3), temp233(3,3), temp333(3,3), temp11(1,1);
			// calculate the derivative of h2
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij

			gmm::dense_matrix<double> Hk(3,1), hk(3,1), NHk(3,3);
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk);

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h2_ver2); //cout<<h2_ver2;// calculate h22ver2
			gmm::mult(NHk, h2_ver2, temp333); gmm::copy(temp333, h2_ver2);

			gmm::copy(Mat_I, h2_ver0); gmm::scale(h2_ver0, -1.0);	gmm::add(eijeijT, h2_ver0); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11);			gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233);  // calculate eijEikT
			gmm::add(temp133, temp233, temp333);   // calcualte dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij));				gmm::scale(temp233, -1.0); // calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h2_ver0, temp133, temp333);	gmm::copy(temp333, h2_ver0); //cout << h2_ver0;// calculate h22ver0
			gmm::mult(NHk, h2_ver0, temp333);		gmm::copy(temp333, h2_ver0);

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij2vj
			gmm::mult(temp333, temp233, h2_ver1); //cout << h2_ver1; 
			gmm::mult(NHk, h2_ver1, temp333);		gmm::copy(temp333, h2_ver1);

			// calculate the derivative of h1
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk);

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h1_ver1); //cout<< h1_ver1;// calculate h12ver1
			gmm::mult(NHk, h1_ver1, temp333); gmm::copy(temp333, h1_ver1);

			gmm::copy(Mat_I, h1_ver2); gmm::scale(h1_ver2, -1.0); gmm::add(eijeijT, h1_ver2); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h1_ver2, temp133, temp333);	gmm::copy(temp333, h1_ver2); //cout << h1_ver2;// calculate h12ver2
			gmm::mult(NHk, h1_ver2, temp333);		gmm::copy(temp333, h1_ver2);

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij2vj
			gmm::mult(temp333, temp233, h1_ver0); //cout << h1_ver0;
			gmm::mult(NHk, h1_ver0, temp333);		gmm::copy(temp333, h1_ver0);

			// calculate the derivative of h0
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
			}
			gmm::copy(pPi, Eij);	gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik);	gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk);

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h0_ver0); //cout << h0_ver0;// calculate h12ver1
			gmm::mult(NHk, h0_ver0, temp333); gmm::copy(temp333, h0_ver0);

			gmm::copy(Mat_I, h0_ver1); gmm::scale(h0_ver1, -1.0); gmm::add(eijeijT, h0_ver1); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calcualte dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij));				gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133);	gmm::scale(temp133, -1.0);	gmm::add(h0_ver1, temp133, temp333);	gmm::copy(temp333, h0_ver1); //cout << h0_ver1;// calculate h12ver2
			gmm::mult(NHk, h0_ver1, temp333);		gmm::copy(temp333, h0_ver1);

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);	gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233);	gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij2vj
			gmm::mult(temp333, temp233, h0_ver2); //cout << h0_ver2;
			gmm::mult(NHk, h0_ver2, temp333);		gmm::copy(temp333, h0_ver2);


			// \partial grad.x/\partial nx \times \partial nx/\partial u
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(0,0) + nux[m_triangles[i].ver1]*h1_ver0(0,0) + nux[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(0,1) + nux[m_triangles[i].ver1]*h1_ver0(0,1) + nux[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(0,2) + nux[m_triangles[i].ver1]*h1_ver0(0,2) + nux[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(0,0) + nux[m_triangles[i].ver1]*h1_ver1(0,0) + nux[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(0,1) + nux[m_triangles[i].ver1]*h1_ver1(0,1) + nux[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(0,2) + nux[m_triangles[i].ver1]*h1_ver1(0,2) + nux[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(0,0) + nux[m_triangles[i].ver1]*h1_ver2(0,0) + nux[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(0,1) + nux[m_triangles[i].ver1]*h1_ver2(0,1) + nux[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(0,2) + nux[m_triangles[i].ver1]*h1_ver2(0,2) + nux[m_triangles[i].ver2]*h2_ver2(0,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); ++ it) {
				mat_J(i*9+0+Start_id , it.index()) += m_TPPIBG[i].v0.x*pen_scale*gradJ(m_triangles[i].ver0*3+0, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); ++ it) {
				mat_J(i*9+0+Start_id , it.index()) += m_TPPIBG[i].v1.x*pen_scale*gradJ(m_triangles[i].ver1*3+0, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); ++ it) {
				mat_J(i*9+0+Start_id , it.index()) += m_TPPIBG[i].v2.x*pen_scale*gradJ(m_triangles[i].ver2*3+0, it.index())/beta_scale;
			}
			mat_f[i*9+0+Start_id ] = pen_scale*(m_triangles[i].grad.x);

			// \partial grad.y/\partial nx \times \partial nx/\partial u
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(1,0) + nux[m_triangles[i].ver1]*h1_ver0(1,0) + nux[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(1,1) + nux[m_triangles[i].ver1]*h1_ver0(1,1) + nux[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(1,2) + nux[m_triangles[i].ver1]*h1_ver0(1,2) + nux[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(1,0) + nux[m_triangles[i].ver1]*h1_ver1(1,0) + nux[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(1,1) + nux[m_triangles[i].ver1]*h1_ver1(1,1) + nux[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(1,2) + nux[m_triangles[i].ver1]*h1_ver1(1,2) + nux[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(1,0) + nux[m_triangles[i].ver1]*h1_ver2(1,0) + nux[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(1,1) + nux[m_triangles[i].ver1]*h1_ver2(1,1) + nux[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(1,2) + nux[m_triangles[i].ver1]*h1_ver2(1,2) + nux[m_triangles[i].ver2]*h2_ver2(1,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); ++ it) {
				mat_J(i*9+1+Start_id , it.index()) += m_TPPIBG[i].v0.y*pen_scale*gradJ(m_triangles[i].ver0*3+0, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); ++ it) {
				mat_J(i*9+1+Start_id , it.index()) += m_TPPIBG[i].v1.y*pen_scale*gradJ(m_triangles[i].ver1*3+0, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); ++ it) {
				mat_J(i*9+1+Start_id , it.index()) += m_TPPIBG[i].v2.y*pen_scale*gradJ(m_triangles[i].ver2*3+0, it.index())/beta_scale;
			}
			mat_f[i*9+1+Start_id ] = pen_scale*(m_triangles[i].grad.y);

			// \partial grad.z/\partial nx \times \partial nx/\partial u
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(2,0) + nux[m_triangles[i].ver1]*h1_ver0(2,0) + nux[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(2,1) + nux[m_triangles[i].ver1]*h1_ver0(2,1) + nux[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver0(2,2) + nux[m_triangles[i].ver1]*h1_ver0(2,2) + nux[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(2,0) + nux[m_triangles[i].ver1]*h1_ver1(2,0) + nux[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(2,1) + nux[m_triangles[i].ver1]*h1_ver1(2,1) + nux[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver1(2,2) + nux[m_triangles[i].ver1]*h1_ver1(2,2) + nux[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(2,0) + nux[m_triangles[i].ver1]*h1_ver2(2,0) + nux[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(2,1) + nux[m_triangles[i].ver1]*h1_ver2(2,1) + nux[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nux[m_triangles[i].ver0]*h0_ver2(2,2) + nux[m_triangles[i].ver1]*h1_ver2(2,2) + nux[m_triangles[i].ver2]*h2_ver2(2,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+0)); ++ it) {
				mat_J(i*9+2+Start_id , it.index()) += m_TPPIBG[i].v0.z*pen_scale*gradJ(m_triangles[i].ver0*3+0, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+0)); ++ it) {
				mat_J(i*9+2+Start_id , it.index()) += m_TPPIBG[i].v1.z*pen_scale*gradJ(m_triangles[i].ver1*3+0, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+0)); ++ it) {
				mat_J(i*9+2+Start_id , it.index()) += m_TPPIBG[i].v2.z*pen_scale*gradJ(m_triangles[i].ver2*3+0, it.index())/beta_scale;
			}
			mat_f[i*9+2+Start_id ] = pen_scale*(m_triangles[i].grad.z);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// \partial grad1.x/\partial ny \times \partial ny/\partial u
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(0,0) + nuy[m_triangles[i].ver1]*h1_ver0(0,0) + nuy[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(0,1) + nuy[m_triangles[i].ver1]*h1_ver0(0,1) + nuy[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(0,2) + nuy[m_triangles[i].ver1]*h1_ver0(0,2) + nuy[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(0,0) + nuy[m_triangles[i].ver1]*h1_ver1(0,0) + nuy[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*1) =  
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(0,1) + nuy[m_triangles[i].ver1]*h1_ver1(0,1) + nuy[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(0,2) + nuy[m_triangles[i].ver1]*h1_ver1(0,2) + nuy[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(0,0) + nuy[m_triangles[i].ver1]*h1_ver2(0,0) + nuy[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(0,1) + nuy[m_triangles[i].ver1]*h1_ver2(0,1) + nuy[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(0,2) + nuy[m_triangles[i].ver1]*h1_ver2(0,2) + nuy[m_triangles[i].ver2]*h2_ver2(0,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); ++ it) {
				mat_J(i*9+3+Start_id , it.index()) += m_TPPIBG[i].v0.x*pen_scale*gradJ(m_triangles[i].ver0*3+1, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); ++ it) {
				mat_J(i*9+3+Start_id , it.index()) += m_TPPIBG[i].v1.x*pen_scale*gradJ(m_triangles[i].ver1*3+1, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); ++ it) {
				mat_J(i*9+3+Start_id , it.index()) += m_TPPIBG[i].v2.x*pen_scale*gradJ(m_triangles[i].ver2*3+1, it.index())/beta_scale;
			}
			mat_f[i*9+3+Start_id ] = pen_scale*(m_triangles[i].grad1.x);

			// \partial grad1.y/\partial ny \times \partial ny/\partial u
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(1,0) + nuy[m_triangles[i].ver1]*h1_ver0(1,0) + nuy[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(1,1) + nuy[m_triangles[i].ver1]*h1_ver0(1,1) + nuy[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(1,2) + nuy[m_triangles[i].ver1]*h1_ver0(1,2) + nuy[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(1,0) + nuy[m_triangles[i].ver1]*h1_ver1(1,0) + nuy[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(1,1) + nuy[m_triangles[i].ver1]*h1_ver1(1,1) + nuy[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(1,2) + nuy[m_triangles[i].ver1]*h1_ver1(1,2) + nuy[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(1,0) + nuy[m_triangles[i].ver1]*h1_ver2(1,0) + nuy[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(1,1) + nuy[m_triangles[i].ver1]*h1_ver2(1,1) + nuy[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(1,2) + nuy[m_triangles[i].ver1]*h1_ver2(1,2) + nuy[m_triangles[i].ver2]*h2_ver2(1,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); ++ it) {
				mat_J(i*9+4+Start_id , it.index()) += m_TPPIBG[i].v0.y*pen_scale*gradJ(m_triangles[i].ver0*3+1, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); ++ it) {
				mat_J(i*9+4+Start_id , it.index()) += m_TPPIBG[i].v1.y*pen_scale*gradJ(m_triangles[i].ver1*3+1, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); ++ it) {
				mat_J(i*9+4+Start_id , it.index()) += m_TPPIBG[i].v2.y*pen_scale*gradJ(m_triangles[i].ver2*3+1, it.index())/beta_scale;
			}
			mat_f[i*9+4+Start_id ] = pen_scale*(m_triangles[i].grad1.y);

			// \partial grad1.z/\partial ny \times \partial ny/\partial u
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(2,0) + nuy[m_triangles[i].ver1]*h1_ver0(2,0) + nuy[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(2,1) + nuy[m_triangles[i].ver1]*h1_ver0(2,1) + nuy[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver0(2,2) + nuy[m_triangles[i].ver1]*h1_ver0(2,2) + nuy[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(2,0) + nuy[m_triangles[i].ver1]*h1_ver1(2,0) + nuy[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(2,1) + nuy[m_triangles[i].ver1]*h1_ver1(2,1) + nuy[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver1(2,2) + nuy[m_triangles[i].ver1]*h1_ver1(2,2) + nuy[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(2,0) + nuy[m_triangles[i].ver1]*h1_ver2(2,0) + nuy[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(2,1) + nuy[m_triangles[i].ver1]*h1_ver2(2,1) + nuy[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nuy[m_triangles[i].ver0]*h0_ver2(2,2) + nuy[m_triangles[i].ver1]*h1_ver2(2,2) + nuy[m_triangles[i].ver2]*h2_ver2(2,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+1)); ++ it) {
				mat_J(i*9+5+Start_id , it.index()) += m_TPPIBG[i].v0.z*pen_scale*gradJ(m_triangles[i].ver0*3+1, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+1)); ++ it) {
				mat_J(i*9+5+Start_id , it.index()) += m_TPPIBG[i].v1.z*pen_scale*gradJ(m_triangles[i].ver1*3+1, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+1)); ++ it) {
				mat_J(i*9+5+Start_id , it.index()) += m_TPPIBG[i].v2.z*pen_scale*gradJ(m_triangles[i].ver2*3+1, it.index())/beta_scale;
			}
			mat_f[i*9+5+Start_id ] = pen_scale*(m_triangles[i].grad1.z);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// \partial grad2.x/\partial nz \times \partial nz/\partial u
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(0,0) + nuz[m_triangles[i].ver1]*h1_ver0(0,0) + nuz[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(0,1) + nuz[m_triangles[i].ver1]*h1_ver0(0,1) + nuz[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(0,2) + nuz[m_triangles[i].ver1]*h1_ver0(0,2) + nuz[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(0,0) + nuz[m_triangles[i].ver1]*h1_ver1(0,0) + nuz[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(0,1) + nuz[m_triangles[i].ver1]*h1_ver1(0,1) + nuz[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(0,2) + nuz[m_triangles[i].ver1]*h1_ver1(0,2) + nuz[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(0,0) + nuz[m_triangles[i].ver1]*h1_ver2(0,0) + nuz[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(0,1) + nuz[m_triangles[i].ver1]*h1_ver2(0,1) + nuz[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(0,2) + nuz[m_triangles[i].ver1]*h1_ver2(0,2) + nuz[m_triangles[i].ver2]*h2_ver2(0,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); ++ it) {
				mat_J(i*9+6+Start_id , it.index()) += m_TPPIBG[i].v0.x*pen_scale*gradJ(m_triangles[i].ver0*3+2, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); ++ it) {
				mat_J(i*9+6+Start_id , it.index()) += m_TPPIBG[i].v1.x*pen_scale*gradJ(m_triangles[i].ver1*3+2, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); ++ it) {
				mat_J(i*9+6+Start_id , it.index()) += m_TPPIBG[i].v2.x*pen_scale*gradJ(m_triangles[i].ver2*3+2, it.index())/beta_scale;
			}
			mat_f[i*9+6+Start_id ] = pen_scale*(m_triangles[i].grad2.x);

			// \partial grad2.y/\partial nz \times \partial nz/\partial u
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(1,0) + nuz[m_triangles[i].ver1]*h1_ver0(1,0) + nuz[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(1,1) + nuz[m_triangles[i].ver1]*h1_ver0(1,1) + nuz[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(1,2) + nuz[m_triangles[i].ver1]*h1_ver0(1,2) + nuz[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(1,0) + nuz[m_triangles[i].ver1]*h1_ver1(1,0) + nuz[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(1,1) + nuz[m_triangles[i].ver1]*h1_ver1(1,1) + nuz[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(1,2) + nuz[m_triangles[i].ver1]*h1_ver1(1,2) + nuz[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(1,0) + nuz[m_triangles[i].ver1]*h1_ver2(1,0) + nuz[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(1,1) + nuz[m_triangles[i].ver1]*h1_ver2(1,1) + nuz[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(1,2) + nuz[m_triangles[i].ver1]*h1_ver2(1,2) + nuz[m_triangles[i].ver2]*h2_ver2(1,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); ++ it) {
				mat_J(i*9+7+Start_id , it.index()) += m_TPPIBG[i].v0.y*pen_scale*gradJ(m_triangles[i].ver0*3+2, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); ++ it) {
				mat_J(i*9+7+Start_id , it.index()) += m_TPPIBG[i].v1.y*pen_scale*gradJ(m_triangles[i].ver1*3+2, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); ++ it) {
				mat_J(i*9+7+Start_id , it.index()) += m_TPPIBG[i].v2.y*pen_scale*gradJ(m_triangles[i].ver2*3+2, it.index())/beta_scale;
			}
			mat_f[i*9+7+Start_id ] = pen_scale*(m_triangles[i].grad2.y);

			// \partial grad2.z/\partial nz \times \partial nz/\partial u
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(2,0) + nuz[m_triangles[i].ver1]*h1_ver0(2,0) + nuz[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(2,1) + nuz[m_triangles[i].ver1]*h1_ver0(2,1) + nuz[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver0(2,2) + nuz[m_triangles[i].ver1]*h1_ver0(2,2) + nuz[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(2,0) + nuz[m_triangles[i].ver1]*h1_ver1(2,0) + nuz[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(2,1) + nuz[m_triangles[i].ver1]*h1_ver1(2,1) + nuz[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver1(2,2) + nuz[m_triangles[i].ver1]*h1_ver1(2,2) + nuz[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(2,0) + nuz[m_triangles[i].ver1]*h1_ver2(2,0) + nuz[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(2,1) + nuz[m_triangles[i].ver1]*h1_ver2(2,1) + nuz[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(nuz[m_triangles[i].ver0]*h0_ver2(2,2) + nuz[m_triangles[i].ver1]*h1_ver2(2,2) + nuz[m_triangles[i].ver2]*h2_ver2(2,2));
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver0*3+2)); ++ it) {
				mat_J(i*9+8+Start_id , it.index()) += m_TPPIBG[i].v0.z*pen_scale*gradJ(m_triangles[i].ver0*3+2, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver1*3+2)); ++ it) {
				mat_J(i*9+8+Start_id , it.index()) += m_TPPIBG[i].v1.z*pen_scale*gradJ(m_triangles[i].ver1*3+2, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(gradJ, m_triangles[i].ver2*3+2)); ++ it) {
				mat_J(i*9+8+Start_id , it.index()) += m_TPPIBG[i].v2.z*pen_scale*gradJ(m_triangles[i].ver2*3+2, it.index())/beta_scale;
			}
			mat_f[i*9+8+Start_id ] = pen_scale*(m_triangles[i].grad2.z);
		}
	} 
	else {
		// the tvu term, \gamma/2 \|(p+\lambda^k/\gamma)-\nabla u\|^2, this is defined on each triangle
		vector<double>ux, uy, uz; ux.resize(m_vnum); uy.resize(m_vnum); uz.resize(m_vnum);
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
			ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
			uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
			uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
		}
		GetPPIGradients (m_vnum,ux);
		GetPPIGradients1(m_vnum,uy);
		GetPPIGradients2(m_vnum,uz);
		for (int i = 0; i < m_trinum; ++ i) {
			// need to calculate the derivative of m_TPPIBG[i] to each vertex;
			gmm::dense_matrix<double> h0_ver0(3,3), h0_ver1(3,3), h0_ver2(3,3), h1_ver0(3,3), h1_ver1(3,3), h1_ver2(3,3), h2_ver0(3,3), h2_ver1(3,3), h2_ver2(3,3), Mat_I(3,3);
			gmm::copy(gmm::identity_matrix(), Mat_I); 

			gmm::dense_matrix<double> pPi(3,1), pPj(3,1), pPk(3,1), Eij(3,1), eij(3,1), Eik(3,1), eijeijT(3,3), m_eijeijT(3,3), temp133(3,3), temp233(3,3), temp333(3,3), temp11(1,1);
			// calculate the derivative of h2
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
			}
			gmm::copy(pPi, Eij); gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik); gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij

			gmm::dense_matrix<double> Hk(3,1), hk(3,1), NHk(3,3);
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk);

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h2_ver2); //cout<<h2_ver2;// calculate h22ver2
			gmm::mult(NHk, h2_ver2, temp333); gmm::copy(temp333, h2_ver2);

			gmm::copy(Mat_I, h2_ver0); gmm::scale(h2_ver0, -1.0); gmm::add(eijeijT, h2_ver0); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233);  // calculate eijEikT
			gmm::add(temp133, temp233, temp333);   // calcualte dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233); gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); gmm::scale(temp233, -1.0); // calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133); gmm::scale(temp133, -1.0);	gmm::add(h2_ver0, temp133, temp333); gmm::copy(temp333, h2_ver0); //cout << h2_ver0;// calculate h22ver0
			gmm::mult(NHk, h2_ver0, temp333); gmm::copy(temp333, h2_ver0);

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333); gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233); gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij2vj
			gmm::mult(temp333, temp233, h2_ver1); //cout << h2_ver1; 
			gmm::mult(NHk, h2_ver1, temp333); gmm::copy(temp333, h2_ver1);

			// calculate the derivative of h1
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
			}
			gmm::copy(pPi, Eij); gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik); gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk);

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h1_ver1); //cout<< h1_ver1;// calculate h12ver1
			gmm::mult(NHk, h1_ver1, temp333); gmm::copy(temp333, h1_ver1);

			gmm::copy(Mat_I, h1_ver2); gmm::scale(h1_ver2, -1.0); gmm::add(eijeijT, h1_ver2); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calculate dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233); gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133); gmm::scale(temp133, -1.0);	gmm::add(h1_ver2, temp133, temp333); gmm::copy(temp333, h1_ver2); //cout << h1_ver2;// calculate h12ver2
			gmm::mult(NHk, h1_ver2, temp333); gmm::copy(temp333, h1_ver2);
			
			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333); gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233); gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij2vj
			gmm::mult(temp333, temp233, h1_ver0); //cout << h1_ver0;
			gmm::mult(NHk, h1_ver0, temp333); gmm::copy(temp333, h1_ver0);

			// calculate the derivative of h0
			for (int idx = 0; idx < 3; ++ idx) {
				pPi(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver1)).data()[idx];
				pPj(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver2)).data()[idx];
				pPk(idx, 0) = T_Mesh.point(MyMesh::VertexHandle(m_triangles[i].ver0)).data()[idx];
			}
			gmm::copy(pPi, Eij); gmm::scale(Eij, -1.0); gmm::add(pPj, Eij); // calculate Eij
			gmm::copy(pPi, Eik); gmm::scale(Eik, -1.0); gmm::add(pPk, Eik); // calculate Eik
			gmm::copy(Eij, eij);	gmm::scale(eij, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij

			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(eij, hk); gmm::scale(hk, -1.0*temp11(0,0)); gmm::add(Eik, hk, Hk);
			gmm::copy(Hk, hk); gmm::scale(hk, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::copy(Mat_I, temp133); gmm::scale(temp133, 1/gmm::mat_euclidean_norm_sqr(Hk));
			gmm::mult(hk, gmm::transposed(hk), eijeijT); gmm::scale(eijeijT, -2.0); gmm::add(temp133, eijeijT, NHk);

			gmm::mult(eij, gmm::transposed(eij), eijeijT);		gmm::copy(eijeijT, m_eijeijT);		gmm::scale(m_eijeijT, -1.0); // calculate eijeijT -eijeijT
			gmm::add(Mat_I, m_eijeijT, h0_ver0); //cout << h0_ver0;// calculate h12ver1
			gmm::mult(NHk, h0_ver0, temp333); gmm::copy(temp333, h0_ver0);

			gmm::copy(Mat_I, h0_ver1); gmm::scale(h0_ver1, -1.0); gmm::add(eijeijT, h0_ver1); //calculate -I + eijeijT
			gmm::mult(gmm::transposed(eij), Eik, temp11); gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333);  // calcualte dot(eij, Eik)I + eijEikT
			gmm::add(Mat_I, m_eijeijT, temp233); gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); gmm::scale(temp233, -1.0);// calculate eij2vi = -1*eij2vj
			gmm::mult(temp333, temp233, temp133); gmm::scale(temp133, -1.0);	gmm::add(h0_ver1, temp133, temp333); gmm::copy(temp333, h0_ver1); //cout << h0_ver1;// calculate h12ver2
			gmm::mult(NHk, h0_ver1, temp333); gmm::copy(temp333, h0_ver1);

			gmm::mult(gmm::transposed(eij), Eik, temp11);  gmm::copy(Mat_I, temp133); gmm::scale(temp133, temp11(0,0)); // calculate dot(eij, Eik)I
			gmm::mult(eij, gmm::transposed(Eik), temp233); // calculate eijEikT
			gmm::add(temp133, temp233, temp333); gmm::scale(temp333, -1.0); // calculate -[dot(eij, Eik)I + eijEikT]
			gmm::add(Mat_I, m_eijeijT, temp233); gmm::scale(temp233, 1/gmm::mat_euclidean_norm(Eij)); // calculate eij2vj
			gmm::mult(temp333, temp233, h0_ver2); //cout << h0_ver2;
			gmm::mult(NHk, h0_ver2, temp333); gmm::copy(temp333, h0_ver2);

			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*0) = pen_scale*m_TPPIBG[i].v0.x + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(0,0) + ux[m_triangles[i].ver1]*h1_ver0(0,0) + ux[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(0,1) + ux[m_triangles[i].ver1]*h1_ver0(0,1) + ux[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(0,2) + ux[m_triangles[i].ver1]*h1_ver0(0,2) + ux[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*0) = pen_scale*m_TPPIBG[i].v1.x + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(0,0) + ux[m_triangles[i].ver1]*h1_ver1(0,0) + ux[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(0,1) + ux[m_triangles[i].ver1]*h1_ver1(0,1) + ux[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(0,2) + ux[m_triangles[i].ver1]*h1_ver1(0,2) + ux[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*0) = pen_scale*m_TPPIBG[i].v2.x + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(0,0) + ux[m_triangles[i].ver1]*h1_ver2(0,0) + ux[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(0,1) + ux[m_triangles[i].ver1]*h1_ver2(0,1) + ux[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+0+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(0,2) + ux[m_triangles[i].ver1]*h1_ver2(0,2) + ux[m_triangles[i].ver2]*h2_ver2(0,2));
			mat_f[i*9+0+Start_id] = pen_scale*(m_triangles[i].grad.x);
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*0) = pen_scale*m_TPPIBG[i].v0.y + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(1,0) + ux[m_triangles[i].ver1]*h1_ver0(1,0) + ux[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(1,1) + ux[m_triangles[i].ver1]*h1_ver0(1,1) + ux[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(1,2) + ux[m_triangles[i].ver1]*h1_ver0(1,2) + ux[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*0) = pen_scale*m_TPPIBG[i].v1.y +
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(1,0) + ux[m_triangles[i].ver1]*h1_ver1(1,0) + ux[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(1,1) + ux[m_triangles[i].ver1]*h1_ver1(1,1) + ux[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(1,2) + ux[m_triangles[i].ver1]*h1_ver1(1,2) + ux[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*0) = pen_scale*m_TPPIBG[i].v2.y + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(1,0) + ux[m_triangles[i].ver1]*h1_ver2(1,0) + ux[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(1,1) + ux[m_triangles[i].ver1]*h1_ver2(1,1) + ux[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+1+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(1,2) + ux[m_triangles[i].ver1]*h1_ver2(1,2) + ux[m_triangles[i].ver2]*h2_ver2(1,2));
			mat_f[i*9+1+Start_id] = pen_scale*(m_triangles[i].grad.y);
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*0) = pen_scale*m_TPPIBG[i].v0.z + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(2,0) + ux[m_triangles[i].ver1]*h1_ver0(2,0) + ux[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(2,1) + ux[m_triangles[i].ver1]*h1_ver0(2,1) + ux[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver0(2,2) + ux[m_triangles[i].ver1]*h1_ver0(2,2) + ux[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*0) = pen_scale*m_TPPIBG[i].v1.z +
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(2,0) + ux[m_triangles[i].ver1]*h1_ver1(2,0) + ux[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(2,1) + ux[m_triangles[i].ver1]*h1_ver1(2,1) + ux[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver1(2,2) + ux[m_triangles[i].ver1]*h1_ver1(2,2) + ux[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*0) = pen_scale*m_TPPIBG[i].v2.z + 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(2,0) + ux[m_triangles[i].ver1]*h1_ver2(2,0) + ux[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(2,1) + ux[m_triangles[i].ver1]*h1_ver2(2,1) + ux[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+2+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(ux[m_triangles[i].ver0]*h0_ver2(2,2) + ux[m_triangles[i].ver1]*h1_ver2(2,2) + ux[m_triangles[i].ver2]*h2_ver2(2,2));
			mat_f[i*9+2+Start_id] = pen_scale*(m_triangles[i].grad.z);
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(0,0) + uy[m_triangles[i].ver1]*h1_ver0(0,0) + uy[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*1) = pen_scale*m_TPPIBG[i].v0.x + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(0,1) + uy[m_triangles[i].ver1]*h1_ver0(0,1) + uy[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(0,2) + uy[m_triangles[i].ver1]*h1_ver0(0,2) + uy[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(0,0) + uy[m_triangles[i].ver1]*h1_ver1(0,0) + uy[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*1) = pen_scale*m_TPPIBG[i].v1.x + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(0,1) + uy[m_triangles[i].ver1]*h1_ver1(0,1) + uy[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(0,2) + uy[m_triangles[i].ver1]*h1_ver1(0,2) + uy[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(0,0) + uy[m_triangles[i].ver1]*h1_ver2(0,0) + uy[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*1) = pen_scale*m_TPPIBG[i].v2.x + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(0,1) + uy[m_triangles[i].ver1]*h1_ver2(0,1) + uy[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+3+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(0,2) + uy[m_triangles[i].ver1]*h1_ver2(0,2) + uy[m_triangles[i].ver2]*h2_ver2(0,2));
			mat_f[i*9+3+Start_id] = pen_scale*(m_triangles[i].grad1.x);
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(1,0) + uy[m_triangles[i].ver1]*h1_ver0(1,0) + uy[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*1) = pen_scale*m_TPPIBG[i].v0.y + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(1,1) + uy[m_triangles[i].ver1]*h1_ver0(1,1) + uy[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(1,2) + uy[m_triangles[i].ver1]*h1_ver0(1,2) + uy[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(1,0) + uy[m_triangles[i].ver1]*h1_ver1(1,0) + uy[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*1) = pen_scale*m_TPPIBG[i].v1.y +
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(1,1) + uy[m_triangles[i].ver1]*h1_ver1(1,1) + uy[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(1,2) + uy[m_triangles[i].ver1]*h1_ver1(1,2) + uy[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(1,0) + uy[m_triangles[i].ver1]*h1_ver2(1,0) + uy[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*1) = pen_scale*m_TPPIBG[i].v2.y + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(1,1) + uy[m_triangles[i].ver1]*h1_ver2(1,1) + uy[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+4+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(1,2) + uy[m_triangles[i].ver1]*h1_ver2(1,2) + uy[m_triangles[i].ver2]*h2_ver2(1,2));
			mat_f[i*9+4+Start_id] = pen_scale*(m_triangles[i].grad1.y);
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(2,0) + uy[m_triangles[i].ver1]*h1_ver0(2,0) + uy[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*1) = pen_scale*m_TPPIBG[i].v0.z + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(2,1) + uy[m_triangles[i].ver1]*h1_ver0(2,1) + uy[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver0+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver0(2,2) + uy[m_triangles[i].ver1]*h1_ver0(2,2) + uy[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(2,0) + uy[m_triangles[i].ver1]*h1_ver1(2,0) + uy[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*1) = pen_scale*m_TPPIBG[i].v1.z +
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(2,1) + uy[m_triangles[i].ver1]*h1_ver1(2,1) + uy[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver1+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver1(2,2) + uy[m_triangles[i].ver1]*h1_ver1(2,2) + uy[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(2,0) + uy[m_triangles[i].ver1]*h1_ver2(2,0) + uy[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*1) = pen_scale*m_TPPIBG[i].v2.z + 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(2,1) + uy[m_triangles[i].ver1]*h1_ver2(2,1) + uy[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+5+Start_id , m_triangles[i].ver2+m_vnum*2) = 
				pen_scale*(uy[m_triangles[i].ver0]*h0_ver2(2,2) + uy[m_triangles[i].ver1]*h1_ver2(2,2) + uy[m_triangles[i].ver2]*h2_ver2(2,2));
			mat_f[i*9+5+Start_id] = pen_scale*(m_triangles[i].grad1.z);
			//////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(0,0) + uz[m_triangles[i].ver1]*h1_ver0(0,0) + uz[m_triangles[i].ver2]*h2_ver0(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(0,1) + uz[m_triangles[i].ver1]*h1_ver0(0,1) + uz[m_triangles[i].ver2]*h2_ver0(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver0+m_vnum*2) = pen_scale*m_TPPIBG[i].v0.x + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(0,2) + uz[m_triangles[i].ver1]*h1_ver0(0,2) + uz[m_triangles[i].ver2]*h2_ver0(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(0,0) + uz[m_triangles[i].ver1]*h1_ver1(0,0) + uz[m_triangles[i].ver2]*h2_ver1(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(0,1) + uz[m_triangles[i].ver1]*h1_ver1(0,1) + uz[m_triangles[i].ver2]*h2_ver1(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver1+m_vnum*2) = pen_scale*m_TPPIBG[i].v1.x + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(0,2) + uz[m_triangles[i].ver1]*h1_ver1(0,2) + uz[m_triangles[i].ver2]*h2_ver1(0,2));

			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(0,0) + uz[m_triangles[i].ver1]*h1_ver2(0,0) + uz[m_triangles[i].ver2]*h2_ver2(0,0));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(0,1) + uz[m_triangles[i].ver1]*h1_ver2(0,1) + uz[m_triangles[i].ver2]*h2_ver2(0,1));
			mat_J(i*9+6+Start_id , m_triangles[i].ver2+m_vnum*2) = pen_scale*m_TPPIBG[i].v2.x + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(0,2) + uz[m_triangles[i].ver1]*h1_ver2(0,2) + uz[m_triangles[i].ver2]*h2_ver2(0,2));
			mat_f[i*9+6+Start_id] = pen_scale*(m_triangles[i].grad2.x);
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(1,0) + uz[m_triangles[i].ver1]*h1_ver0(1,0) + uz[m_triangles[i].ver2]*h2_ver0(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(1,1) + uz[m_triangles[i].ver1]*h1_ver0(1,1) + uz[m_triangles[i].ver2]*h2_ver0(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver0+m_vnum*2) = pen_scale*m_TPPIBG[i].v0.y + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(1,2) + uz[m_triangles[i].ver1]*h1_ver0(1,2) + uz[m_triangles[i].ver2]*h2_ver0(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(1,0) + uz[m_triangles[i].ver1]*h1_ver1(1,0) + uz[m_triangles[i].ver2]*h2_ver1(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(1,1) + uz[m_triangles[i].ver1]*h1_ver1(1,1) + uz[m_triangles[i].ver2]*h2_ver1(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver1+m_vnum*2) = pen_scale*m_TPPIBG[i].v1.y +
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(1,2) + uz[m_triangles[i].ver1]*h1_ver1(1,2) + uz[m_triangles[i].ver2]*h2_ver1(1,2));

			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(1,0) + uz[m_triangles[i].ver1]*h1_ver2(1,0) + uz[m_triangles[i].ver2]*h2_ver2(1,0));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(1,1) + uz[m_triangles[i].ver1]*h1_ver2(1,1) + uz[m_triangles[i].ver2]*h2_ver2(1,1));
			mat_J(i*9+7+Start_id , m_triangles[i].ver2+m_vnum*2) = pen_scale*m_TPPIBG[i].v2.y + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(1,2) + uz[m_triangles[i].ver1]*h1_ver2(1,2) + uz[m_triangles[i].ver2]*h2_ver2(1,2));
			mat_f[i*9+7+Start_id] = pen_scale*(m_triangles[i].grad2.y);
			//////////////////////////////////////////////////////////////////////////

			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(2,0) + uz[m_triangles[i].ver1]*h1_ver0(2,0) + uz[m_triangles[i].ver2]*h2_ver0(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(2,1) + uz[m_triangles[i].ver1]*h1_ver0(2,1) + uz[m_triangles[i].ver2]*h2_ver0(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver0+m_vnum*2) = pen_scale*m_TPPIBG[i].v0.z + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver0(2,2) + uz[m_triangles[i].ver1]*h1_ver0(2,2) + uz[m_triangles[i].ver2]*h2_ver0(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(2,0) + uz[m_triangles[i].ver1]*h1_ver1(2,0) + uz[m_triangles[i].ver2]*h2_ver1(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(2,1) + uz[m_triangles[i].ver1]*h1_ver1(2,1) + uz[m_triangles[i].ver2]*h2_ver1(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver1+m_vnum*2) = pen_scale*m_TPPIBG[i].v1.z +
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver1(2,2) + uz[m_triangles[i].ver1]*h1_ver1(2,2) + uz[m_triangles[i].ver2]*h2_ver1(2,2));

			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*0) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(2,0) + uz[m_triangles[i].ver1]*h1_ver2(2,0) + uz[m_triangles[i].ver2]*h2_ver2(2,0));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*1) = 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(2,1) + uz[m_triangles[i].ver1]*h1_ver2(2,1) + uz[m_triangles[i].ver2]*h2_ver2(2,1));
			mat_J(i*9+8+Start_id , m_triangles[i].ver2+m_vnum*2) = pen_scale*m_TPPIBG[i].v2.z + 
				pen_scale*(uz[m_triangles[i].ver0]*h0_ver2(2,2) + uz[m_triangles[i].ver1]*h1_ver2(2,2) + uz[m_triangles[i].ver2]*h2_ver2(2,2));
			mat_f[i*9+8+Start_id] = pen_scale*(m_triangles[i].grad2.z);
			////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////
		}
	}
	//add u-f term
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		int vertex_id = v_it.handle().idx();
		mat_J(vertex_id*3+0+m_trinum*9, vertex_id+m_vnum*0) = 1.0;
		mat_J(vertex_id*3+1+m_trinum*9, vertex_id+m_vnum*1) = 1.0;
		mat_J(vertex_id*3+2+m_trinum*9, vertex_id+m_vnum*2) = 1.0;

		mat_f[vertex_id*3+0+m_trinum*9] = (T_Mesh.point(v_it.handle())[0] - m_vertices[vertex_id].tar_x);
		mat_f[vertex_id*3+1+m_trinum*9] = (T_Mesh.point(v_it.handle())[1] - m_vertices[vertex_id].tar_y);
		mat_f[vertex_id*3+2+m_trinum*9] = (T_Mesh.point(v_it.handle())[2] - m_vertices[vertex_id].tar_z);
	}
}

void TriangularMesh::GradientTesting(double h = 0.0005, bool TestTV = false, int choice = 0)
{
	MyMesh T_Mesh = this->m_ObjTriMesh;
	RowSparseMatrix J, J_, J_estimate;
	DenseMatrix fx, fx_e, fx_me, fx_diff;
	int ntri = T_Mesh.n_faces();
	vector<VECTOR3D> px, py, pz; px.resize(ntri); py.resize(ntri); pz.resize(ntri);
	vector<VECTOR3D> lambda_x, lambda_y, lambda_z; lambda_x.resize(ntri); lambda_y.resize(ntri); lambda_z.resize(ntri);
	for(int i=0;i<ntri;i++)
	{
		px[i].x=0; lambda_x[i].x=0;		px[i].y=0; lambda_x[i].y=0;		px[i].z=0; lambda_x[i].z=0;
		py[i].x=0; lambda_y[i].x=0;		py[i].y=0; lambda_y[i].y=0;		py[i].z=0; lambda_y[i].z=0;
		pz[i].x=0; lambda_z[i].x=0;		pz[i].y=0; lambda_z[i].y=0;		pz[i].z=0; lambda_z[i].z=0;
	}
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		m_vertices[v_it.handle().idx()].light_x = m_vertices[v_it.handle().idx()].x;
		m_vertices[v_it.handle().idx()].light_y = m_vertices[v_it.handle().idx()].y;
		m_vertices[v_it.handle().idx()].light_z = m_vertices[v_it.handle().idx()].z;
	}
	double fidParam = 0,  pld_eta = 0,  pcd_eta = 0,  fcd_eta = 0,  pnd_eta = 0, fnd_eta = 0, varsigma = 100,  pc_eta = 0, penParam = 0.01,  lapParam = 0;
	bool UseTVU = false, UseTVNorm = false;

	if (TestTV) {
		cout << "Testing the tv gradient matrix: " << choice << endl;
		switch (choice)
		{
		case 0:
			fidParam = pld_eta = pcd_eta = fcd_eta = lapParam = 1.0; UseTVU = true;
			break;
		case 1:
			fidParam = pld_eta = pcd_eta = fcd_eta = lapParam = 1.0; UseTVNorm = true;
			break;
		default:
			break;
		}
	} else {
		cout << "Testing the normal gradient matrix: " << choice << endl;
		switch (choice)
		{
		case 0:
			pcd_eta = 1.0;
			break;
		case 1:
			fcd_eta = 1.0; 
			break;
		case 2:
			pcd_eta = 1.0; fcd_eta = 1.0; pld_eta = 1.0; lapParam = 1.0;
			break;
		case 3:
			lapParam = 1.0; fcd_eta = 1.0; pcd_eta = 1.0; pld_eta = 1.0; fidParam = 1.0;
			break;
		default:
			break;
		}

	}
	this->TV_JacobianMatrix_Construction(T_Mesh, J, fx, fidParam, pld_eta, pcd_eta, fcd_eta, pnd_eta, fnd_eta, varsigma, pc_eta, penParam, lapParam, px, py, pz, lambda_x, lambda_y, lambda_z, UseTVU, UseTVNorm, false);
	char buffer[255]; fstream nof;
	sprintf(buffer, "Results\\mat_J_%s%d.txt", TestTV?"TV":"_",choice);
	nof.open(buffer, ios::out);
	if (!nof) {
		cout << "Can not open txt file to save mat_J data..." << endl;
		return;
	}
	int r = 0;
	for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(J);
		rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(J); ++ rIt, ++ r) {
			nof << r << ": ";
			gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
			gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
			for (; it != ite; ++ it) {
				nof << it.index() << "," << J(r, it.index()) << "; ";
			}
			nof << endl;
	}
	nof.close();
	J_estimate.resize(J.nrows(), J.ncols());	gmm::clear(J_estimate); fx_diff.resize(fx.nrows(), fx.ncols()); 
	for (int idx = 0; idx < m_vnum; ++ idx) {
		for (int j = 0; j < 3; ++ j) {
			int dj = idx+m_vnum*j;
			OpenMesh::Vec3f vi_ori = T_Mesh.point(MyMesh::VertexHandle(idx));
			OpenMesh::Vec3f veci = vi_ori;

			veci[j] += h;
			T_Mesh.set_point(MyMesh::VertexHandle(idx), veci);
			this->TV_JacobianMatrix_Construction(T_Mesh, J_, fx_e, fidParam, pld_eta, pcd_eta, fcd_eta, pnd_eta, fnd_eta, varsigma, pc_eta,penParam, lapParam, px, py, pz, lambda_x, lambda_y, lambda_z, UseTVU, UseTVNorm, false);

			veci[j] -= 2*h;
			T_Mesh.set_point(MyMesh::VertexHandle(idx), veci);
			this->TV_JacobianMatrix_Construction(T_Mesh, J_, fx_me, fidParam, pld_eta, pcd_eta, fcd_eta, pnd_eta, fnd_eta, varsigma, pc_eta,penParam, lapParam, px, py, pz, lambda_x, lambda_y, lambda_z, UseTVU, UseTVNorm, false);
			gmm::scale(fx_me, -1.0);  gmm::scale(fx_diff, 0.0);
			gmm::add(fx_e, fx_me, fx_diff);
			gmm::scale(fx_diff,1./(2.*h));
			for( int fi=0;fi<J_estimate.nrows();fi++){
				J_estimate(fi,dj) = fx_diff(fi,0);
			}
			T_Mesh.set_point(MyMesh::VertexHandle(idx), vi_ori);
		}
	}
	double norm_J = gmm::mat_euclidean_norm(J);				cout << "The norm of prop gradient matrix is: " << norm_J << endl;
	double norm_Je = gmm::mat_euclidean_norm(J_estimate);	cout << "The norm of real gradient matrix is: " << norm_Je << endl;
	sprintf(buffer, "Results\\mat_J_estimate_%s%d.txt", TestTV?"TV":"_", choice);
	fstream enof(buffer,std::ios::out);
	if (!enof) {
		cout << "Can not open txt file to save mat_J data..." << endl;
		return;
	}
	r = 0;
	for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(J_estimate);
		rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(J_estimate); ++ rIt, ++ r) {
			enof << r << ": ";
			gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
			gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
			for (; it != ite; ++ it) {
				enof << it.index() << "," << J_estimate(r, it.index()) << "; ";
			}
			enof << endl;
	}
	enof.close();

	gmm::scale(J_estimate, -1.0);
	gmm::add(J, J_estimate, J_);
	double norm_J_ = gmm::mat_euclidean_norm(J_);			cout << "The norm of the difference matrix is: " << norm_J_ << endl;
	sprintf(buffer, "Results\\mat_J_difference_%s%d.txt", TestTV?"TV":"_", choice);
	enof.open(buffer,std::ios::out);
	if (!enof) {
		cout << "Can not open txt file to save mat_J_difference data..." << endl;
		return;
	}
	r = 0;
	for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(J_);
		rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(J_); ++ rIt, ++ r) {
			enof << r << ": ";
			gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
			gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
			for (; it != ite; ++ it) {
				enof << it.index() << "," << J_(r, it.index()) << "; ";
			}
			enof << endl;
	}
	enof.close();
}


void TriangularMesh::ALM_MeshSmooth()
{
/*  Solves vectorial TV L2 denoising problem using ALM.

    The problem:
	    min_u regParam*\sum|\nabla u|s_tau + 1/2 \|u-f\|_2.

    Here u is a multi-channel image.

    Parameters:
	   penParam : the penalty parameter used in ALM;
	   innerL   : the number of inner iteration, usually L=1;
	   outTole  : the tolerance for outer iteration;
 */
	MyMesh T_Mesh = this->m_ObjTriMesh;
	
    unsigned long nver,ntri,i;
	double *ur,*ug,*ub,*b;
	double *urold,*ugold,*ubold;
	VECTOR3D *pr,*pg,*pb,*lambdar,*lambdag,*lambdab,wr,wg,wb,w;
	int itol,itmax,iter;
    double tol,err;
	
	nver = m_vnum;
	ntri = m_trinum;
	itmax = 300;
	itol = 1;//1,2,3, or 4.
	tol = 1.0e-10;
	
	ur = new double[nver];
	ug = new double[nver];
	ub = new double[nver];
	urold = new double[nver];
	ugold = new double[nver];
	ubold = new double[nver];
	b = new double[nver];
	pr = new VECTOR3D[ntri];
	pg = new VECTOR3D[ntri];
	pb = new VECTOR3D[ntri];
	lambdar = new VECTOR3D[ntri];
	lambdag = new VECTOR3D[ntri];
	lambdab = new VECTOR3D[ntri];

    unsigned int innerL,l;
	double regParam, fidParam, penParam, outTole, stoppingCond;
/*   Parameters setting  */
	innerL = 1;
	fidParam = 100;
	regParam = 1;
	penParam = 0.2;
	outTole = 1.0e-10;

	numc::RowMatSym<double> alphaMinusrLap(nver,nver);
	ONERING or;
	double diag;
	unsigned long j;
	for(i=0;i<nver;i++)
	{
		or = m_ONERINGS[i];
		diag = 0;
		for(j=0;j<or.size();j++)
		{
			alphaMinusrLap(i,or[j].iver) = -penParam*or[j].co;
			diag-= or[j].co;
		}
		alphaMinusrLap(i,i) = fidParam*m_vertices[i].BCDArea - penParam*diag;
	}
    numc::SparseSolver solver;
    solver.getMatA() = alphaMinusrLap;
    solver.getMatA().mMtype = numc::CSRMatrix<double>::RealSymmIndef;
    solver.init();

	vector<double> vec_n_energy; vec_n_energy.resize(m_vnum); char buffer[255]; 
	for (int i = 0; i < vec_n_energy.size(); ++ i) {
		vec_n_energy[i] = 0.0;
	}
	double pld_eta = 0;
	vector<double>ux, uy, uz; ux.resize(nver); uy.resize(nver); uz.resize(nver);
	int outL = -1;

/*   Initialization   */
	for(i=0;i<nver;i++)
	{	
		ur[i] = T_Mesh.point(MyMesh::VertexHandle(i))[0]; ug[i] = T_Mesh.point(MyMesh::VertexHandle(i))[1]; ub[i] = T_Mesh.point(MyMesh::VertexHandle(i))[2];
	}
	for(i=0;i<ntri;i++)
	{
		pr[i].x=0; lambdar[i].x=0;
		pr[i].y=0; lambdar[i].y=0;
		pr[i].z=0; lambdar[i].z=0;

		pg[i].x=0; lambdag[i].x=0;
		pg[i].y=0; lambdag[i].y=0;
		pg[i].z=0; lambdag[i].z=0;

		pb[i].x=0; lambdab[i].x=0;
		pb[i].y=0; lambdab[i].y=0;
		pb[i].z=0; lambdab[i].z=0;
	}
/*   Iteration        */
	do
	{
		outL ++;
		for(i=0;i<nver;i++)
		{
			urold[i] = ur[i];
			ugold[i] = ug[i];
			ubold[i] = ub[i];
		}
		/*  Inner iteration             */
		for(l=0;l<innerL;l++)
		{
			for(i=0;i<nver;i++){	
				ux[i] = ur[i]; uy[i] = ug[i]; uz[i] = ub[i];
			}
			OpenMesh::Vec3f CurV;
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				CurV[0] = ux[v_it.handle().idx()]; CurV[1] = uy[v_it.handle().idx()]; CurV[2] = uz[v_it.handle().idx()];
				T_Mesh.set_point(v_it.handle(), CurV);
			}
			double v_energy = CalculateVEnergy(T_Mesh);
			double n_energy = CalculateNEnergy(T_Mesh, vec_n_energy);
#ifdef _DEBUG
			sprintf(buffer, "Results\\NEnergy_%d_%d_%.0f_%.0f_%.3f.txt", outL, iter, fidParam, pld_eta, n_energy);
			fstream nof(buffer,std::ios::out);
			if (!nof) {
				cout << "Can not open txt file to save N energy data..." << endl;
				return;
			}
			for (int i = 0; i < vec_n_energy.size(); ++ i) {
				nof << vec_n_energy[i] << endl;
			}
			nof.close();
#endif
			double tv_energy = CalculateVTVEnergy(T_Mesh, ux, uy, uz);
			cout << "The current energy is V:" << 0.5*fidParam*v_energy << ", N:" << pld_eta*n_energy << ", TV: " << penParam*tv_energy << "; Sum: " 
				<< 0.5*fidParam*v_energy + pld_eta * n_energy + penParam * tv_energy << endl;
			OpenMesh::IO::Options write_options;
			sprintf(buffer, "Results\\ALMMeshResult_%d_%.0f_%.0f_%.3f.off", outL, fidParam, pld_eta, penParam);
			if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
				std::cerr << "Cannot write mesh to file " << buffer << std::endl;
			}
		/*  Solve the u-sub problem  (solve for ur,ug,ub).   */

		   for(i=0;i<ntri;i++)
		   {
		    	pr[i].x = lambdar[i].x+penParam*pr[i].x;
			    pr[i].y = lambdar[i].y+penParam*pr[i].y;
			    pr[i].z = lambdar[i].z+penParam*pr[i].z;
		   }
		   GetDivergence(ntri,pr);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].r-m_vertices[i].divergence);
		   solver.solve(b,ur);
//         ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,ur,fidParam,penParam,itol,tol,itmax,&iter,&err);

		   for(i=0;i<ntri;i++)
		   {
		    	pg[i].x = lambdag[i].x+penParam*pg[i].x;
			    pg[i].y = lambdag[i].y+penParam*pg[i].y;
			    pg[i].z = lambdag[i].z+penParam*pg[i].z;
		   }
		   GetDivergence(ntri,pg);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].g-m_vertices[i].divergence);
		   solver.solve(b,ug);
//         ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,ug,fidParam,penParam,itol,tol,itmax,&iter,&err);

		   for(i=0;i<ntri;i++)
		   {
		    	pb[i].x = lambdab[i].x+penParam*pb[i].x;
			    pb[i].y = lambdab[i].y+penParam*pb[i].y;
			    pb[i].z = lambdab[i].z+penParam*pb[i].z;
		   }
		   GetDivergence(ntri,pb);
		   for(i=0;i<nver;i++)
		    	b[i] = m_vertices[i].BCDArea*(fidParam*m_tempRGB[i].b-m_vertices[i].divergence);
		   solver.solve(b,ub);
//         ALM_vTVL2Denoising_usub_nr_linbcg(nver,b,ub,fidParam,penParam,itol,tol,itmax,&iter,&err);

        /*  Solve the p-sub problem  (solve for pr,pg,pb).   */
           
		   GetPPIGradients(nver,ur);
		   GetPPIGradients1(nver,ug);
		   GetPPIGradients2(nver,ub);

		   for(i=0;i<ntri;i++)
		   {
			    wr.x = m_triangles[i].grad.x - lambdar[i].x/penParam;
		    	wr.y = m_triangles[i].grad.y - lambdar[i].y/penParam;
			    wr.z = m_triangles[i].grad.z - lambdar[i].z/penParam;

			    wg.x = m_triangles[i].grad1.x - lambdag[i].x/penParam;
		    	wg.y = m_triangles[i].grad1.y - lambdag[i].y/penParam;
			    wg.z = m_triangles[i].grad1.z - lambdag[i].z/penParam;

			    wb.x = m_triangles[i].grad2.x - lambdab[i].x/penParam;
		    	wb.y = m_triangles[i].grad2.y - lambdab[i].y/penParam;
			    wb.z = m_triangles[i].grad2.z - lambdab[i].z/penParam;

				if(DotProduct(wr,wr)+DotProduct(wg,wg)+DotProduct(wb,wb) <= POWER(regParam/penParam))
				{
				     pr[i].x = 0; pg[i].x = 0; pb[i].x = 0;
				     pr[i].y = 0; pg[i].y = 0; pb[i].y = 0;
				     pr[i].z = 0; pg[i].z = 0; pb[i].z = 0;
				}
			    else
				{
					 double tempNorm;
					 tempNorm = sqrt(DotProduct(wr,wr)+DotProduct(wg,wg)+DotProduct(wb,wb));

					 pr[i].x = (1-regParam/penParam/tempNorm)*wr.x; pg[i].x = (1-1/penParam/tempNorm)*wg.x; pb[i].x = (1-1/penParam/tempNorm)*wb.x;
				     pr[i].y = (1-regParam/penParam/tempNorm)*wr.y; pg[i].y = (1-1/penParam/tempNorm)*wg.y; pb[i].y = (1-1/penParam/tempNorm)*wb.y;
				     pr[i].z = (1-regParam/penParam/tempNorm)*wr.z; pg[i].z = (1-1/penParam/tempNorm)*wg.z; pb[i].z = (1-1/penParam/tempNorm)*wb.z;
				}
		   } 
		}
        /*   Update Lagrange multipliers  (lambdar,lambdag,lambdab).      */
//		GetPPIGradients(nver,ur);
//		GetPPIGradients1(nver,ug);
//		GetPPIGradients2(nver,ub);

		for(i=0;i<ntri;i++)
		{
            lambdar[i].x+= penParam*(pr[i].x-m_triangles[i].grad.x);
			lambdar[i].y+= penParam*(pr[i].y-m_triangles[i].grad.y);
			lambdar[i].z+= penParam*(pr[i].z-m_triangles[i].grad.z);

            lambdag[i].x+= penParam*(pg[i].x-m_triangles[i].grad1.x);
			lambdag[i].y+= penParam*(pg[i].y-m_triangles[i].grad1.y);
			lambdag[i].z+= penParam*(pg[i].z-m_triangles[i].grad1.z);

            lambdab[i].x+= penParam*(pb[i].x-m_triangles[i].grad2.x);
			lambdab[i].y+= penParam*(pb[i].y-m_triangles[i].grad2.y);
			lambdab[i].z+= penParam*(pb[i].z-m_triangles[i].grad2.z);
		}
		/*   Compute the stopping condition     */
        stoppingCond = 0;
/*		VECTOR3D tempr,tempg,tempb;
		for(i=0;i<ntri;i++)
		{
			tempr.x = pr[i].x-m_triangles[i].grad.x;
			tempr.y = pr[i].y-m_triangles[i].grad.y;
			tempr.z = pr[i].z-m_triangles[i].grad.z;
			
			tempg.x = pg[i].x-m_triangles[i].grad1.x;
			tempg.y = pg[i].y-m_triangles[i].grad1.y;
			tempg.z = pg[i].z-m_triangles[i].grad1.z;
			
			tempb.x = pb[i].x-m_triangles[i].grad2.x;
			tempb.y = pb[i].y-m_triangles[i].grad2.y;
			tempb.z = pb[i].z-m_triangles[i].grad2.z;

			stoppingCond+= m_trianglesArea[i]*(POWER(tempr.x)+POWER(tempr.y)+POWER(tempr.z)
				+POWER(tempg.x)+POWER(tempg.y)+POWER(tempg.z)
				+POWER(tempb.x)+POWER(tempb.y)+POWER(tempb.z));
		}
*/
		for(i=0;i<nver;i++)
			stoppingCond+= (POWER(ur[i]-urold[i])+POWER(ug[i]-ugold[i])+POWER(ub[i]-ubold[i])) * m_vertices[i].BCDArea;

	}
	while(stoppingCond>outTole);

	for(i=0;i<nver;i++)
	{
		m_vertices[i].r = ur[i]; m_vertices[i].g = ug[i]; m_vertices[i].b = ub[i];
	}


	delete ur; delete ug; delete ub;
	delete urold; delete ugold; delete ubold;
	delete b;
	delete pr; delete pg; delete pb;
	delete lambdar; delete lambdag; delete lambdab;
}

//void TriangularMesh::ALM_TVNorm_MeshRefinement(double fidParam, double pld_eta, double penParam, double regParam = 1.0, bool IsComplexVersion = false)
//{
//	unsigned long nver = m_vnum, ntri = m_trinum; int itmax = 300, itol = 1;//1,2,3, or 4.
//	double tol = 1.0e-15;
//	vector<double>ux, uy, uz; ux.resize(nver); uy.resize(nver); uz.resize(nver);
//	vector<double>uxold, uyold, uzold; uxold.resize(nver); uyold.resize(nver); uzold.resize(nver);
//	vector<double>nux, nuy, nuz; nux.resize(nver); nuy.resize(nver); nuz.resize(nver);
//	vector<VECTOR3D> px, py, pz; px.resize(ntri); py.resize(ntri); pz.resize(ntri);
//	vector<VECTOR3D> lambda_x, lambda_y, lambda_z; lambda_x.resize(ntri); lambda_y.resize(ntri); lambda_z.resize(ntri);
//	/*   Parameters setting  */
//	unsigned int innerL = 1; int outL = -1;
//	double outTole = 1.0e-10, stoppingCond;// fidParam = 1000, pld_eta = 0, penParam = 0.01, regParam = 1,
//	double var = 2.0, epsilon1 = 1.0e-5, epsilon2 = 1.0e-5, minEdgeLength = 5.0e-3;
//	numc::RowMat<double> alphaMinusrLap(nver*3,nver*3);
//	numc::RowMat<double> NormVMatrix(nver*3,nver*3);
//	double *b = new double[nver*3]; double *vx = new double[nver*3]; vector<double> vec_b; vec_b.resize(nver*3);
//	numc::SparseSolver solver;	
//	MyMesh T_Mesh = this->m_ObjTriMesh; OpenMesh::Vec3f CurV; OpenMesh::IO::Options write_options;
//	cout << "Start the mesh refinement process: " << endl;
//	for(int i=0;i<nver;i++) {	
//		ux[i] = m_vertices[i].x; uy[i] = m_vertices[i].y; uz[i] = m_vertices[i].z;
//		//ux[i] = 0; uy[i] = 0; uz[i] = 0;
//	}
//	for(int i=0;i<ntri;i++)
//	{
//		px[i].x=0; lambda_x[i].x=0;		px[i].y=0; lambda_x[i].y=0;		px[i].z=0; lambda_x[i].z=0;
//		py[i].x=0; lambda_y[i].x=0;		py[i].y=0; lambda_y[i].y=0;		py[i].z=0; lambda_y[i].z=0;
//		pz[i].x=0; lambda_z[i].x=0;		pz[i].y=0; lambda_z[i].y=0;		pz[i].z=0; lambda_z[i].z=0;
//	}
//	vector<vector<double>> energy_result; energy_result.clear();
//	vector<double> vec_n_energy; vec_n_energy.resize(m_vnum); char buffer[255]; 
//	for (int i = 0; i < vec_n_energy.size(); ++ i) {
//		vec_n_energy[i] = 0.0;
//	}
//	do { // start the iteration
//		outL ++;
//		for(int i=0;i<m_vnum;i++) {
//			uxold[i] = ux[i];	uyold[i] = uy[i];	uzold[i] = uz[i];	//b[i] = 0;  b[i+m_vnum] = 0; b[i+m_vnum*2] = 0;	vec_b[i] = 0;  vec_b[i+m_vnum] = 0; vec_b[i+m_vnum*2] = 0;
//		}
//		// inner iteration
//		for (int iter = 0; iter < innerL; ++ iter) {
//			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//				CurV[0] = ux[v_it.handle().idx()]; CurV[1] = uy[v_it.handle().idx()]; CurV[2] = uz[v_it.handle().idx()];
//				T_Mesh.set_point(v_it.handle(), CurV);
//			}
//			double v_energy = CalculateVEnergy(T_Mesh);
//			double n_energy = CalculateNEnergy(T_Mesh, vec_n_energy);
//			double tv_energy = CalculateTVEnergy(ux, uy, uz);
//			vector<double> t_energy; t_energy.clear(); t_energy.push_back(0.5*fidParam*v_energy); t_energy.push_back(pld_eta*n_energy);
//			t_energy.push_back(regParam*tv_energy); t_energy.push_back(0.5*fidParam * v_energy + pld_eta * n_energy + regParam * tv_energy);
//			energy_result.push_back(t_energy);
//			cout << "The current energy is V:" << 0.5*fidParam*v_energy << ", N:" << pld_eta*n_energy << ", TV: " << regParam*tv_energy << "; Sum: " 
//				<< 0.5*fidParam * v_energy + pld_eta * n_energy + regParam * tv_energy << endl;
//			OpenMesh::IO::Options write_options;
//			write_options.set(OpenMesh::IO::Options::VertexNormal); 
//			sprintf(buffer, "Results\\TVNorm_MeshResult_%d_%.0f_%.0f_%.3f_%s.off", outL, fidParam, pld_eta, penParam, IsComplexVersion?"Complex":"Simple");
//			if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
//				std::cerr << "Cannot write mesh to file " << buffer << std::endl;
//			}
//
//			// solve the u-sub problem, using Levenberg–Marquardt Algorithm
//			double sigma, eOld, eNew;
//			int dimension = m_vnum * 3;
//			std::vector<double> m_b(dimension, 0);
//			std::vector<double> m_x(dimension, 0);
//			std::vector<double> tau(dimension, 0);
//			std::vector<OpenMesh::Vec3f> originalPhi(m_vnum);
//			RowSparseMatrix mat_J;	std::vector<double> m_f;
//
//			bool OnlyNormTerm = true;
//			TVNorm_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, fidParam, pld_eta, penParam, px, py, pz, lambda_x, lambda_y, lambda_z, IsComplexVersion);
//
//			eOld = gmm::vect_norm2_sqr(m_f);
//			RowSparseMatrix mat_JTJ(dimension, dimension);
//			RowSparseMatrix mat_A(dimension, dimension);
//			gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
//			for (int i = 0; i < dimension; i++) {
//				tau[i] = mat_JTJ(i, i);
//			}
//			gmm::scale(tau, 1.e-5);
//			gmm::mult(gmm::transposed(mat_J), m_f, m_b);
//			gmm::scale(m_b, -1.0);
//			double normG = gmm::vect_norminf(m_b);
//
//			numc::SparseSolver solver;	double *b = new double[dimension];	double *vx = new double[dimension];
//			int iteration = 0; int maxStep = 100;  //100
//			while (iteration++ < maxStep) {
//				cout << ".";
//				if (normG < epsilon1) {
//					break;
//				}
//				gmm::copy(mat_JTJ, mat_A);
//				for (int i = 0; i < dimension; i++) {
//					mat_A(i, i) += tau[i];
//				}
//				// solve mat_A*m_x = m_b
//				numc::RowMat<double> RM_A(dimension, dimension); int r = 0;
//				for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(mat_A);
//					rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(mat_A); ++ rIt, ++ r) {
//						gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
//						gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
//						for (; it != ite; ++ it) {
//							RM_A(r, it.index()) = mat_A(r, it.index());
//						}
//						b[r] = m_b[r];
//				}
//				solver.getMatA() = RM_A;
//				solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
//				solver.init();
//				// solve the nonlinear equation to calculate the new vertex position in vx;
//				solver.solve(b, vx);
//				solver.clear();
//				for (int r = 0; r < dimension; ++ r) {
//					m_x[r] = vx[r];
//				}
//
//				double normH = 0.0, normX = 0.0, l0h = 0;
//				normH = gmm::vect_norm2(m_x);
//				for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//					int idx = v_it.handle().idx();
//					originalPhi[v_it.handle().idx()] = T_Mesh.point(v_it.handle());
//					normX += originalPhi[v_it.handle().idx()].sqrnorm();
//					for (int i = 0; i < 3; ++i) {
//						l0h += 0.5*m_x[i*m_vnum + idx]*(tau[i*m_vnum + idx]*m_x[i*m_vnum + idx] + m_b[i*m_vnum + idx]);
//					}
//				}
//				if (normH < epsilon2*(sqrt(normX) + epsilon2)) {
//					break;
//				} else {
//					//change the coordinate value of each vertex...
//					OpenMesh::Vec3f CurV;
//					for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//						int vIndex = v_it.handle().idx();
//						T_Mesh.point(v_it.handle()) += OpenMesh::Vec3f(m_x[vIndex], m_x[m_vnum + vIndex], m_x[2*m_vnum + vIndex]);
//					}
//					if (iteration%1 == 0) {
//						write_options.set(OpenMesh::IO::Options::VertexNormal); 
//						sprintf(buffer, "Results\\LMResults\\TVNorm_LMMeshResult_%d_%d_%s.off", outL, iteration, IsComplexVersion?"Complex":"Simple");
//						if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
//							std::cerr << "Cannot write mesh to file " << buffer << std::endl;
//						}
//					}
//					TVNorm_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, fidParam, pld_eta, penParam, px, py, pz, lambda_x, lambda_y, lambda_z, IsComplexVersion);
//					eNew = gmm::vect_norm2_sqr(m_f);
//					sigma = (eOld - eNew)/l0h;
//					if (sigma >= 0) {
//						gmm::scale(tau, std::max(1./3., 1. - pow(2.*sigma - 1., 3.)));
//						var = 2.0;
//						eOld = eNew;
//						gmm::mult(gmm::transposed(mat_J), m_f, m_b);
//						gmm::scale(m_b, -1.0);
//						normG = gmm::vect_norminf(m_b);
//						gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
//					} else {
//						//cout << " Sigma = " << sigma ;
//						gmm::scale(tau, var);
//						var *= 2.0;
//						for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//							T_Mesh.point(v_it) = originalPhi[v_it.handle().idx()];
//						}
//					}
//				}
//			}
//			write_options.set(OpenMesh::IO::Options::VertexNormal); 
//			sprintf(buffer, "Results\\LMResults\\TVNorm_LMMeshResult_Final_%d_%s.off", outL, IsComplexVersion?"Complex":"Simple");
//			if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
//				std::cerr << "Cannot write mesh to file " << buffer << std::endl;
//			}
//			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//				ux[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[0];
//				uy[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[1];
//				uz[v_it.handle().idx()] = T_Mesh.point(v_it.handle())[2];
//			}
//			epsilon1 *= 0.001; epsilon2 *= 0.001;
//
//			/*  Solve the p-sub problem  (solve for px,py,pz).   */
//			T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
//			T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
//			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
//				nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle())[0];
//				nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle())[1];
//				nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle())[2];
//			}
//			GetPPIGradients (m_vnum,nux);
//			GetPPIGradients1(m_vnum,nuy);
//			GetPPIGradients2(m_vnum,nuz);
//			VECTOR3D wx, wy, wz;
//			double pdiff = 0.0, avg_wxy = 0.0; int i1 = 0, i2 = 0;
//			for(int i=0;i<ntri;i++)
//			{
//				wx.x = m_triangles[i].grad.x - lambda_x[i].x/penParam;
//				wx.y = m_triangles[i].grad.y - lambda_x[i].y/penParam;
//				wx.z = m_triangles[i].grad.z - lambda_x[i].z/penParam;
//
//				wy.x = m_triangles[i].grad1.x - lambda_y[i].x/penParam;
//				wy.y = m_triangles[i].grad1.y - lambda_y[i].y/penParam;
//				wy.z = m_triangles[i].grad1.z - lambda_y[i].z/penParam;
//
//				wz.x = m_triangles[i].grad2.x - lambda_z[i].x/penParam;
//				wz.y = m_triangles[i].grad2.y - lambda_z[i].y/penParam;
//				wz.z = m_triangles[i].grad2.z - lambda_z[i].z/penParam;
//				avg_wxy += sqrt(NormSquare(wx)+NormSquare(wy)+NormSquare(wz));
//
//				if(DotProduct(wx,wx)+DotProduct(wy,wy)+DotProduct(wz,wz) <= POWER(regParam/penParam)) {
//					px[i].x = 0; py[i].x = 0; pz[i].x = 0;
//					px[i].y = 0; py[i].y = 0; pz[i].y = 0;
//					px[i].z = 0; py[i].z = 0; pz[i].z = 0;
//					i1 ++;
//				} else {
//					double tempNorm; i2++;
//					tempNorm = sqrt(DotProduct(wx,wx)+DotProduct(wy,wy)+DotProduct(wz,wz));
//					px[i].x = (1-regParam/penParam/tempNorm)*wx.x; py[i].x = (1-1/penParam/tempNorm)*wy.x; pz[i].x = (1-1/penParam/tempNorm)*wz.x;
//					px[i].y = (1-regParam/penParam/tempNorm)*wx.y; py[i].y = (1-1/penParam/tempNorm)*wy.y; pz[i].y = (1-1/penParam/tempNorm)*wz.y;
//					px[i].z = (1-regParam/penParam/tempNorm)*wx.z; py[i].z = (1-1/penParam/tempNorm)*wy.z; pz[i].z = (1-1/penParam/tempNorm)*wz.z;
//				}
//				pdiff += Norm(px[i]) + Norm(py[i]) + Norm(pz[i]);
//			}
//			cout << endl << "Solve p sub-problem: Pdiff is: " << pdiff << "  " << i1 << "/" << i2 << " under " << avg_wxy/ntri << endl;
//		}
//		/*   Update Lagrange multipliers  (lambda_x,lambda_y,lambda_z).      */
//
//		for(int i=0;i<ntri;i++)
//		{
//			lambda_x[i].x+= penParam*(px[i].x-m_triangles[i].grad.x);
//			lambda_x[i].y+= penParam*(px[i].y-m_triangles[i].grad.y);
//			lambda_x[i].z+= penParam*(px[i].z-m_triangles[i].grad.z);
//
//			lambda_y[i].x+= penParam*(py[i].x-m_triangles[i].grad1.x);
//			lambda_y[i].y+= penParam*(py[i].y-m_triangles[i].grad1.y);
//			lambda_y[i].z+= penParam*(py[i].z-m_triangles[i].grad1.z);
//
//			lambda_z[i].x+= penParam*(pz[i].x-m_triangles[i].grad2.x);
//			lambda_z[i].y+= penParam*(pz[i].y-m_triangles[i].grad2.y);
//			lambda_z[i].z+= penParam*(pz[i].z-m_triangles[i].grad2.z);
//		}
//		/*   Compute the stopping condition     */
//		stoppingCond = 0;
//		for(int i=0;i<nver;i++){
//			stoppingCond+= (POWER(ux[i]-uxold[i])+POWER(uy[i]-uyold[i])+POWER(uz[i]-uzold[i])) * m_vertices[i].BCDArea;
//		}
//		cout << "The error is: " << stoppingCond << "; the threshold is: " << outTole << endl;
//	} while (stoppingCond>outTole);
//	for(int i=0;i<nver;i++){
//		m_vertices[i].x = ux[i]; m_vertices[i].y = uy[i]; m_vertices[i].z = uz[i];
//		//m_vertices[i].ref_x = ux[i]; m_vertices[i].ref_y = uy[i]; m_vertices[i].ref_z = uz[i];
//	}
//	T_Mesh = this->m_ObjTriMesh;
//	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//		CurV[0] = m_vertices[v_it.handle().idx()].x; CurV[1] = m_vertices[v_it.handle().idx()].y; CurV[2] = m_vertices[v_it.handle().idx()].z;
//		T_Mesh.set_point(v_it.handle(), CurV);
//	}
//	write_options.set(OpenMesh::IO::Options::VertexNormal); 
//	if (IsComplexVersion) {
//		if ( !OpenMesh::IO::write_mesh(T_Mesh, "Results\\ALM_TVNorm_FinalMesh_Complex.off", write_options) ) {
//			std::cerr << "Cannot write mesh to file " << "ALM_TVNorm_FinalMesh_Complex.off" << std::endl;
//		}
//	} else {
//		if ( !OpenMesh::IO::write_mesh(T_Mesh, "Results\\ALM_TVNorm_FinalMesh_Simple.off", write_options) ) {
//			std::cerr << "Cannot write mesh to file " << "ALM_TVNorm_FinalMesh_Simple.off" << std::endl;
//		}
//	}
//	fstream nof("Results\\TVNorm_EnergyResult.txt",std::ios::out);
//	if (!nof) {
//		cout << "Can not open TVNorm_EnergyResult.txt to save data..." << endl;
//		return;
//	}
//	for (int i = 0; i < energy_result.size(); ++ i) {
//		for (int j = 0; j < energy_result[i].size(); ++ j) {
//			nof << setiosflags(ios::fixed)<<setprecision(6)<<setw(10)<<setiosflags(ios::left) << energy_result[i][j] << "  ";
//		}
//		nof << endl;
//	}
//	nof.close();
//	delete b; delete vx;
//	return;
//}
//
//
//void TriangularMesh::TVNorm_JacobianMatrix_Construction(MyMesh& T_Mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, double fidParam, double pld_eta, double penParam, vector<VECTOR3D> &px, vector<VECTOR3D> &py, vector<VECTOR3D> &pz, vector<VECTOR3D> &lambda_x, vector<VECTOR3D> &lambda_y, vector<VECTOR3D> &lambda_z, bool IsComplexVersion = false)
//{
//	// construct the jacobian matrix according to the new model
//	mat_J.resize(m_vnum*6, m_vnum*3);	mat_f.resize(m_vnum*6);
//	if (IsComplexVersion) {
//		mat_J.resize(m_vnum*6+m_trinum*9, m_vnum*3);	mat_f.resize(m_vnum*6+m_trinum*9);
//	}
//	gmm::clear(mat_J);							gmm::clear(mat_f);
//	double fid_scale = sqrt(fidParam*0.5), beta_scale = sqrt(pld_eta*0.5), pen_scale = sqrt(penParam*0.5);
//
//	OpenMesh::Vec3f pointA , pointB , pointC; int ida, idb, idc;
//	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
//		const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
//		// first for all neighbor triangle, store its vertex in counter-clock wise, A(current vi), B, C
//		vector< OpenMesh::Vec3f >	NeiVertexList;		NeiVertexList.clear();
//		vector< int >				NeiVertexIDList;	NeiVertexIDList.clear();
//		vector< double >			NeiFaceAreaList;	NeiFaceAreaList.clear();
//		for (MyMesh::ConstVertexVertexIter vv_it = T_Mesh.cvv_iter(v_it.handle()); vv_it; ++ vv_it) { // the iterator is clockwise~~~~~~~
//			NeiVertexIDList.push_back(vv_it.handle().idx());
//			pointA = T_Mesh.point(vv_it.handle());
//			NeiVertexList.push_back(pointA); 
//		}
//		std::reverse(NeiVertexIDList.begin(), NeiVertexIDList.end());
//		std::reverse(NeiVertexList.begin(), NeiVertexList.end());
//
//		//n(A) = pb*pc + pc*pd + pd*pe + pe*pf + pf*pg + pg*pb
//		OpenMesh::Vec3f NA; NA.vectorize(0.0);
//		for (int i = 0; i < NeiVertexList.size(); ++ i) { // need to add boundary condition judgment~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			NA += OpenMesh::cross(NeiVertexList[i], NeiVertexList[(i+1+NeiVertexList.size())%NeiVertexList.size()]);
//		}
//		double NAx = NA[0]; double NAy = NA[1]; double NAz = NA[2]; 
//		double NAx_2 = NAx*NAx; double NAy_2 = NAy*NAy; double NAz_2 = NAz*NAz;
//		double LNA = NA.norm(); double LNA_2 = std::pow(LNA,2); double LNA_3 = std::pow(LNA,3); double LNA_4 = std::pow(LNA, 4);
//		OpenMesh::Vec3f N_NA = NA; N_NA.normalize(); // normalized NA;
//		OpenMesh::Vec3f T_NA(m_vertices[vertex_id].ps_normal_x,m_vertices[vertex_id].ps_normal_y,m_vertices[vertex_id].ps_normal_z); T_NA.normalize();
//
//		// the normal term, \beta/2 \|n(u)-n(u)'\|^2
//		for (int i = 0; i < NeiVertexList.size(); ++ i) {
//			// to calculate the coefficient of point c based on the normal of point a, pb*pc+pc*pd
//			OpenMesh::Vec3f pb = NeiVertexList[(i-1+NeiVertexList.size())%NeiVertexList.size()];
//			OpenMesh::Vec3f pd = NeiVertexList[(i+1+NeiVertexList.size())%NeiVertexList.size()];
//			OpenMesh::Vec3f pc = NeiVertexList[(i+0+NeiVertexList.size())%NeiVertexList.size()];
//			double pbx = pb[0]; double pby = pb[1]; double pbz = pb[2];
//			double pdx = pd[0]; double pdy = pd[1]; double pdz = pd[2];
//			double pcx = pc[0]; double pcy = pc[1]; double pcz = pc[2];
//			int cur_vid = NeiVertexIDList[i]; // for the center point pc
//
//			OpenMesh::Vec3f Re, Rest_Norm; Re.vectorize(0.0);
//			for (int j = 0; j < NeiVertexList.size(); ++ j) {
//				if (j != i && (j+1+NeiVertexList.size())%NeiVertexList.size() != i) {
//					Re += OpenMesh::cross(NeiVertexList[j], NeiVertexList[(j+1+NeiVertexList.size())%NeiVertexList.size()]);
//				}
//			}
//			double Rex = Re[0]; double Rey = Re[1]; double Rez = Re[2];
//
//			// construct the jacobian matrix mat_J and f matrix mat_f, using complex version
//			// JtN =(jacobian(N,C)./LNA - transpose(NA)*NA*jacobian(N,C)./LNA._3); // this is the correct normal
//			mat_J(vertex_id*3+0, cur_vid+m_vnum*0) = beta_scale*((NAx*NAz*(pby - pdy) - NAx*NAy*(pbz - pdz))/LNA_3);
//			mat_J(vertex_id*3+0, cur_vid+m_vnum*1) = beta_scale*(((pbz - pdz)*NAx_2 - NAz*(pbx - pdx)*NAx)/LNA_3 - (pbz - pdz)/LNA);
//			mat_J(vertex_id*3+0, cur_vid+m_vnum*2) = beta_scale*((pby - pdy)/LNA - ((pby - pdy)*NAx_2 - NAy*(pbx - pdx)*NAx)/LNA_3);
//
//			mat_J(vertex_id*3+1, cur_vid+m_vnum*0) = beta_scale*((pbz - pdz)/LNA - ((pbz - pdz)*NAy_2 - NAz*(pby - pdy)*NAy)/LNA_3);
//			mat_J(vertex_id*3+1, cur_vid+m_vnum*1) = beta_scale*(-(NAy*NAz*(pbx - pdx) - NAx*NAy*(pbz - pdz))/LNA_3);
//			mat_J(vertex_id*3+1, cur_vid+m_vnum*2) = beta_scale*(((pbx - pdx)*NAy_2 - NAx*(pby - pdy)*NAy)/LNA_3 - (pbx - pdx)/LNA);
//
//			mat_J(vertex_id*3+2, cur_vid+m_vnum*0) = beta_scale*(((pby - pdy)*NAz_2 - NAy*(pbz - pdz)*NAz)/LNA_3 - (pby - pdy)/LNA);
//			mat_J(vertex_id*3+2, cur_vid+m_vnum*1) = beta_scale*((pbx - pdx)/LNA - ((pbx - pdx)*NAz_2 - NAx*(pbz - pdz)*NAz)/LNA_3);
//			mat_J(vertex_id*3+2, cur_vid+m_vnum*2) = beta_scale*((NAy*NAz*(pbx - pdx) - NAx*NAz*(pby - pdy))/LNA_3);
//
//		}
//		mat_f[vertex_id*3+0] = beta_scale*(N_NA[0] - m_vertices[vertex_id].ps_normal_x);
//		mat_f[vertex_id*3+1] = beta_scale*(N_NA[1] - m_vertices[vertex_id].ps_normal_y);
//		mat_f[vertex_id*3+2] = beta_scale*(N_NA[2] - m_vertices[vertex_id].ps_normal_z);
//
//		// the fidelity term, \alpha/2 \|u-f\|^2
//		mat_J(vertex_id*3+0+m_vnum*3, vertex_id+m_vnum*0) = fid_scale*1.0;
//		mat_J(vertex_id*3+1+m_vnum*3, vertex_id+m_vnum*1) = fid_scale*1.0;
//		mat_J(vertex_id*3+2+m_vnum*3, vertex_id+m_vnum*2) = fid_scale*1.0;
//
//		mat_f[vertex_id*3+0+m_vnum*3] = fid_scale*(T_Mesh.point(v_it.handle())[0] - m_vertices[vertex_id].ref_x);
//		mat_f[vertex_id*3+1+m_vnum*3] = fid_scale*(T_Mesh.point(v_it.handle())[1] - m_vertices[vertex_id].ref_y);
//		mat_f[vertex_id*3+2+m_vnum*3] = fid_scale*(T_Mesh.point(v_it.handle())[2] - m_vertices[vertex_id].ref_z);
//
//	}
//	if (IsComplexVersion) {
//		// the tvn term, \gamma/2 \|(p+\lambda^k/\gamma)-\nabla n(u)\|^2, this is defined on each triangle
//		vector<double>nux, nuy, nuz; nux.resize(m_vnum); nuy.resize(m_vnum); nuz.resize(m_vnum);
//		T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
//		T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
//		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) { 
//			nux[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[0];
//			nuy[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[1];
//			nuz[v_it.handle().idx()] = T_Mesh.normal(v_it.handle()).data()[2];
//		}
//		GetPPIGradients (m_vnum,nux);
//		GetPPIGradients1(m_vnum,nuy);
//		GetPPIGradients2(m_vnum,nuz);
//		gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
//		for (int i = 0; i < m_trinum; ++ i) {
//			// \partial grad.x/\partial nx \times \partial nx/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+0)); ++ it) {
//				if (mat_J(i*9+0+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+0+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.x*pen_scale*mat_J(m_triangles[i].ver0*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+0+m_vnum*6, it.index()) += m_TPPIBG[i].v0.x*pen_scale*mat_J(m_triangles[i].ver0*3+0, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+0)); ++ it) {
//				if (mat_J(i*9+0+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+0+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.x*pen_scale*mat_J(m_triangles[i].ver1*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+0+m_vnum*6, it.index()) += m_TPPIBG[i].v1.x*pen_scale*mat_J(m_triangles[i].ver1*3+0, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+0)); ++ it) {
//				if (mat_J(i*9+0+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+0+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.x*pen_scale*mat_J(m_triangles[i].ver2*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+0+m_vnum*6, it.index()) += m_TPPIBG[i].v2.x*pen_scale*mat_J(m_triangles[i].ver2*3+0, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+0+m_vnum*6] = pen_scale*(m_triangles[i].grad.x  - (lambda_x[i].x/penParam + px[i].x));
//
//			// \partial grad.y/\partial nx \times \partial nx/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+0)); ++ it) {
//				if (mat_J(i*9+1+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+1+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.y*pen_scale*mat_J(m_triangles[i].ver0*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+1+m_vnum*6, it.index()) += m_TPPIBG[i].v0.y*pen_scale*mat_J(m_triangles[i].ver0*3+0, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+0)); ++ it) {
//				if (mat_J(i*9+1+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+1+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.y*pen_scale*mat_J(m_triangles[i].ver1*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+1+m_vnum*6, it.index()) += m_TPPIBG[i].v1.y*pen_scale*mat_J(m_triangles[i].ver1*3+0, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+0)); ++ it) {
//				if (mat_J(i*9+1+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+1+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.y*pen_scale*mat_J(m_triangles[i].ver2*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+1+m_vnum*6, it.index()) += m_TPPIBG[i].v2.y*pen_scale*mat_J(m_triangles[i].ver2*3+0, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+1+m_vnum*6] = pen_scale*(m_triangles[i].grad.y  - (lambda_x[i].y/penParam + px[i].y));
//
//			// \partial grad.z/\partial nx \times \partial nx/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+0)); ++ it) {
//				if (mat_J(i*9+2+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+2+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.z*pen_scale*mat_J(m_triangles[i].ver0*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+2+m_vnum*6, it.index()) += m_TPPIBG[i].v0.z*pen_scale*mat_J(m_triangles[i].ver0*3+0, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+0)); ++ it) {
//				if (mat_J(i*9+2+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+2+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.z*pen_scale*mat_J(m_triangles[i].ver1*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+2+m_vnum*6, it.index()) += m_TPPIBG[i].v1.z*pen_scale*mat_J(m_triangles[i].ver1*3+0, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+0)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+0)); ++ it) {
//				if (mat_J(i*9+2+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+2+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.z*pen_scale*mat_J(m_triangles[i].ver2*3+0, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+2+m_vnum*6, it.index()) += m_TPPIBG[i].v2.z*pen_scale*mat_J(m_triangles[i].ver2*3+0, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+2+m_vnum*6] = pen_scale*(m_triangles[i].grad.z  - (lambda_x[i].z/penParam + px[i].z));
//			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//			// \partial grad1.x/\partial ny \times \partial ny/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+1)); ++ it) {
//				if (mat_J(i*9+3+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+3+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.x*pen_scale*mat_J(m_triangles[i].ver0*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+3+m_vnum*6, it.index()) += m_TPPIBG[i].v0.x*pen_scale*mat_J(m_triangles[i].ver0*3+1, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+1)); ++ it) {
//				if (mat_J(i*9+3+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+3+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.x*pen_scale*mat_J(m_triangles[i].ver1*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+3+m_vnum*6, it.index()) += m_TPPIBG[i].v1.x*pen_scale*mat_J(m_triangles[i].ver1*3+1, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+1)); ++ it) {
//				if (mat_J(i*9+3+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+3+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.x*pen_scale*mat_J(m_triangles[i].ver2*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+3+m_vnum*6, it.index()) += m_TPPIBG[i].v2.x*pen_scale*mat_J(m_triangles[i].ver2*3+1, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+3+m_vnum*6] = pen_scale*(m_triangles[i].grad1.x - (lambda_y[i].x/penParam + py[i].x));
//
//			// \partial grad1.y/\partial ny \times \partial ny/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+1)); ++ it) {
//				if (mat_J(i*9+4+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+4+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.y*pen_scale*mat_J(m_triangles[i].ver0*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+4+m_vnum*6, it.index()) += m_TPPIBG[i].v0.y*pen_scale*mat_J(m_triangles[i].ver0*3+1, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+1)); ++ it) {
//				if (mat_J(i*9+4+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+4+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.y*pen_scale*mat_J(m_triangles[i].ver1*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+4+m_vnum*6, it.index()) += m_TPPIBG[i].v1.y*pen_scale*mat_J(m_triangles[i].ver1*3+1, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+1)); ++ it) {
//				if (mat_J(i*9+4+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+4+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.y*pen_scale*mat_J(m_triangles[i].ver2*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+4+m_vnum*6, it.index()) += m_TPPIBG[i].v2.y*pen_scale*mat_J(m_triangles[i].ver2*3+1, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+4+m_vnum*6] = pen_scale*(m_triangles[i].grad1.y - (lambda_y[i].y/penParam + py[i].y));
//
//			// \partial grad1.z/\partial ny \times \partial ny/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+1)); ++ it) {
//				if (mat_J(i*9+5+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+5+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.z*pen_scale*mat_J(m_triangles[i].ver0*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+5+m_vnum*6, it.index()) += m_TPPIBG[i].v0.z*pen_scale*mat_J(m_triangles[i].ver0*3+1, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+1)); ++ it) {
//				if (mat_J(i*9+5+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+5+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.z*pen_scale*mat_J(m_triangles[i].ver1*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+5+m_vnum*6, it.index()) += m_TPPIBG[i].v1.z*pen_scale*mat_J(m_triangles[i].ver1*3+1, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+1)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+1)); ++ it) {
//				if (mat_J(i*9+5+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+5+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.z*pen_scale*mat_J(m_triangles[i].ver2*3+1, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+5+m_vnum*6, it.index()) += m_TPPIBG[i].v2.z*pen_scale*mat_J(m_triangles[i].ver2*3+1, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+5+m_vnum*6] = pen_scale*(m_triangles[i].grad1.z - (lambda_y[i].z/penParam + py[i].z));
//			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//			// \partial grad2.x/\partial nz \times \partial nz/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+2)); ++ it) {
//				if (mat_J(i*9+6+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+6+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.x*pen_scale*mat_J(m_triangles[i].ver0*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+6+m_vnum*6, it.index()) += m_TPPIBG[i].v0.x*pen_scale*mat_J(m_triangles[i].ver0*3+2, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+2)); ++ it) {
//				if (mat_J(i*9+6+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+6+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.x*pen_scale*mat_J(m_triangles[i].ver1*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+6+m_vnum*6, it.index()) += m_TPPIBG[i].v1.x*pen_scale*mat_J(m_triangles[i].ver1*3+2, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+2)); ++ it) {
//				if (mat_J(i*9+6+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+6+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.x*pen_scale*mat_J(m_triangles[i].ver2*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+6+m_vnum*6, it.index()) += m_TPPIBG[i].v2.x*pen_scale*mat_J(m_triangles[i].ver2*3+2, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+6+m_vnum*6] = pen_scale*(m_triangles[i].grad2.x - (lambda_z[i].x/penParam + pz[i].x));
//
//			// \partial grad2.y/\partial nz \times \partial nz/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+2)); ++ it) {
//				if (mat_J(i*9+7+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+7+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.y*pen_scale*mat_J(m_triangles[i].ver0*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+7+m_vnum*6, it.index()) += m_TPPIBG[i].v0.y*pen_scale*mat_J(m_triangles[i].ver0*3+2, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+2)); ++ it) {
//				if (mat_J(i*9+7+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+7+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.y*pen_scale*mat_J(m_triangles[i].ver1*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+7+m_vnum*6, it.index()) += m_TPPIBG[i].v1.y*pen_scale*mat_J(m_triangles[i].ver1*3+2, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+2)); ++ it) {
//				if (mat_J(i*9+7+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+7+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.y*pen_scale*mat_J(m_triangles[i].ver2*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+7+m_vnum*6, it.index()) += m_TPPIBG[i].v2.y*pen_scale*mat_J(m_triangles[i].ver2*3+2, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+7+m_vnum*6] = pen_scale*(m_triangles[i].grad2.y - (lambda_z[i].y/penParam + pz[i].y));
//
//			// \partial grad2.z/\partial nz \times \partial nz/\partial u
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver0*3+2)); ++ it) {
//				if (mat_J(i*9+8+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+8+m_vnum*6, it.index()) =  m_TPPIBG[i].v0.z*pen_scale*mat_J(m_triangles[i].ver0*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+8+m_vnum*6, it.index()) += m_TPPIBG[i].v0.z*pen_scale*mat_J(m_triangles[i].ver0*3+2, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver1*3+2)); ++ it) {
//				if (mat_J(i*9+8+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+8+m_vnum*6, it.index()) =  m_TPPIBG[i].v1.z*pen_scale*mat_J(m_triangles[i].ver1*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+8+m_vnum*6, it.index()) += m_TPPIBG[i].v1.z*pen_scale*mat_J(m_triangles[i].ver1*3+2, it.index())/beta_scale;
//				}
//			}
//			for (it = vect_const_begin(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+2)); it != vect_const_end(gmm::mat_const_row(mat_J, m_triangles[i].ver2*3+2)); ++ it) {
//				if (mat_J(i*9+8+m_vnum*6, it.index()) < 0.000001) {
//					mat_J(i*9+8+m_vnum*6, it.index()) =  m_TPPIBG[i].v2.z*pen_scale*mat_J(m_triangles[i].ver2*3+2, it.index())/beta_scale;
//				} else {
//					mat_J(i*9+8+m_vnum*6, it.index()) += m_TPPIBG[i].v2.z*pen_scale*mat_J(m_triangles[i].ver2*3+2, it.index())/beta_scale;
//				}
//			}
//			mat_f[i*9+8+m_vnum*6] = pen_scale*(m_triangles[i].grad2.z - (lambda_z[i].z/penParam + pz[i].z));
//		}
//	}
//	
//	return;
//}

void TriangularMesh::MeshRefinement(double alpha, double beta, double eta, double varsigma, bool IsComplexVersion = false)
{
	// mesh refine only by position fidelity term, norm fidelity term and norm difference term
	/*   Parameters setting  */
	double outTole = 1.0e-10;// fidParam = 1000, pld_eta = 0, penParam = 0.01, regParam = 1,
	double var = 2.0, epsilon1 = 1.0e-10, epsilon2 = 1.0e-10, minEdgeLength = 5.0e-3; char buffer[255]; 
	vector<double>ux, uy, uz; ux.resize(m_vnum); uy.resize(m_vnum); uz.resize(m_vnum);
	double *b = new double[m_vnum*3]; double *vx = new double[m_vnum*3]; vector<double> vec_b; vec_b.resize(m_vnum*3);
	numc::SparseSolver solver;	
	MyMesh T_Mesh = this->m_ObjTriMesh; OpenMesh::Vec3f CurV; OpenMesh::IO::Options write_options;
	cout << "Start the mesh refinement process: " << endl;
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		CurV[0] = ux[v_it.handle().idx()]; CurV[1] = uy[v_it.handle().idx()]; CurV[2] = uz[v_it.handle().idx()];
		T_Mesh.set_point(v_it.handle(), CurV);
	}
	vector<vector<double>> energy_result; energy_result.clear();
	double v_energy = CalculateVEnergy(T_Mesh);
	vector<double> vec_n_energy; vec_n_energy.resize(m_vnum); 
	double n_energy = CalculateNEnergy(T_Mesh, vec_n_energy);
	double tv_energy = CalculateVTVEnergy(T_Mesh, ux, uy, uz);
	vector<double> t_energy; t_energy.clear(); t_energy.push_back(0.5*alpha*v_energy); t_energy.push_back(beta*n_energy);
	t_energy.push_back(eta*tv_energy); t_energy.push_back(0.5*alpha* v_energy + beta * n_energy + eta * tv_energy);
	energy_result.push_back(t_energy);
	cout << "The current energy is V:" << 0.5*alpha*v_energy << ", N:" << beta*n_energy << ", TV: " << alpha*tv_energy << "; Sum: " 
		<< 0.5*alpha * v_energy + beta* n_energy + eta * tv_energy << endl;

	sprintf(buffer, "Results\\TVU_MeshResult_%.0f_%.0f_%.3f_%s.off", alpha, beta, eta, IsComplexVersion?"Complex":"Simple");
	if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
		std::cerr << "Cannot write mesh to file " << buffer << std::endl;
	}

	// solve the u-sub problem, using Levenberg–Marquardt Algorithm
	double sigma, eOld, eNew;
	int dimension = m_vnum * 3;
	std::vector<double> m_b(dimension, 0);
	std::vector<double> m_x(dimension, 0);
	std::vector<double> tau(dimension, 0);
	std::vector<OpenMesh::Vec3f> originalPhi(m_vnum);
	RowSparseMatrix mat_J;	std::vector<double> m_f;

	bool OnlyNormTerm = true;
	MR_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, alpha, beta, eta, varsigma);

	eOld = gmm::vect_norm2_sqr(m_f);
	RowSparseMatrix mat_JTJ(dimension, dimension);
	RowSparseMatrix mat_A(dimension, dimension);
	gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
	for (int i = 0; i < dimension; i++) {
		tau[i] = mat_JTJ(i, i);
	}
	gmm::scale(tau, 1.e-5);
	gmm::mult(gmm::transposed(mat_J), m_f, m_b);
	gmm::scale(m_b, -1.0);
	double normG = gmm::vect_norminf(m_b);

	int iteration = 0; int maxStep =100;  // 100
	while (iteration++ < maxStep) {
		cout << ".";
		if (normG < epsilon1) {
			break;
		}
		gmm::copy(mat_JTJ, mat_A);
		for (int i = 0; i < dimension; i++) {
			mat_A(i, i) += tau[i];
		}
		// solve mat_A*m_x = m_b
		numc::RowMat<double> RM_A(dimension, dimension); int r = 0;
		for (gmm::linalg_traits<RowSparseMatrix>::row_iterator rIt = gmm::linalg_traits<RowSparseMatrix>::row_begin(mat_A);
			rIt != gmm::linalg_traits<RowSparseMatrix>::row_end(mat_A); ++ rIt, ++ r) {
				gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it = vect_const_begin(*rIt);
				gmm::linalg_traits< gmm::wsvector<double> >::const_iterator ite = vect_const_end(*rIt);
				for (; it != ite; ++ it) {
					RM_A(r, it.index()) = mat_A(r, it.index());
				}
				b[r] = m_b[r];
		}
		solver.getMatA() = RM_A;
		solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
		solver.init();
		// solve the nonlinear equation to calculate the new vertex position in vx;
		solver.solve(b, vx);
		solver.clear();
		for (int r = 0; r < dimension; ++ r) {
			m_x[r] = vx[r];
		}
		////modify m_x, do not change too much according to the value of MinEdgeLength
		//double maxChange = 0.0;
		//for (int i = 0; i < m_vnum; ++ i) {
		//	OpenMesh::Vec3f deltaP(m_x[i], m_x[i+m_vnum], m_x[i+m_vnum*2]);
		//	maxChange = (deltaP.length() > maxChange)?deltaP.length():maxChange;
		//}
		//double stepSize = (minEdgeLength/maxChange > 1.0)?1.0:(minEdgeLength/maxChange);
		//stepSize = 0.2;
		//gmm::scale(m_x, stepSize);
		////end of modify m_x

		double normH = 0.0, normX = 0.0, l0h = 0;
		normH = gmm::vect_norm2(m_x);
		for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
			int idx = v_it.handle().idx();
			originalPhi[v_it.handle().idx()] = T_Mesh.point(v_it.handle());
			normX += originalPhi[v_it.handle().idx()].sqrnorm();
			for (int i = 0; i < 3; ++i) {
				l0h += 0.5*m_x[i*m_vnum + idx]*(tau[i*m_vnum + idx]*m_x[i*m_vnum + idx] + m_b[i*m_vnum + idx]);
			}
		}
		if (normH < epsilon2*(sqrt(normX) + epsilon2)) {
			break;
		} else {
			//change the coordinate value of each vertex...
			OpenMesh::Vec3f CurV;
			for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
				int vIndex = v_it.handle().idx();
				T_Mesh.point(v_it.handle()) += OpenMesh::Vec3f(m_x[vIndex], m_x[m_vnum + vIndex], m_x[2*m_vnum + vIndex]);
			}
			if (iteration%10 == 0) {
				sprintf(buffer, "Results\\LMResults\\TVU_LMMeshResult_%d_%s.off", iteration, IsComplexVersion?"Complex":"Simple");
				if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
					std::cerr << "Cannot write mesh to file " << buffer << std::endl;
				}
			}
			MR_JacobianMatrix_Construction(T_Mesh, mat_J, m_f, alpha, beta, eta, varsigma);
			eNew = gmm::vect_norm2_sqr(m_f);
			sigma = (eOld - eNew)/l0h;
			if (sigma >= 0) {
				gmm::scale(tau, std::max(1./3., 1. - pow(2.*sigma - 1., 3.)));
				var = 2.0;
				eOld = eNew;
				gmm::mult(gmm::transposed(mat_J), m_f, m_b);
				gmm::scale(m_b, -1.0);
				normG = gmm::vect_norminf(m_b);
				gmm::mult(gmm::transposed(mat_J), mat_J, mat_JTJ);
			} else {
				//cout << " Sigma = " << sigma ;
				gmm::scale(tau, var);
				var *= 2.0;
				for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
					T_Mesh.point(v_it) = originalPhi[v_it.handle().idx()];
				}
			}
		}
	}
	//epsilon1 *= 0.001; epsilon2 *= 0.001;
	sprintf(buffer, "Results\\LMResults\\TVU_LMMeshResult_Final_%s.off", IsComplexVersion?"Complex":"Simple");
	if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
		std::cerr << "Cannot write mesh to file " << buffer << std::endl;
	}
}

void TriangularMesh::MR_JacobianMatrix_Construction(MyMesh& T_Mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, double alpha, double beta, double eta, double varsigma)
{
	// construct the jacobian matrix according to the new model
	mat_J.resize(m_vnum*6+T_Mesh.n_edges()*3, m_vnum*3);	mat_f.resize(m_vnum*6+T_Mesh.n_edges()*3);
	gmm::clear(mat_J);				gmm::clear(mat_f);
	double fid_scale = sqrt(alpha*0.5), beta_scale = sqrt(beta*0.5), nd_scale = sqrt(eta*0.5);

	OpenMesh::Vec3f pointA , pointB , pointC; int ida, idb, idc;
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
		// first for all neighbor triangle, store its vertex in counter-clock wise, A(current vi), B, C
		vector< OpenMesh::Vec3f >	NeiVertexList;		NeiVertexList.clear();
		vector< int >				NeiVertexIDList;	NeiVertexIDList.clear();
		vector< double >			NeiFaceAreaList;	NeiFaceAreaList.clear();
		for (MyMesh::ConstVertexVertexIter vv_it = T_Mesh.cvv_iter(v_it.handle()); vv_it; ++ vv_it) { // the iterator is clockwise~~~~~~~
			NeiVertexIDList.push_back(vv_it.handle().idx());
			pointA = T_Mesh.point(vv_it.handle());
			NeiVertexList.push_back(pointA); 
		}
		std::reverse(NeiVertexIDList.begin(), NeiVertexIDList.end());
		std::reverse(NeiVertexList.begin(), NeiVertexList.end());

		//n(A) = pb*pc + pc*pd + pd*pe + pe*pf + pf*pg + pg*pb
		OpenMesh::Vec3f NA; NA.vectorize(0.0);
		for (int i = 0; i < NeiVertexList.size(); ++ i) { // need to add boundary condition judgment~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			NA += OpenMesh::cross(NeiVertexList[i], NeiVertexList[(i+1+NeiVertexList.size())%NeiVertexList.size()]);
		}
		double NAx = NA[0]; double NAy = NA[1]; double NAz = NA[2]; 
		double NAx_2 = NAx*NAx; double NAy_2 = NAy*NAy; double NAz_2 = NAz*NAz;
		double LNA = NA.norm(); double LNA_2 = std::pow(LNA,2); double LNA_3 = std::pow(LNA,3); double LNA_4 = std::pow(LNA, 4);
		OpenMesh::Vec3f N_NA = NA; N_NA.normalize(); // normalized NA;
		OpenMesh::Vec3f T_NA(m_vertices[vertex_id].ps_normal_x,m_vertices[vertex_id].ps_normal_y,m_vertices[vertex_id].ps_normal_z); T_NA.normalize();

		// the normal term, \beta/2 \|n(u)-n(u)'\|^2
		for (int i = 0; i < NeiVertexList.size(); ++ i) {
			// to calculate the coefficient of point c based on the normal of point a, pb*pc+pc*pd
			OpenMesh::Vec3f pb = NeiVertexList[(i-1+NeiVertexList.size())%NeiVertexList.size()];
			OpenMesh::Vec3f pd = NeiVertexList[(i+1+NeiVertexList.size())%NeiVertexList.size()];
			OpenMesh::Vec3f pc = NeiVertexList[(i+0+NeiVertexList.size())%NeiVertexList.size()];
			double pbx = pb[0]; double pby = pb[1]; double pbz = pb[2];
			double pdx = pd[0]; double pdy = pd[1]; double pdz = pd[2];
			double pcx = pc[0]; double pcy = pc[1]; double pcz = pc[2];
			int cur_vid = NeiVertexIDList[i]; // for the center point pc

			OpenMesh::Vec3f Re, Rest_Norm; Re.vectorize(0.0);
			for (int j = 0; j < NeiVertexList.size(); ++ j) {
				if (j != i && (j+1+NeiVertexList.size())%NeiVertexList.size() != i) {
					Re += OpenMesh::cross(NeiVertexList[j], NeiVertexList[(j+1+NeiVertexList.size())%NeiVertexList.size()]);
				}
			}
			double Rex = Re[0]; double Rey = Re[1]; double Rez = Re[2];

			// construct the jacobian matrix mat_J and f matrix mat_f, using complex version
			// JtN =(jacobian(N,C)./LNA - transpose(NA)*NA*jacobian(N,C)./LNA._3); // this is the correct normal
			mat_J(vertex_id*3+0, cur_vid+m_vnum*0) = beta_scale*((NAx*NAz*(pby - pdy) - NAx*NAy*(pbz - pdz))/LNA_3);
			mat_J(vertex_id*3+0, cur_vid+m_vnum*1) = beta_scale*(((pbz - pdz)*NAx_2 - NAz*(pbx - pdx)*NAx)/LNA_3 - (pbz - pdz)/LNA);
			mat_J(vertex_id*3+0, cur_vid+m_vnum*2) = beta_scale*((pby - pdy)/LNA - ((pby - pdy)*NAx_2 - NAy*(pbx - pdx)*NAx)/LNA_3);

			mat_J(vertex_id*3+1, cur_vid+m_vnum*0) = beta_scale*((pbz - pdz)/LNA - ((pbz - pdz)*NAy_2 - NAz*(pby - pdy)*NAy)/LNA_3);
			mat_J(vertex_id*3+1, cur_vid+m_vnum*1) = beta_scale*(-(NAy*NAz*(pbx - pdx) - NAx*NAy*(pbz - pdz))/LNA_3);
			mat_J(vertex_id*3+1, cur_vid+m_vnum*2) = beta_scale*(((pbx - pdx)*NAy_2 - NAx*(pby - pdy)*NAy)/LNA_3 - (pbx - pdx)/LNA);

			mat_J(vertex_id*3+2, cur_vid+m_vnum*0) = beta_scale*(((pby - pdy)*NAz_2 - NAy*(pbz - pdz)*NAz)/LNA_3 - (pby - pdy)/LNA);
			mat_J(vertex_id*3+2, cur_vid+m_vnum*1) = beta_scale*((pbx - pdx)/LNA - ((pbx - pdx)*NAz_2 - NAx*(pbz - pdz)*NAz)/LNA_3);
			mat_J(vertex_id*3+2, cur_vid+m_vnum*2) = beta_scale*((NAy*NAz*(pbx - pdx) - NAx*NAz*(pby - pdy))/LNA_3);

		}
		mat_f[vertex_id*3+0] = beta_scale*(N_NA[0] - m_vertices[vertex_id].ps_normal_x);
		mat_f[vertex_id*3+1] = beta_scale*(N_NA[1] - m_vertices[vertex_id].ps_normal_y);
		mat_f[vertex_id*3+2] = beta_scale*(N_NA[2] - m_vertices[vertex_id].ps_normal_z);

		// the fidelity term, \alpha/2 \|u-f\|^2
		mat_J(vertex_id*3+0+m_vnum*3, vertex_id+m_vnum*0) = fid_scale*1.0;
		mat_J(vertex_id*3+1+m_vnum*3, vertex_id+m_vnum*1) = fid_scale*1.0;
		mat_J(vertex_id*3+2+m_vnum*3, vertex_id+m_vnum*2) = fid_scale*1.0;

		mat_f[vertex_id*3+0+m_vnum*3] = fid_scale*(T_Mesh.point(v_it.handle())[0] - m_vertices[vertex_id].ref_x);
		mat_f[vertex_id*3+1+m_vnum*3] = fid_scale*(T_Mesh.point(v_it.handle())[1] - m_vertices[vertex_id].ref_y);
		mat_f[vertex_id*3+2+m_vnum*3] = fid_scale*(T_Mesh.point(v_it.handle())[2] - m_vertices[vertex_id].ref_z);

	}
	
	// the normal difference term, eta*\omega_ij\|n_i - n_j\|^2
	T_Mesh.request_face_normals(); T_Mesh.update_face_normals();
	T_Mesh.request_vertex_normals(); T_Mesh.update_vertex_normals();
	gmm::linalg_traits< gmm::wsvector<double> >::const_iterator it, ite;
	for (MyMesh::EdgeIter e_it = T_Mesh.edges_begin(); e_it != T_Mesh.edges_end(); ++ e_it) {
		MyMesh::HalfedgeHandle h1 = T_Mesh.halfedge_handle(e_it,0);
		MyMesh::HalfedgeHandle h2 = T_Mesh.halfedge_handle(e_it,1);
		int iid = T_Mesh.to_vertex_handle(h1).idx();
		int jid = T_Mesh.to_vertex_handle(h2).idx();

		double omegaij = std::exp(-std::pow((m_vertices[iid].intensity - m_vertices[jid].intensity), 2.0)/varsigma);

		for (int k = 0; k < 3; k++) {
			for (it = vect_const_begin(gmm::mat_const_row(mat_J, iid*3+k)); it != vect_const_end(gmm::mat_const_row(mat_J, iid*3+k)); ++ it) {
				mat_J(e_it.handle().idx()*3+k+m_vnum*6, it.index()) += nd_scale*omegaij*mat_J(iid*3+k, it.index())/beta_scale;
			}
			for (it = vect_const_begin(gmm::mat_const_row(mat_J, jid*3+k)); it != vect_const_end(gmm::mat_const_row(mat_J, jid*3+k)); ++ it) {
				mat_J(e_it.handle().idx()*3+k+m_vnum*6, it.index()) += nd_scale*omegaij*mat_J(jid*3+k, it.index())/beta_scale;
			}
			mat_f[e_it.handle().idx()*3+k+m_vnum*6] = nd_scale*omegaij*(T_Mesh.normal(MyMesh::VertexHandle(iid)).data()[k]- T_Mesh.normal(MyMesh::VertexHandle(jid)).data()[k]);
		}

	}
}

void TriangularMesh::UpdateMeshInfo(MyMesh &T_Mesh)
{
	T_Mesh.update_face_normals();  T_Mesh.update_vertex_normals();
	for (MyMesh::VertexIter v_it = T_Mesh.vertices_begin(); v_it != T_Mesh.vertices_end(); ++ v_it) {
		m_vertices[v_it.handle().idx()].x = T_Mesh.point(v_it.handle()).data()[0];
		m_vertices[v_it.handle().idx()].y = T_Mesh.point(v_it.handle()).data()[1];
		m_vertices[v_it.handle().idx()].z = T_Mesh.point(v_it.handle()).data()[2];
	}
	BuildTPPIBG();//to compute the gradient of basis functions of primal-primal interpolation.
	//  WriteTPPIBG();
	ComputeBCDualArea();
	if (AnisotropicLaplace) {
		CalculateVertexVoronoiArea();
	}
	BuildTrianglesArea();
}

void TriangularMesh::CalculateVertexVoronoiArea()
{
	// calculate the voronoi area for each vertex
	for (MyMesh::VertexIter v_it = m_ObjTriMesh.vertices_begin(); v_it != m_ObjTriMesh.vertices_end(); ++ v_it) {
		int vertex_id = v_it.handle().idx(); m_vertices[vertex_id].Voronoi_Area = 0.0;
		for (MyMesh::VertexEdgeIter ve_it = m_ObjTriMesh.ve_iter(v_it); ve_it; ++ ve_it) {
			MyMesh::FaceHandle fh1 = m_ObjTriMesh.face_handle(m_ObjTriMesh.halfedge_handle(ve_it,0)); 
			MyMesh::FaceHandle fh2 = m_ObjTriMesh.face_handle(m_ObjTriMesh.halfedge_handle(ve_it,1));
			int iid = m_ObjTriMesh.to_vertex_handle(m_ObjTriMesh.halfedge_handle(ve_it,0)).idx();
			int jid = m_ObjTriMesh.to_vertex_handle(m_ObjTriMesh.halfedge_handle(ve_it,1)).idx();

			MyMesh::ConstFaceVertexIter cfv_it = m_ObjTriMesh.cfv_iter(fh1);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			}  int k1id = cfv_it.handle().idx();
			cfv_it = m_ObjTriMesh.cfv_iter(fh2);
			while (cfv_it.handle().idx() == iid || cfv_it.handle().idx() == jid) {
				++cfv_it;
			}  int k2id = cfv_it.handle().idx();

			OpenMesh::Vec3f pPi  = m_ObjTriMesh.point(MyMesh::VertexHandle(iid)),  pPj  = m_ObjTriMesh.point(MyMesh::VertexHandle(jid));
			OpenMesh::Vec3f pPk1 = m_ObjTriMesh.point(MyMesh::VertexHandle(k1id)), pPk2 = m_ObjTriMesh.point(MyMesh::VertexHandle(k2id));
			double cotaij = OpenMesh::dot(pPi-pPk1, pPj-pPk1)/OpenMesh::cross(pPi-pPk1, pPj-pPk1).norm();
			double cotbij = OpenMesh::dot(pPi-pPk2, pPj-pPk2)/OpenMesh::cross(pPi-pPk2, pPj-pPk2).norm();
			m_vertices[vertex_id].Voronoi_Area += (cotaij + cotbij)*(pPi-pPj).sqrnorm()/8.0;
		}
	}
	return;
}