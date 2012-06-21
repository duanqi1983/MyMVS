#ifndef __MESHTVREFINE_H__
#define __MESHTVREFINE_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <conio.h>
#include <Windows.h>
#include <vector>
//#include <cgal/Taucs_solver_traits.h>
//#pragma comment(lib,"libmetis.lib")
//#pragma comment(lib,"libtaucs.lib")
//#pragma comment(lib,"vcf2c.lib")
//typedef CGAL::Taucs_solver_traits<float> TaucsSolver;
//typedef TaucsSolver::Matrix TaucsMatrix;
//typedef TaucsSolver::Vector TaucsVector;
using namespace std;
#include <engine.h>
#include <matrix.h>
#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libmex.lib")
#pragma comment(lib,"libmat.lib")
#pragma comment(lib,"libut.lib")
#pragma comment(lib,"libeng.lib")
#pragma comment(lib,"mclcommain.lib")

#include "mkl_addon.h"
#include "speigen.h"
#include <OpenMesh/Core/io/MeshIO.hh>
//#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

#define USE_FITTED_COLOR

struct MyMeshTraits : OpenMesh::DefaultTraits
{
	//typedef OpenMesh::Vec4f Color;
	VertexAttributes( OpenMesh::Attributes::Normal |
		OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
};

//typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyPolyMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<MyMeshTraits>  MyMesh;

class MeshTVRefine
{
public:
	MeshTVRefine();
	//MeshTVRefine(double _alpha, double _beta, double _gamma, double _sigma):m_alpha(_alpha),m_beta(_beta),m_gamma(_gamma),m_sigma(_sigma){
	//	ObjTriMesh.clear();
	//	Vertex_Intensity.clear();
	//	Vertex_PSNormal.clear();
	//	m_BCDArea.clear();
	//	m_TriangleArea.clear();
	//	m_divergence.clear();
	//	m_vertexweight.clear();
	//	m_vertexfaces.clear();
	//	m_edgelaplaceweight.clear();
	//	m_TPPIBG.clear();
	//	m_TriGradient.clear();
	//};
	~MeshTVRefine();
	MeshTVRefine& operator=(const MeshTVRefine &inmtvr);
	MeshTVRefine(const MeshTVRefine &inmtvr);
	bool LoadMeshFile(const char* filename);
	bool LoadColorFile(const char* filename, double m_sigma);
	bool LoadPSNormalFile(const char* filename);
	bool LoadFaceNormalFile(const char* filename);
	void UpdateTriangleArea(); // for each triangular face
	void UpdateBCDArea(); // for each vertex
	bool UpdateVertexWeight(double m_sigma);
	void UpdateEdgeLaplaceWeight();
	float CalculateHalfEdgeLaplaceWeight(MyMesh::HalfedgeHandle hehd);
	bool CalDivengence(const vector< vector<OpenMesh::Vec3f> > &vf); // Given a vector field vf for each triangular on mesh, calculate its divergence for each vertex
	bool CalGradient(const vector<OpenMesh::Vec3f> &vertexvec, vector< vector<OpenMesh::Vec3f> > &out_TriGradient);
	void ALMTVMeshRefine(double m_alpha, double m_beta, double m_gamma, int iter_step, double m_eta);
	void ALMTVMeshRefineByFaceNormal(double m_alpha, double m_beta, double m_gamma, int iter_step);
	void BuildTPPIBG(MyMesh &InMesh, vector< map<int, OpenMesh::Vec3f> > &vec_TPPIBG);
	void FindPPIBGradient(MyMesh &InMesh, const int ver_id, const MyMesh::FaceHandle f_handle, OpenMesh::Vec3f &vec_g, float &h);
	void MeshRefineByNormal(double ialpha, double ibeta);
	void MeshRefineByFaceNormal(double ialpha, double ibeta);

	void LocalNormalRefinement(const double ialpha, const double ibeta, const double igamma);
	void LocalNormalTVRefinement(const double ialpha, const double igamma, const double isigma);

	double CalculateNEnergy(MyMesh &CurMesh, const std::vector<OpenMesh::Vec3f> InputN);
	double CalculateTVEnergy(MyMesh &CurMesh);
	double CalculateVEnergy(MyMesh &CurMesh, const std::vector<OpenMesh::Vec3f> InputV);
//protected:
private:
	MyMesh ObjTriMesh;
	vector<double> Vertex_Intensity;				// used to calculate the weight for each vertex to its 1-ring neighbors.
	vector<OpenMesh::Vec3f> Vertex_PSNormal;
	vector<OpenMesh::Vec3f> Face_Normal;
	vector<double>	m_BCDArea;
	vector<double>	m_TriangleArea;
	vector< vector<double> > m_divergence;			// x,y,z channels
	vector< map<int, OpenMesh::Vec3f> > m_TPPIBG;	// three vertex of each triangular face
	vector< vector<OpenMesh::Vec3f> > m_TriGradient;// for the x, y, z channels
	vector<double>	m_vertexweight;
	vector<int>		m_vertexfaces;
	vector<double>	m_edgelaplaceweight;

	//double m_alpha;					// parameter for fidelity term;
	//double m_beta;					// parameter for normal term
	//double m_gamma;					// parameter for penalty term
	//double m_sigma;					// parameter to calculate vertex weight;
};


#endif __MESHTVREFINE_H__