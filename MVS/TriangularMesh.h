// TriangularMesh.h: interface for the TriangularMesh class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRIANGULARMESH_H__487DD18F_542D_435B_961C_376988EE3253__INCLUDED_)
#define AFX_TRIANGULARMESH_H__487DD18F_542D_435B_961C_376988EE3253__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//#include <gl/GLU.h>
#define _SCL_SECURE_NO_DEPRECATE
#include <gmm/gmm.h>
#include <gmm/gmm_matrix.h>
#include <gmm/gmm_superlu_interface.h>
typedef gmm::row_matrix< gmm::wsvector<double> > RowSparseMatrix;
typedef gmm::row_matrix< std::vector<double> > RowDenseMatrix;
typedef gmm::dense_matrix<double> DenseMatrix;

#include "POINT3D.h"
#include "TypeDef.h"
#include "mkl_addon.h"
#include "speigen.h"
#include <OpenMesh/Core/io/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
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
#pragma comment(lib,"OpenMeshCored.lib")
#pragma comment(lib,"OpenMeshToolsd.lib")

#include <engine.h>
#include <matrix.h>
#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libmex.lib")
#pragma comment(lib,"libmat.lib")
#pragma comment(lib,"libut.lib")
#pragma comment(lib,"libeng.lib")
#pragma comment(lib,"mclcommain.lib")
extern Engine *m_ep;
extern bool ScaleDelta;
extern bool AnisotropicLaplace;
extern bool RecordColor;

class TriangularMesh  
{
public:
	TriangularMesh();
	virtual ~TriangularMesh();

private:
	void ALM_MulRegionLabelling_Binarization(double *u, unsigned int nRegion, unsigned long nver);
	void ALM_MulRegionLabelling_Project2K(double *u, unsigned int nRegion, unsigned long nver);
	void ALM_MulRegionLabelling_sDev_colorimage(double *u, double *s, unsigned long nRegion, unsigned long nver);
	void ALM_MulRegionLabelling_sDev_grayimage(double * u, double *  s, unsigned long nRegion, unsigned long nver);
	void ALM_MulRegionLabelling_sDev_colorimage_givenMean(double *meanr, double *meang, double *meanb, double *s, unsigned long nRegion, unsigned long nver);
	void ALM_MulRegionLabelling_sDev_grayimage_givenMean(double * mean, double *  s, unsigned long nRegion, unsigned long nver);
	void AddSaltPepperNoise(float level);
	bool HasObtuseAngle();
	void AddGaussianNoise2VertexNormal(double variance);
	void AddGaussianNoise2VertexGeometry(double variance);
	void AddGaussianNoise(double variance);
	double GaussianRandom();
	void GetMeshInformation();
	double GetLargestEdgeLength();
	double GetSmallestLargestTriAreaRatio();
	double GetSmallestEdgeLength();
	double GetSmallestEdgeLengthRatio();
	double GetSmallestAngleRatio();
	double g2HeatCoeff(double s);
	double g1HeatCoeff(double s);
	double gHeatCoeff(double s);
	void Contour2Intensity(double imin, double imax);
	double SignMin2(double a, double b);
	void ResetContourTimestep_3();
	void ResetContourTimestep_2();
	void WriteWeightedCurveLength_Forobservation();
	bool BackwardsContour();
	void ResetContourTimestep();
	void Design_e1e2(TRIANGLE t, VECTOR3D *e1, VECTOR3D *e2);
	void Translate(VERTEX3D *v, VECTOR3D d, double len);
	double LengthVector3D(VECTOR3D v);
	VERTEX3D GetFootInATriangle(TRIANGLEWITHCOOR twc, int ver);
	void Unfolding(TRIANGLEWITHCOOR s,TRIANGLEWITHCOOR *t);
	void FMM_LocalSolver(int vertex, int triangle, double F);
	void FMM_LocalSolver_Acute(TRIANGLEWITHCOOR twc, double u0, double u1, double *u2, double F);
	void Reinitialize_Contour(double time);
	double ComputeWeightedCurveLength();
	void WriteMeshToPovray();
	void WriteContourToPovray();
	void WriteWeightedCurveLength();
	void ComputeAveragedVertexAbGradient();
	void BuildEdgeIndicator(unsigned char type);
	bool ContouringATriangle(POINT3d tmp1,POINT3d tmp2,POINT3d tmp3,POINT3d *linend1,POINT3d *linend2);
	void DrawContour_MulRegionLabelling();
	void DrawContour();
	void DrawMesh();
	void BuildONEDISKSANISOTROPIC();
	void BuildPrincipleDirections();
	void Cramer3(MATRIX3BY3 m,double r1,double r2,double r3,double *x1,double *x2,double *x3);
	double Determinant3(MATRIX3BY3 m);
	void Cramer3(MATRIX3BY3 m,double r[],double x[]);
	VECTOR3D GetCoProjection(VECTOR3D v,VECTOR3D r);
	double GetProjection(VECTOR3D v,VECTOR3D r);
	void GetPrincipleDirection(int iver,VECTOR3D *e1,VECTOR3D *e2);
	void NormalizeVector3D(double *xx,double *yy,double *zz);
	void BuildTrianglesArea();
	void GetPPIGradients(long n,double q[]);
	void GetPPIGradients1(long n,double q[]);//
	void GetPPIGradients2(long n,double q[]);//
	void GetPPIGradients(double q[],VECTOR3D *g);//
	void GetDivergence(long ntri,VECTOR3D vf[]);
	void GetDivergence(VECTOR3D vf[],double *div);//
	void ImplicitL2RDTMeinhardtStrips_Mein();
	void ImplicitL2RDTTuring_GiererMein();
	void ImplicitL2RDTTuring_GrayScott();
	void ImplicitL2RDTTuring_Brusselator();
	void ImplicitL2RDTTuring_Turing();
	double RDTexturing_BrusselatorG(double a, double b, double s, double alpha, double beta);
	double RDTexturing_BrusselatorF(double a,double b,double  s,double alpha,double beta);
	void AddUniformNoise(float mag);
	void WriteData();
	double RDTexturing_TurkG(double a,double b,double rm);
	double RDTexturing_TurkF(double a,double b,double rm);

	void ALM_TVL2Denoising_usub_nr_linbcg(unsigned long n,double b[],double x[],double fidParam, double penParam,int itol,double tol,int itmax,int* iter,double* err);
	void ALM_TVL2Denoising_usub_nr_atimes(unsigned long n,double x[],double r[],double fidParam, double penParam,int itrnsp);
	void ALM_TVL2Denoising_usub_nr_dsprsax(unsigned long n,double x[],double b[],double fidParam, double penParam);
   	void ALM_TVL2Denoising_usub_nr_asolve(unsigned long n,double b[],double x[],double fidParam, double penParam,int itrnsp);

    void ALM_vTVL2Denoising_usub_nr_linbcg(unsigned long n,double b[],double x[],double fidParam, double penParam,int itol,double tol,int itmax,int* iter,double* err);
	void ALM_vTVL2Denoising_usub_nr_atimes(unsigned long n,double x[],double r[],double fidParam, double penParam,int itrnsp);
	void ALM_vTVL2Denoising_usub_nr_dsprsax(unsigned long n,double x[],double b[],double fidParam, double penParam);
   	void ALM_vTVL2Denoising_usub_nr_asolve(unsigned long n,double b[],double x[],double fidParam, double penParam,int itrnsp);

    void ALM_TVgL1_usub_nr_linbcg(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itol,double tol,int itmax,int* iter,double* err);
	void ALM_TVgL1_usub_nr_atimes(unsigned long n,double x[],double r[],double penParam_p, double penParam_z,int itrnsp);
	void ALM_TVgL1_usub_nr_dsprsax(unsigned long n,double x[],double b[],double penParam_p, double penParam_z);
   	void ALM_TVgL1_usub_nr_asolve(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itrnsp);

    void ALM_MulRegionLabelling_usub_nr_linbcg(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itol,double tol,int itmax,int* iter,double* err);
	void ALM_MulRegionLabelling_usub_nr_atimes(unsigned long n,double x[],double r[],double penParam_p, double penParam_z,int itrnsp);
	void ALM_MulRegionLabelling_usub_nr_dsprsax(unsigned long n,double x[],double b[],double penParam_p, double penParam_z);
   	void ALM_MulRegionLabelling_usub_nr_asolve(unsigned long n,double b[],double x[],double penParam_p, double penParam_z,int itrnsp);

	void nr_linbcg(unsigned long n,double b[],double x[],double tstep,int itol,double tol,int itmax,int* iter,double* err);
	double nr_snrm(unsigned long n,double sx[],int itol);
	void nr_asolve(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_L2Diff(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_L2Denois(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_L2Inpaint(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_TVDiff(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_NonlinearDiff(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_GeodesicCurvatureFlow(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_GeodesicWeightedCurvatureFlow(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_TVDenois(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_TVInpaint(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_L2RDTextur(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_TVRDTextur(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_L2DiffDirectionalData(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_TVDiffDirectionalData(unsigned long n,double b[],double x[],double tstep,int itrnsp);
    void nr_asolve_AnisotropicRDTextur(unsigned long n,double b[],double x[],double tstep,int itrnsp);
	void nr_asolve_AnisotropicDAntialiasing(unsigned long n,double b[],double x[],double tstep,int itrnsp);//added on 4/4/08.

	void nr_dsprstx(unsigned long n,double x[],double b[],double tstep);
	//attention: the coefficient matrices of gcf are not symmetric.
	void nr_dsprstx_GeodesicCurvatureFlow(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprstx_WeightedGeodesicCurvatureFlow(unsigned long n,double x[],double b[],double tstep);
	void nr_atimes(unsigned long n,double x[],double r[],double tstep,int itrnsp);
	void nr_dsprsax(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_L2Diff(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_L2Denois(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_L2Inpaint(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_TVDiff(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_NonlinearDiff(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_TVDenois(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_TVInpaint(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_L2RDTextur(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_TVRDTextur(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_L2DiffDirectionalData(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_TVDiffDirectionalData(unsigned long n,double x[],double b[],double tstep);
    void nr_dsprsax_AnisotropicRDTextur(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_AnisotropicDAntialiasing(unsigned long n,double x[],double b[],double tstep);//added on 4/4/08.
	void nr_dsprsax_GeodesicCurvatureFlow(unsigned long n,double x[],double b[],double tstep);
	void nr_dsprsax_GeodesicWeightedCurvatureFlow(unsigned long n,double x[],double b[],double tstep);

	bool HasOrvInOnering(ONERING or,ONERINGVERTEX orv);
	void BuildOnerings();
	void ComputeBCDualArea();
	void BuildONEDISKS();
	VECTOR3D GetVector3Dfrom2Vertices(VERTEX3D begin,VERTEX3D end);
	VECTOR3D GetBCDNormal(int v,EDGE e,TRIANGLE t);
	VECTOR3D GetBCDNormal(VERTEX3D v,EDGEWITHCOOR ewc,TRIANGLEWITHCOOR twc);
	VERTEX3D GetBaryCenter(VERTEX3D v);
	void WriteTPPIBG();
	void BuildTPPIBG();
	VECTOR3D FindPPIGradient(int it, unsigned int dataType);//dataType指示对该顶点的什么量求梯度，如颜色信息的r分量？g分量？b分量？
	VECTOR3D FindPPIGradient(int it);
	void FindPPIGradient(TRIANGLE t, VECTOR3D *g, unsigned int dataType);//dataType指示对该顶点的什么量求梯度，如颜色信息的r分量？g分量？b分量？
	void FindPPIGradient(TRIANGLE t, VECTOR3D *g);
	void FindPPIBGradient(int ver, TRIANGLE t, VECTOR3D *g, double *h);
	void BuildBaryCentersForTriangles();
	VERTEX3D GetBaryCenter(TRIANGLEWITHCOOR twc);
	VERTEX3D GetBaryCenter(TRIANGLE t);
	COLORRGB GetColorFromPoint3D(POINT3d p);
	VERTEX3D GetVertex3DFromPoint3D(POINT3d p);
	VERTEX3D GetBaryCenter(EDGEWITHCOOR ewc);
	VERTEX3D GetBaryCenter(EDGE e);
	void BuildCCDualVerticesForTriangles();
	double Norm(VECTOR3D v);
	double NormSquare(VECTOR3D v);
	double Norm(VECTOR3D V,double alpha);
	VERTEX3D GetCCDual(TRIANGLE t);
	bool IsAdjointTriangle2(TRIANGLE t,EDGE e,TRIANGLE at);
	bool IsAdjointTriangle(TRIANGLE t,EDGE e,TRIANGLE at);
	bool HasThisEdge2(TRIANGLE t,EDGE e);
	EDGE FindAdjointCoTriEdge(TRIANGLE t,EDGE e,int v0);
	int FindOppositeVertex(TRIANGLE t,EDGE e);
	int FindAnotherVertex(TRIANGLE t,int v0);
	bool HasThisEdge(TRIANGLE t,EDGE e);
	bool IsSharingEdge(TRIANGLE t1,TRIANGLE t2);
	void BuildOneDisks();
	bool WriteNormalsToFile();
	void NormalizeVector3D(VECTOR3D *v);
	void FindNormals();
	
	int m_vnum;
	int m_trinum;
	int m_nTime;
	double m_Lambda;
	double m_tStep;
	double m_alpha;
	unsigned char m_ObjectDrawing;// showing which object to be shown.
	                             // 0: the mesh (maybe colored);
	                             // 1: the contour;
	                             // 2: the mesh and the contour;
	vecDoubles m_WeightedCurveLength;
	vecDoubles m_tSteps;
	vecDoubles m_WeightedCurveLength_Forobservation;// no pop-up for observation.
	vecDoubles m_tSteps_Forobservation;// no pop-up for observation.
	vecDoubles m_CPUtime;

	// the following items describes the mesh quality. The angleratio and edgelengthratio are defined triangle by triangle.
	double m_SmallestAngleRatio;
	double m_SmallestEdgeLengthRatio;
	double m_SmallestEdgeLength;
	double m_LargestEdgeLength;
	double m_SmallestLargestEdgeLengthRatio;
	double m_SmallestLargestTriAreaRatio;// the ratio between the areas of smallest and largest triangle.

	RDTEXTUREPARAM m_RDTp;
	double m_RDTDiffRate;//
	vecPoints m_vertices;
	vecTriangles2 m_triangles;
	vecVertex3D m_trianglesBC;
	vecOneDisks m_onedisks;//store onedisks for all vertices of the mesh: counterclockwise or clockwise order.
    vecONEDISKS m_ONEDISKS;
	vecONEDISKSANISOTROPIC m_ONEDISKSANISOTROPIC;
	vecONERINGS m_ONERINGS;//store one-rings for all vertices of the mesh: for each vertex, store the index and
	                       //lb operator coefficient of vertices of its one-ring.
	                       //here the one-ring does not contain the vertex itself,so the lb operator
	                       //coefficient of itself does not included in the one-ring data. it can be determined
	                       //by the equallity that all coefficients sum to zero.
    vecTPPIBG m_TPPIBG;
    vecDoubles m_tempContour;// used to store current contour data for one-step backwards when neceesory.
	vecCOLORRGB m_tempRGB;//used for source term in denoising algorithms.
	vecDoubles m_trianglesArea;//areas of triangles in the mesh.
	unsigned short m_ProcessingMethod;//1 L2diffusion; 2 TVdiffusion;
	                                  //3 L2denoising; 4 TVdenoising;
	                                  //5 L2inpainting;6 TVinpainting;
	                                  //7 L2R-D texturing;8 TVR-D texturing;
	                                  //9 L2diffusion of unit directional data;
	                                  //10 TVdiffusion of unit directional data;
	                                  //11 Anisotropic R-D texturing;
	                                  //12 Geodesic curvature flow of a contour;
	                                  //13 Geodesic Weighted curvature flow with an edge indicator of a contour;
	                                  //14 Scale Space of images by Geodesic curvature flow (in fact the same as 12);
	                                  //16 Anisotropic diffusion for antialiasing.
	                                  //17 General nonlinear diffusion with coefficient function g=g(s).

	void UnifyData();
	VECTOR3D CrossProduct(VECTOR3D v1, VECTOR3D v2);
	double DotProduct(VECTOR3D v1, VECTOR3D v2);
public:

	void ALM_vTVL2NormalFiltering();
    void ALM_TVL2Denoising();
	void ALM_vTVL2Denoising();
	void ALM_TVgL1();
	void ALM_MulRegionLabelling();
	void ALM_MulRegionLabelling_givenMean();

	void SetInitialContourFunction();
	void SetMarchingSteps(int ms);
	void ImplicitL2Diff_Contour(double smoothtime);
	void SetObjectDrawing(unsigned char od);
	void ShiftObjectDrawing();
	void ImplicitScaleSpacebyGCF();
	void ImplicitGeodesicWeightedCurvatureFlow();
	void ImplicitGeodesicCurvatureFlow();
	void ImplicitGeodesicWeightedCurvatureFlow_AdaptiveTimeSteps();// added on 16/04/08.
	void ImplicitGeodesicCurvatureFlow_AdaptiveTimeSteps();// added on 16/04/08.
	void WriteDataToPovray();
	void ImplicitAnisotropicRDTextur();
	void ImplicitAnisotropicDAntialiasing();//added on 4/4/08.
	void ImplicitL2DiffDirecData();
	void ImplicitTVDiffDirecData();
	void SetRDTextureParam(double Da,double Db,double ar,double gr,double decay);
	void ImplicitTVRDTextur();
	void ImplicitTVInpaint();
	void ImplicitTVDenois();
	void ImplicitTVDiff();
	void ImplicitNonlinearDiff();
	void ImplicitL2Inpaint();
	void ImplicitL2Diff();
	void SetProcessingMethod(unsigned short pm);
	void ImplicitL2RDTextur();
	void ImplicitL2Denoising();
	void LBDiffusion2();
	void LBDiffusion();
	void Draw();
	bool LoadFromFile(const char* filename);
	bool SaveToFile(const char* filename);

	// add for alm_mesh_refinement
//#define TEST_MESHREFINE
//#define TVNSMOOTH
	void GetDivergence(long ntri,vector<VECTOR3D>vf);
	void GetPPIGradients(long n, vector<double> q);
	void GetPPIGradients1(long n, vector<double> q);
	void GetPPIGradients2(long n, vector<double> q);
	bool LoadMeshFile(const char* filename);
	bool LoadTargetMeshFile(const char* filename);
	bool LoadVertexColorFile(const char* filename, double m_sigma);
	bool UpdateVertexWeight(double m_sigma);
	bool LoadPSNormalFile(const char* filename);

	void ALM_TVU_MeshRefinement(string meshname, double PosfidParam, double LitfidParam, double pld_eta, double pcd_eta, double fcd_eta, double pnd_eta, double fnd_eta, double varsigma, double pc_eta, double penParam, double regParam, double lapParam, bool UseTVU, bool UseTVNorm, int iter_step, bool UseFaceArea, bool UseMatlabSolver);
	void TV_JacobianMatrix_Construction(MyMesh& T_mesh, RowSparseMatrix& mat_J, DenseMatrix& mat_f, double PosfidParam, double LitfidParam, double pld_eta, double pcd_eta, double fcd_eta,
		double pnd_eta, double fnd_eta, double varsigma, double pc_eta, double penParam, double lapParam, 
		vector<VECTOR3D> &px, vector<VECTOR3D> &py, vector<VECTOR3D> &pz, vector<VECTOR3D> &lambda_x, vector<VECTOR3D> &lambda_y, vector<VECTOR3D> &lambda_z, 
		bool UseTVU, bool UseTVNorm, bool UseFaceArea);
	void CalculateEnergyTerm(DenseMatrix mat_f, vector<double>& energy, bool UsePositionFidelity, bool UseLightFidelity, bool UsePointLightDiff, bool UsePointColor, bool UsePointColorDiff, bool UseFaceColorDiff, bool UsePointNormalDiff, bool UseFaceNormalDiff, bool UseLaplace, bool UseTVNorm, bool UseTVU);
	void CalculateInitialLight(string meshname, double pc_param, double pld_param);
	void Light_JacobianMatrix_Construction(MyMesh& T_mesh, RowSparseMatrix& mat_J, DenseMatrix& mat_f, double pc_param, double pld_param);

	double CalculateVEnergy(MyMesh &T_Mesh, bool UseFaceArea);
	double CalculateNEnergy(MyMesh &T_Mesh, vector<double> &vec_pervertex_energy, bool UseFaceArea);
	double CalculateNDPEnergy(MyMesh &T_Mesh, double varsigma, bool UseFaceArea);
	double CalculateNDFEnergy(MyMesh &T_Mesh, double varsigma, bool UseFaceArea);
	double CalculateVTVEnergy(MyMesh &T_Mesh, vector<double> &ux, vector<double> &uy, vector<double> &uz, bool UseFaceArea);
	double CalculateVTVNormEnergy(MyMesh &T_Mesh, vector<double> &nux, vector<double> &nuy, vector<double> &nuz, bool UseFaceArea);
	double CalculateLaplaceEnergy(MyMesh &T_Mesh, bool UseFaceArea);
	void UpdateMeshInfo(MyMesh &T_Mesh);

	void LM_Testing(bool TestTV, int choice);
	void GradientTesting(double h, bool TestTV, int choice);
	void LM_JacobianMatrix_Construction(MyMesh& T_mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, int choice);
	void LM_TV_JacobianMatrix_Construction(MyMesh& T_mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, bool UseTVNorm);
	void ALM_MeshSmooth();
	void MeshRefinement(double alpha, double beta, double eta, double varsigma, bool IsComplexVersion);
	void MR_JacobianMatrix_Construction(MyMesh& T_mesh, RowSparseMatrix& mat_J, std::vector<double>& mat_f, double alpha, double beta, double eta, double varsigma);
	//void matlabSolvingLinearSystem(gmm::row_matrix< gmm::wsvector<double> >& m_A, std::vector<Scalar>& m_x, std::vector<Scalar>& m_b);
	void CalculateVertexVoronoiArea();
public:
	MyMesh m_ObjTriMesh;
};

#endif // !defined(AFX_TRIANGULARMESH_H__487DD18F_542D_435B_961C_376988EE3253__INCLUDED_)
