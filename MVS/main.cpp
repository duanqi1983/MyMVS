#include "daisy/daisy.h"
#include "ModelingOptions.h"
#include "TriangularMesh.h"
#include "parameters.h"
#include "matrix.h"
#include "cminpack.h"
#include "viewPoint.h"
#include "track.h"
#include "lbfgs.h"
#include "readData.h"
#include "LMpart.h"
#include "LBFGSpart.h"
//#include "MeshTVRefine.h"
#include "chooseStereoPairs.h"
#include "testMatchingResult.h"
#include "output.h"

//typedef CGAL::Simple_cartesian<double>     Kernel;
//typedef Kernel::Point_3                    Point_3;
//typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
//typedef Polyhedron::Vertex_iterator        Vertex_iterator;

int nViews;
viewPoint** listViewPoints;
double *** Rij = NULL;
double *** Tij = NULL;
vector<track> tracks;
bool **paired = NULL;
float **** daisyMEM = NULL;
int cntF;
double* camera_dis = NULL;
double *thetas = NULL;
double medianD;
ModelingOptions MOptions;
//pcl::PointCloud<pcl::PointNormal> my_pcl_cloud;
//pcl::PolygonMesh my_pcl_mesh;

MyMesh ObjTriMesh;
char buffer[255];
double SplineFittingDecreasePara = 0.4;
string psnormalfile;
string intensityfile;
bool ScaleDelta = true;
bool AnisotropicLaplace = true;

Engine *m_ep = NULL;

int main(int argc, char** argv)
{
	//OpenMesh::IO::Options r_options, w_options; 
	//string tmeshfile = "C:\\Users\\duan_qi\\Desktop\\Reconstruction\\PhotoSynthToolkit11\\templeRing\\pmvs\\models\\pmvs_options.txt.ply";
	//string fmeshfile = "C:\\Users\\duan_qi\\Desktop\\Reconstruction\\PhotoSynthToolkit11\\templeRing\\pmvs\\models\\pmvs_options.txt.filtered.off";
	//r_options.set(OpenMesh::IO::Options::VertexColor); w_options.set(OpenMesh::IO::Options::VertexColor);
	//OpenMesh::IO::read_mesh(ObjTriMesh, tmeshfile, r_options);
	//if ( !r_options.check( OpenMesh::IO::Options::VertexColor ) ) {
	//	cout << "Color is not loaded.." << endl;
	//}
	//int rcount = 0;
	//for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
	//	OpenMesh::Vec3f tColor;
	//	tColor[0] = ObjTriMesh.color(v_it).data()[0];
	//	tColor[1] = ObjTriMesh.color(v_it).data()[1];
	//	tColor[2] = ObjTriMesh.color(v_it).data()[2];
	//	if (tColor.norm() < 40) {
	//		ObjTriMesh.delete_vertex(v_it, false);
	//		rcount ++;
	//	}
	//}
	//ObjTriMesh.garbage_collection();
	//OpenMesh::IO::write_mesh(ObjTriMesh, fmeshfile, w_options);

	ParseParam(argc,argv, MOptions);
	ScaleDelta = MOptions.ScaleDelta;
	AnisotropicLaplace = MOptions.AnisotropicLaplace;
	if (MOptions.UseMatlabSolver) {
		cout << "Use matlab solver for linear equations." << endl;
		if (!(m_ep = engOpen("\0"))) {
			std::cout << "Can not start Matlab engine" << std::endl;
			return false;
		}
		engSetVisible(m_ep, false);
	}
	ScaleDelta?cout<<"Scale delta P each time. ":cout<<" "; AnisotropicLaplace?cout<<"Using anisotropic laplacian term.":cout<<" "; cout << endl;
	path = MOptions.DirName;
	printf("Number of threads %d\n",omp_get_num_procs());
	omp_set_num_threads(omp_get_num_procs());
	omp_set_num_threads(8);
	double timer_start = (double)cv::getTickCount();
	if (!FileExisted( (MOptions.DirName + "InitialPoissonModel.ply").c_str() )) {
		if (!LoadMVSResult()) 
		{
			readMiddleBuryData2(MOptions.DirName);
			chooseStereoPairs();	

			stereoMatching();
			printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());

			buildTracks();
			printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());

			//writeToOBJ();
			//printf("\nTime = %lfs\n",((double)getTickCount()-timer_start)/getTickFrequency());

			calNormals();
			printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());

			verifyTracks();
			printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());
			
			SaveMVSResult();
		}
		writeToNPTS2(tracks,(MOptions.DirName + "PointInfo.npts"));
		//outputVerticesWithNormals(tracks,(MOptions.DirName + "PointModel.ply"));
		printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());
		//printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)w/cv::getTickFrequency());
		string PlyModelName = MOptions.DirName + "PoissonModel";
		PoissonReconstruction((MOptions.DirName + "PointInfo.npts"), PlyModelName);
		MyCopyFile(PlyModelName, (MOptions.DirName + "InitialPoissonModel.ply"));
	}

	readMiddleBuryData2(MOptions.DirName);

	//read the initial generated Poisson object model
	OpenMesh::IO::Options read_options, write_options;

	string ObjName = MOptions.DirName; 
	if (ObjName.find_last_of("\\") == ObjName.length()-1) {
		ObjName.erase(ObjName.end()-1);
	} 
	ObjName = ObjName.substr(ObjName.find_last_of("\\")+1, ObjName.length());   MOptions.meshname = ObjName;
	string meshfile = (MOptions.DirName + ObjName +"-remeshed.off");
	if (!FileExisted(meshfile.c_str())) {
		string cmd = "meshfix.exe "; cmd += MOptions.DirName + "InitialPoissonModel.ply";
		WinExec(cmd.c_str(),0);
		::Sleep(5000);
		MyMoveFile(MOptions.DirName + "InitialPoissonModel_fixed.off", MOptions.DirName + "temp.off");
		for (int i = 0; i < 5; ++ i) {
			cmd = "meshfix.exe "; cmd += MOptions.DirName + "temp.off";
			WinExec(cmd.c_str(),0);
			::Sleep(5000);
			MyMoveFile(MOptions.DirName + "temp_fixed.off", MOptions.DirName + "temp.off");
		}
		MyMoveFile(MOptions.DirName + "temp.off", meshfile);
	}
	fstream fin(meshfile,ios::in); string tag, temp_str;
	fin>>tag; fin>>temp_str;
	if (temp_str[0] == '#') {
		// need to load and rewrite the off file
		char buffer[100];
		fin.getline(buffer, 100);	fin.getline(buffer, 100);
		fstream fout(meshfile+"tmp",ios::out);
		int vnum, trinum, flag; double x, y, z; int a1, v0, v1, v2;
		fout<<tag<<endl;
		fin>>vnum>>trinum>>flag; fout<<vnum<<" "<<trinum<<" "<<flag<<endl;
		for (int i = 0; i < vnum; ++ i) {
			fin>>x>>y>>z; fout<<x<<" "<<y<<" "<<z<<endl;
		}
		for (int i = 0; i < trinum; ++i) {
			fin>>a1>>v0>>v1>>v2; fout<<a1<<" "<<v0<<" "<<v1<<" "<<v2<<endl;
		}
		fout.close();
		MyMoveFile(meshfile+"tmp", meshfile);
	} fin.close();

//	int iter_step = 0;
//	sprintf(buffer, "_%.2d", iter_step);
//	psnormalfile = (MOptions.DirName+"VertexPSNormal"+string(buffer)+".txt");
//	intensityfile = (MOptions.DirName+"VertexColor"+string(buffer)+".txt");
//	if (!FileExisted(intensityfile.c_str()) || !FileExisted(psnormalfile.c_str())) {
//		bool ok = OpenMesh::IO::read_mesh(ObjTriMesh, meshfile, read_options);
//		if (!ok) {
//			cout << "Error in load the off model " << meshfile << endl;
//		}
//		cout << "#Vertex: " << ObjTriMesh.n_vertices() << ",  #Edges: " << ObjTriMesh.n_edges() << ",  #Faces: " << ObjTriMesh.n_faces() << endl;
//		//FilterModelProcess();
//		//cout << "#Vertex after filtering: " << ObjTriMesh.n_vertices() << ",  #Edges: " << ObjTriMesh.n_edges() << ",  #Faces: " << ObjTriMesh.n_faces() << endl;
//
//		if ( !read_options.check( OpenMesh::IO::Options::VertexNormal ) ) {
//			// we need face normals to update the vertex normals
//			ObjTriMesh.request_face_normals();
//			ObjTriMesh.update_face_normals();
//			ObjTriMesh.request_vertex_normals();
//			// let the mesh update the normals
//			ObjTriMesh.update_vertex_normals();
//		}
//#ifdef USING_MATLAB
//		engEvalString(m_ep, "cd c:\\MATLAB");
//		CalculateVertexIntensity(iter_step);
//		RecoverVertexPhotometricNormal(iter_step);
//#endif
//
//#ifdef USING_NAGC
//		CubicSplineFittingPSNormal(iter_step, MOptions.fitting_choice, MOptions.range_value);
//#endif
//	}
	//CubicSplineFittingPSNormal(iter_step, MOptions.fitting_choice, MOptions.range_value);

#ifdef TEST_MESHREFINE
	// test the ps normal using original vertex normal
	// Save the output ps normal for vertex;
	meshfile = "Models\\torus-smooth.off";		MOptions.meshname = "torus";
	const char* remeshfile = "Models\\torus-remeshed.off";
	int iter_step = 0; sprintf(buffer, "_%.2d", iter_step);
	psnormalfile = (MOptions.DirName+"VertexPSNormal"+string(buffer)+".txt");
	intensityfile = (MOptions.DirName+"VertexColor"+string(buffer)+".txt");
	bool ok = OpenMesh::IO::read_mesh(ObjTriMesh, remeshfile, read_options);//bunny2-remeshed   pyramid2    Fandisk  torus  Fandisk2
	cout << "#Vertex: " << ObjTriMesh.n_vertices() << ",  #Edges: " << ObjTriMesh.n_edges() << ",  #Faces: " << ObjTriMesh.n_faces() << endl;
	ObjTriMesh.request_face_normals();
	ObjTriMesh.update_face_normals();
	ObjTriMesh.request_vertex_normals();
	ObjTriMesh.update_vertex_normals();
	fstream	fout(psnormalfile, ios::out);
	if (!fout) {
		cout << "Can not open " << psnormalfile << " to save data..." << endl; 
	}
	fout << ObjTriMesh.n_vertices() << endl; double x, y, z;
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		fout << ObjTriMesh.normal(v_it).data()[0] << "   " << ObjTriMesh.normal(v_it).data()[1] << "   " << ObjTriMesh.normal(v_it).data()[2] << endl;
	}
	fout.close();
	fout.open(intensityfile.c_str(),ios::out);
	if (!fout) {
		cout << "Can not open " << intensityfile << " to save data..." << endl; 
	}
	fout << ObjTriMesh.n_vertices() << endl; 
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		fout << 100+rand()%20 << "   " << 105+rand()%20 << endl;
	}
	fout.close();

	//store the face normal
	string facenormalfile = (MOptions.DirName+"FaceNormal"+string(buffer)+".txt");
	fout.open(facenormalfile, ios::out);
	if (!fout) {
		cout << "Can not open " << facenormalfile << " to save data..." << endl; 
	}
	fout << ObjTriMesh.n_faces() << endl; 
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		fout << ObjTriMesh.normal(f_it).data()[0] << "   " << ObjTriMesh.normal(f_it).data()[1] << "   " << ObjTriMesh.normal(f_it).data()[2] << endl;
	}
	fout.close();
	TriangularMesh testTVTM;
	testTVTM.LoadMeshFile(meshfile.c_str());
	if (!testTVTM.LoadPSNormalFile(psnormalfile.c_str()) || !testTVTM.LoadVertexColorFile(intensityfile.c_str(), MOptions.varsigma)) {
		OpenMesh::IO::read_mesh(ObjTriMesh, meshfile, read_options);
		CubicSplineFittingPSNormal(iter_step, MOptions.fitting_choice, MOptions.range_value);
		testTVTM.LoadPSNormalFile(psnormalfile.c_str());
		testTVTM.LoadVertexColorFile(intensityfile.c_str(), MOptions.varsigma);
	}
	//testTVTM.GradientTesting(0.0001, false, 0);
	//testTVTM.LM_Testing(false, 0);
	//testTVTM.LM_Testing(false, 2);
	//testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, 100, 0.0, 0.0, 0.0, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, false, false, 100);
	//testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, 0.0, 100, 0.0, 0.0, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, false, false, 101);
	//testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, 0.0, 0.0, 100, 0.0, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, false, false, 102);
	//testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, 0.0, 0.0, 0.0, 100, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, false, false, 103);
	//testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, 0.0, 0.0, 0.0, 0.0, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, true,  false, 104);
	//testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, 0.0, 0.0, 0.0, 0.0, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, false, true,  105);
	testTVTM.GradientTesting(0.0001, false, 0);
	testTVTM.GradientTesting(0.0001, false, 1);
	testTVTM.GradientTesting(0.0001, false, 2);
	//testTVTM.GradientTesting(0.0001, false, 2);
	//testTVTM.GradientTesting(0.0001, false, 3);
	//testTVTM.GradientTesting(0.0001, true, 0);
	//testTVTM.GradientTesting(0.0001, true, 1);
	testTVTM.ALM_TVU_MeshRefinement(MOptions.meshname, MOptions.fidParam, MOptions.m_beta, MOptions.ndp_eta, MOptions.ndf_eta, MOptions.varsigma, MOptions.penParam, MOptions.regParam, MOptions.lapParam, MOptions.ALM_TVU, MOptions.ALM_TVNorm, 200, MOptions.UseFaceArea, MOptions.UseMatlabSolver);
	//ReleaseResources();
	engClose(m_ep);
	return 0;
#endif

	timer_start = (double)cv::getTickCount();
	TriangularMesh TVTM;
	TVTM.LoadMeshFile(meshfile.c_str());//bunny2-smooth  pyramid2  Apple
	for (int iter_step = 0; iter_step < 3; iter_step ++) {
		cout << endl << "The " << iter_step << " iteration step of mesh refinement: " << meshfile << endl;
		sprintf(buffer, "_%.2d", iter_step);
		psnormalfile = (MOptions.DirName+"VertexPSNormal"+string(buffer)+".txt");
		intensityfile = (MOptions.DirName+"VertexColor"+string(buffer)+".txt");
		while (!TVTM.LoadPSNormalFile(psnormalfile.c_str()) || !TVTM.LoadVertexColorFile(intensityfile.c_str(), MOptions.varsigma)) {
			ObjTriMesh = TVTM.m_ObjTriMesh;
			UpdateMeshVertexIntensity(intensityfile.c_str());
			UpdateVertexPSNormal(psnormalfile.c_str());
			//CubicSplineFittingPSNormal(iter_step, MOptions.fitting_choice, MOptions.range_value);
		}
		TVTM.ALM_TVU_MeshRefinement(MOptions.meshname, MOptions.fidParam, MOptions.m_beta, MOptions.ndp_eta, MOptions.ndf_eta, MOptions.varsigma, MOptions.pcParam,
			MOptions.penParam, MOptions.regParam, MOptions.lapParam, MOptions.ALM_TVU, MOptions.ALM_TVNorm, iter_step, MOptions.UseFaceArea, MOptions.UseMatlabSolver);
	}
	
	printf("\nALM operation time = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());

	
	printf("\n\nClearing memory...\n");
	//ReleaseResources();
	printf("\n\nFinished!\n");
	printf("\nTime = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());

	//_getch();
	if (MOptions.UseMatlabSolver) {
		engClose(m_ep);
	}
	return 0;
}