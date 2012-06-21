#include "MeshTVRefine.h"
extern Engine *m_ep;

MeshTVRefine::MeshTVRefine()
{
	ObjTriMesh.clear();
	Vertex_Intensity.clear();
	Vertex_PSNormal.clear();
	m_BCDArea.clear();
	m_TriangleArea.clear();
	m_divergence.clear();
	m_vertexweight.clear();
	m_vertexfaces.clear();
	m_edgelaplaceweight.clear();
	m_TPPIBG.clear();
	m_TriGradient.clear();
	//m_alpha = 500;
	//m_beta = 100;
	//m_gamma = 0.01;
	//m_sigma = 0.5;
}

MeshTVRefine::~MeshTVRefine()
{
	ObjTriMesh.clear();
	Vertex_Intensity.clear();
	Vertex_PSNormal.clear();
	m_BCDArea.clear();
	m_TriangleArea.clear();
	m_divergence.clear();
	m_vertexweight.clear();
	m_vertexfaces.clear();
	m_edgelaplaceweight.clear();
	m_TPPIBG.clear();
	m_TriGradient.clear();
}

MeshTVRefine& MeshTVRefine::operator=(const MeshTVRefine &inmtvr)
{
	if (this == &inmtvr) {
		return *this;
	}
	ObjTriMesh.clear();			ObjTriMesh = inmtvr.ObjTriMesh;
	Vertex_Intensity.clear();	Vertex_Intensity = inmtvr.Vertex_Intensity;
	Vertex_PSNormal.clear();	Vertex_PSNormal = inmtvr.Vertex_PSNormal;
	m_BCDArea.clear();			m_BCDArea = inmtvr.m_BCDArea;
	m_TriangleArea.clear();		m_TriangleArea = inmtvr.m_TriangleArea;
	m_divergence.clear();		m_divergence = inmtvr.m_divergence;
	m_vertexweight.clear();		m_vertexweight = inmtvr.m_vertexweight;
	m_vertexfaces.clear();		m_vertexfaces = inmtvr.m_vertexfaces;
	m_edgelaplaceweight.clear();m_edgelaplaceweight = inmtvr.m_edgelaplaceweight;
	m_TPPIBG.clear();			m_TPPIBG = inmtvr.m_TPPIBG;
	m_TriGradient.clear();		m_TriGradient = inmtvr.m_TriGradient;
	//m_sigma = inmtvr.m_sigma;
	//m_alpha = inmtvr.m_alpha;
	//m_beta = inmtvr.m_beta;
	//m_gamma = inmtvr.m_sigma;
	return *this;
}

MeshTVRefine::MeshTVRefine(const MeshTVRefine &inmtvr)
{
	ObjTriMesh.clear();			ObjTriMesh = inmtvr.ObjTriMesh;
	Vertex_Intensity.clear();	Vertex_Intensity = inmtvr.Vertex_Intensity;
	Vertex_PSNormal.clear();	Vertex_PSNormal = inmtvr.Vertex_PSNormal;
	m_BCDArea.clear();			m_BCDArea = inmtvr.m_BCDArea;
	m_TriangleArea.clear();		m_TriangleArea = inmtvr.m_TriangleArea;
	m_divergence.clear();		m_divergence = inmtvr.m_divergence;
	m_vertexweight.clear();		m_vertexweight = inmtvr.m_vertexweight;
	m_vertexfaces.clear();		m_vertexfaces = inmtvr.m_vertexfaces;
	m_edgelaplaceweight.clear();m_edgelaplaceweight = inmtvr.m_edgelaplaceweight;
	m_TPPIBG.clear();			m_TPPIBG = inmtvr.m_TPPIBG;
	m_TriGradient.clear();		m_TriGradient = inmtvr.m_TriGradient;
	//m_sigma = inmtvr.m_sigma;
	//m_alpha = inmtvr.m_alpha;
	//m_beta = inmtvr.m_beta;
	//m_gamma = inmtvr.m_sigma;
}

bool MeshTVRefine::LoadMeshFile(const char* filename)
{
	//numc::RowMat<double> AlphaNormVLap(5,5); double b[5] = {2,3,4,5,6}; double x[5];
	//AlphaNormVLap(0,0) = 1.0;	AlphaNormVLap(0,1) = -1.0;	AlphaNormVLap(0,3) = -3.0;
	//AlphaNormVLap(1,0) = -2.0;	AlphaNormVLap(1,1) = 5.0;	
	//AlphaNormVLap(2,2) = 4.0;	AlphaNormVLap(2,3) = 6.0;	AlphaNormVLap(2,4) = 4.0;
	//AlphaNormVLap(3,0) = -4.0;	AlphaNormVLap(3,2) = 2.0;	AlphaNormVLap(3,3) = 7.0;
	//AlphaNormVLap(4,1) = 8.0;	AlphaNormVLap(4,4) = -5.0;
	//numc::SparseSolver solver;
	//solver.getMatA() = AlphaNormVLap;
	//solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
	//solver.init();
	//solver.solve(b,x);
	//cout << "MKL result: "<< endl;
	//for (int i = 0; i < 5; ++ i) {
	//	cout << x[i] << endl;
	//}
	//// test the matlab sparse matrix
	//numc::CSRMatrix<double> testat = AlphaNormVLap;
	//numc::CSRMatrix<double> testa = AlphaNormVLap;
	////CSRMatrixTranspose(testat, testa);
	////testa.ChangeBase(false);
	//double prdata[100]; double irdata[100]; double jcdata[100]; 
	//for (int i = 0; i < testa.mAv.size(); ++ i) {
	//	prdata[i] = testa.mAv[i];
	//	irdata[i] = testa.mAj[i];
	//}
	//for (int i = 0; i < testa.mAi.size(); ++ i) {
	//	jcdata[i] = testa.mAi[i];
	//}
	//mxArray* spmatrix = mxCreateSparse(5,5,15,mxREAL);
	//double* start_of_pr = (double*)mxGetPr(spmatrix);
	//memcpy(start_of_pr, a, 13*sizeof(double));
	//int* start_of_ir = (int *)mxGetIr(spmatrix); 
	//memcpy(start_of_ir, jc, 13*sizeof(int));
	//int* start_of_jc = (int *)mxGetJc(spmatrix); 
	//memcpy(start_of_jc, ia, 6*sizeof(int));
	//mxArray* barray = mxCreateDoubleMatrix(5,1,mxREAL);
	//double* bptr = mxGetPr(barray);
	//bptr[0] = 1; bptr[1] = 2; bptr[2] = 3; bptr[3] = 4; bptr[4] = 5;
	//engPutVariable(m_ep, "SB", barray);
	//engPutVariable(m_ep, "SA", spmatrix);
	//engEvalString(m_ep, "SX=SA\\SB");
	//mxArray* xarray = engGetVariable(m_ep, "SX");
	//double* xprt = mxGetPr(xarray);
	//cout << "Matlab result: "<< endl;
	//for (int i = 0; i < 5; ++ i) {
	//	cout << xprt[i] << endl;
	//}
	OpenMesh::IO::Options read_options;
	bool ok = OpenMesh::IO::read_mesh(this->ObjTriMesh,string(filename), read_options);
	if (!ok) {
		cout << "Error in load the off model " << filename << endl;
		return false;
	}
	cout << "#Vertex: " << ObjTriMesh.n_vertices() << ",  #Edges: " << ObjTriMesh.n_edges() << ",  #Faces: " << ObjTriMesh.n_faces() << endl;

	if ( !read_options.check( OpenMesh::IO::Options::VertexNormal ) ) {
		// we need face normals to update the vertex normals
		ObjTriMesh.request_face_normals();
		// let the mesh update the normals
		ObjTriMesh.update_normals();
		// dispose the face normals, as we don't need them anymore
		ObjTriMesh.release_face_normals();
	}
	cout << "Load mesh file: " << filename << " Done..." << endl;
	return true;
}

bool MeshTVRefine::LoadColorFile(const char* filename, double m_sigma)
{
	fstream fin(filename, ios::in);
	if (!fin) {
		cout << "Can not open " << filename<< " to read color data..." << endl;
		return false;
	}
	int c_size; fin>>c_size;
	if (c_size != ObjTriMesh.n_vertices()&&ObjTriMesh.n_vertices()>0) {
		cout << "The size of color data " << c_size << " is not match with mesh vertex number " << ObjTriMesh.n_vertices() << endl;
		return false;
	}
	Vertex_Intensity.clear(); Vertex_Intensity.resize(c_size);
	double original_intensity, fitted_intensity;
	for (int i = 0; i < c_size; i ++) {
		 fin >> fitted_intensity >> original_intensity;
#ifdef USE_FITTED_COLOR
		 original_intensity = fitted_intensity;
#endif
		 Vertex_Intensity[i] = original_intensity;
	}
	fin.close();
	cout << "Load color file: " << filename << " Done..." << endl;
	this->UpdateVertexWeight(m_sigma);
	return true;
}

bool MeshTVRefine::LoadPSNormalFile(const char* filename)
{
	fstream fin(filename, ios::in);
	if (!fin) {
		cout << "Can not open " << filename<< " to read photometric normal data..." << endl;
		return false;
	}
	int c_size; fin>>c_size;
	if (c_size != ObjTriMesh.n_vertices()&&ObjTriMesh.n_vertices()>0) {
		cout << "The size of PS normal data "<< c_size <<" is not match with mesh vertex number " << ObjTriMesh.n_vertices() << endl;
		return false;
	}
	Vertex_PSNormal.clear(); Vertex_PSNormal.resize(c_size);	OpenMesh::Vec3f v_normal;
	for (int i = 0; i < c_size; i ++) {
		fin >> v_normal[0] >> v_normal[1] >> v_normal[2]; 
		Vertex_PSNormal[i] = v_normal.normalize();
	}
	fin.close();
	cout << "Load PS normal file: " << filename << " Done..." << endl;
	return true;
}

bool MeshTVRefine::LoadFaceNormalFile(const char* filename)
{
	fstream fin(filename, ios::in);
	if (!fin) {
		cout << "Can not open " << filename << " to read face normal data..." << endl;
		return false;
	}
	int c_size; fin>>c_size;
	if (c_size != ObjTriMesh.n_faces()&&ObjTriMesh.n_faces()>0) {
		cout << "The size of face normal data "<< c_size <<" is not match with mesh vertex number " << ObjTriMesh.n_faces() << endl;
		return false;
	}
	Face_Normal.clear(); Face_Normal.resize(c_size);	OpenMesh::Vec3f v_normal;
	for (int i = 0; i < c_size; i ++) {
		fin >> v_normal[0] >> v_normal[1] >> v_normal[2]; 
		Face_Normal[i] = v_normal.normalize();
	}
	fin.close();
	cout << "Load face normal file: " << filename << " Done..." << endl;
	return true;
}

void MeshTVRefine::UpdateBCDArea()
{
	//calculate the area of the bc dual for each vertex
	if (m_BCDArea.size() != ObjTriMesh.n_vertices()) {
		m_BCDArea.clear();		m_BCDArea.resize(ObjTriMesh.n_vertices());
		m_vertexfaces.clear();	m_vertexfaces.resize(ObjTriMesh.n_vertices());// store the num of faces with certain vertex
	}
	OpenMesh::Vec3f pointA , pointB , pointC, pointThis, tri_barycenter; int ida, idb, idc;
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		pointThis = ObjTriMesh.point(v_it.handle());		
		double bcarea = 0.0; 	int num_faces = 0;
		for (MyMesh::VertexFaceIter vf_it = ObjTriMesh.vf_iter(v_it.handle()); vf_it; ++ vf_it) {
			MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(vf_it.handle());
			pointA = ObjTriMesh.point(cfvIt.handle());			ida = cfvIt.handle().idx();		++cfvIt;
			pointB = ObjTriMesh.point(cfvIt.handle());			idb = cfvIt.handle().idx();		++cfvIt;
			pointC = ObjTriMesh.point(cfvIt.handle());			idc = cfvIt.handle().idx();
			tri_barycenter = (pointA+pointB+pointC)/3.0;
			if (ida == v_it.handle().idx()) {
				bcarea += OpenMesh::cross((pointB - pointThis)/2.0, (tri_barycenter - pointThis)).norm()/2.0;
				bcarea += OpenMesh::cross((pointC - pointThis)/2.0, (tri_barycenter - pointThis)).norm()/2.0;
			}
			if (idb == v_it.handle().idx()) {
				bcarea += OpenMesh::cross((pointA - pointThis)/2.0, (tri_barycenter - pointThis)).norm()/2.0;
				bcarea += OpenMesh::cross((pointC - pointThis)/2.0, (tri_barycenter - pointThis)).norm()/2.0;
			}
			if (idc == v_it.handle().idx()) {
				bcarea += OpenMesh::cross((pointB - pointThis)/2.0, (tri_barycenter - pointThis)).norm()/2.0;
				bcarea += OpenMesh::cross((pointA - pointThis)/2.0, (tri_barycenter - pointThis)).norm()/2.0;
			}
			num_faces ++;
		}
		m_BCDArea[v_it.handle().idx()] = bcarea;	m_vertexfaces[v_it.handle().idx()] = num_faces;
	}
	cout << "Update mesh BCD Area Done..." << endl;
#ifdef _DEBUG
	fstream of("VertexBCDArea.txt",std::ios::out);
	if (!of) {
		cout << "Can not open VertexBCDArea.txt to save data..." << endl;
		return;
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		of << v_it.handle().idx() << ": " << m_BCDArea[v_it.handle().idx()] << endl;
	}
	of.close();
#endif
#ifdef _DEBUG
	fstream of1("VertexFaceNumber.txt",std::ios::out);
	if (!of1) {
		cout << "Can not open VertexFaceNumber.txt to save data..." << endl;
		return;
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		of1 << v_it.handle().idx() << ": " << m_vertexfaces[v_it.handle().idx()] << endl;
	}
	of1.close();
#endif
	return;
}

void MeshTVRefine::UpdateTriangleArea()
{
	m_TriangleArea.clear(); m_TriangleArea.resize(ObjTriMesh.n_faces());
	OpenMesh::Vec3f pointA , pointB , pointC;
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(f_it.handle());
		pointA = ObjTriMesh.point(cfvIt.handle());	++cfvIt;
		pointB = ObjTriMesh.point(cfvIt.handle());	++cfvIt;
		pointC = ObjTriMesh.point(cfvIt.handle());
		m_TriangleArea[f_it.handle().idx()] = OpenMesh::cross((pointB-pointA), (pointC-pointA)).norm()/2.0;
	}
	cout << "Update mesh triangular face area Done..." << endl;
#ifdef _DEBUG
	fstream of("TrianguleArea.txt",std::ios::out);
	if (!of) {
		cout << "Can not open TrianguleArea.txt to save data..." << endl;
		return;
	}
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		of << f_it.handle().idx() << ": " << m_TriangleArea[f_it.handle().idx()] << endl;
	}
	of.close();
#endif
	return;
}

bool MeshTVRefine::UpdateVertexWeight(double m_sigma)
{
	//calculate the vertex weight for tv term
	if (Vertex_Intensity.size() != ObjTriMesh.n_vertices()) {
		cout << "The intensity size is not equal to vertex number." << endl;
		return false;
	}
	double max_weight = -100;
	if (m_vertexweight.size() != ObjTriMesh.n_vertices()) {
		m_vertexweight.clear(); m_vertexweight.resize(ObjTriMesh.n_vertices()); 
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		int count = 0; int avg_intensity = 0.0;
		for (MyMesh::VertexVertexIter vv_it = ObjTriMesh.vv_iter(v_it.handle()); vv_it; ++vv_it) {
			count ++;
			avg_intensity += Vertex_Intensity[vv_it.handle().idx()];
		}
		avg_intensity = avg_intensity/count;
		double cur_weight = std::exp(-std::pow((Vertex_Intensity[v_it.handle().idx()] - avg_intensity),2.0)/(m_sigma*m_sigma) );
		m_vertexweight[v_it.handle().idx()] = cur_weight;
		if (cur_weight > max_weight) {
			max_weight = cur_weight;
		}
	}
	for (int i = 0; i < m_vertexweight.size(); ++ i) {
		m_vertexweight[i] = m_vertexweight[i]/max_weight;
	}
	cout << "Update mesh vertex weight for PS normal term Done..." << endl;
#ifdef _DEBUG
	fstream of("VertexWeight.txt",std::ios::out);
	if (!of) {
		cout << "Can not open VertexWeight.txt to save data..." << endl;
		return false;
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		of << v_it.handle().idx() << ": " << m_vertexweight[v_it.handle().idx()] << endl;
	}
	of.close();
#endif
	return true;
}

float MeshTVRefine::CalculateHalfEdgeLaplaceWeight(MyMesh::HalfedgeHandle hehd)
{
	OpenMesh::Vec3f pp, pq, pr, v1, v2; float result;
	MyMesh::HalfedgeHandle prev_hehd, next_hehd;
	prev_hehd = ObjTriMesh.prev_halfedge_handle(hehd);
	next_hehd = ObjTriMesh.next_halfedge_handle(hehd);
	pp = ObjTriMesh.point(ObjTriMesh.to_vertex_handle(hehd));
	pq = ObjTriMesh.point(ObjTriMesh.to_vertex_handle(prev_hehd));
	pr = ObjTriMesh.point(ObjTriMesh.to_vertex_handle(next_hehd));
	v1 = pp - pr; v2 = pq - pr;
	double cosan = OpenMesh::dot(v1, v2)/(v1.norm()*v2.norm());
	double angle = std::acos(cosan);
	if (abs(std::tan(angle))<0.00001)
	{
		result = 10000;
		printf(" 0 angle!!!!!!\n");
	}
	else
	{
		result = 1.0/(tan(angle));
	}
	return result;
}

void MeshTVRefine::UpdateEdgeLaplaceWeight()
{ 
	if (m_edgelaplaceweight.size() != ObjTriMesh.n_edges()) {
		m_edgelaplaceweight.clear(); m_edgelaplaceweight.resize(ObjTriMesh.n_edges());
	}
	for (MyMesh::EdgeIter e_it = ObjTriMesh.edges_begin(); e_it != ObjTriMesh.edges_end(); ++ e_it) {
		double result = 0.0; int ncount = 0;
		MyMesh::HalfedgeHandle hehd = ObjTriMesh.halfedge_handle(e_it.handle(), 0);
		if (!ObjTriMesh.is_boundary(hehd)) {
			result += CalculateHalfEdgeLaplaceWeight(hehd);
			ncount ++;
		}
		MyMesh::HalfedgeHandle ohehd = ObjTriMesh.halfedge_handle(e_it.handle(), 1); // = ObjTriMesh.opposite_halfedge_handle(hehd);
		if (!ObjTriMesh.is_boundary(ohehd)) {
			result += CalculateHalfEdgeLaplaceWeight(ohehd);
			ncount ++;
		}
		if (ncount < 2) {
			cout << e_it.handle().idx() << " has only one adjacent face." << endl;
		}
		m_edgelaplaceweight[e_it.handle().idx()] = result/ncount;
	}
	cout << "Update mesh edge laplacian weight Done..." << endl;
#ifdef _DEBUG
	fstream of("EdgeLaplaceWeight.txt",std::ios::out);
	if (!of) {
		cout << "Can not open EdgeLaplaceWeight.txt to save data..." << endl;
	}
	for (MyMesh::EdgeIter e_it = ObjTriMesh.edges_begin(); e_it != ObjTriMesh.edges_end(); ++ e_it) {
		of << e_it.handle().idx() << ": " << m_edgelaplaceweight[e_it.handle().idx()] << endl;
	}
	of.close();
#endif
	return;
}

void MeshTVRefine::FindPPIBGradient(MyMesh &InMesh, const int ver_id, const MyMesh::FaceHandle f_handle, OpenMesh::Vec3f &vec_g, float &h)
{
	//calculate the gradient of the primal-primal interpolation function of certain point in certain triangle
	double alpha; OpenMesh::Vec3f v1, v2, pf;
	OpenMesh::Vec3f pointA , pointB , pointC; int ida, idb, idc;
	MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(f_handle);
	pointA = ObjTriMesh.point(cfvIt.handle());		ida = cfvIt.handle().idx(); ++cfvIt;
	//cout << "PointA: " << pointA[0] << ", " << pointA[1] << ", " << pointA[2] << endl;
	pointB = ObjTriMesh.point(cfvIt.handle());		idb = cfvIt.handle().idx(); ++cfvIt;
	//cout << "PointB: " << pointB[0] << ", " << pointB[1] << ", " << pointB[2] << endl;
	pointC = ObjTriMesh.point(cfvIt.handle());		idc = cfvIt.handle().idx();
	//cout << "PointC: " << pointC[0] << ", " << pointC[1] << ", " << pointC[2] << endl;
	if (ver_id == ida) {
		v1 = pointB - pointC;  v2 = pointA - pointC;
		//cout << "v1: " << v1[0] << ", " << v1[1] << ", " << v1[2] << endl;
		//cout << "v2: " << v2[0] << ", " << v2[1] << ", " << v2[2] << endl;
		alpha = OpenMesh::dot(v2, v1)/OpenMesh::dot(v1, v1);
		pf = pointB*alpha + pointC*(1-alpha);
		//cout << "pf: " << pf[0] << ", " << pf[1] << ", " << pf[2] << endl;
		vec_g = pointA - pf;
		//cout << "vec_g: " << vec_g[0] << ", " << vec_g[1] << ", " << vec_g[2] << endl;
		h = vec_g.norm();
		vec_g = vec_g/(h*h);
		//cout << "vec_g: " << vec_g[0] << ", " << vec_g[1] << ", " << vec_g[2] << endl;
	} else if (ver_id == idb) {
		v1 = pointC - pointA;  v2 = pointB - pointA;
		alpha = OpenMesh::dot(v2, v1)/OpenMesh::dot(v1, v1);
		pf = pointC*alpha + pointA*(1-alpha);
		vec_g = pointB - pf;
		h = vec_g.norm();
		vec_g = vec_g/(h*h);
	} else {
		v1 = pointA - pointB;  v2 = pointC - pointB;
		alpha = OpenMesh::dot(v2, v1)/OpenMesh::dot(v1, v1);
		pf = pointA*alpha + pointB*(1-alpha);
		vec_g = pointC - pf;
		h = vec_g.norm();
		vec_g = vec_g/(h*h);
	}
	return;
}

void MeshTVRefine::BuildTPPIBG(MyMesh &InMesh, vector< map<int, OpenMesh::Vec3f> > &vec_TPPIBG)
{ 
	if (vec_TPPIBG.size() != InMesh.n_faces()) {
		vec_TPPIBG.clear(); vec_TPPIBG.resize(InMesh.n_faces());
	}
	map<int, OpenMesh::Vec3f> t_TPPIBG; OpenMesh::Vec3f t_grad;  float h;
	for (MyMesh::FaceIter f_it = InMesh.faces_begin(); f_it != InMesh.faces_end(); ++ f_it) {
		t_TPPIBG.clear();
		for (MyMesh::ConstFaceVertexIter cfv_it = InMesh.cfv_iter(f_it.handle()); cfv_it; ++ cfv_it) {
			FindPPIBGradient(InMesh, cfv_it.handle().idx(), f_it.handle(), t_grad, h);
			t_TPPIBG.insert(pair<int, OpenMesh::Vec3f>(cfv_it.handle().idx(), t_grad));
		}
		vec_TPPIBG[f_it.handle().idx()] = t_TPPIBG;
	}
	cout << "Build mesh TPPIB gradient for each face Done..." << endl;
#ifdef _DEBUG
	fstream of("TPPIBG.txt",std::ios::out);
	if (!of) {
		cout << "Can not open TPPIBG.txt to save data..." << endl;
		return;
	}
	for (MyMesh::FaceIter f_it = InMesh.faces_begin(); f_it != InMesh.faces_end(); ++ f_it) {
		t_TPPIBG = m_TPPIBG[f_it.handle().idx()]; 
		of << f_it.handle().idx() << ":  ";
		for (map<int, OpenMesh::Vec3f>::iterator iter = t_TPPIBG.begin(); iter != t_TPPIBG.end(); iter ++) {
			of << iter->first << ", " << iter->second[0] << ", " << iter->second[1] << ", " << iter->second[2] << "; ";
		}
		of << endl;
	}
	of.close();
#endif
	return;
}

bool MeshTVRefine::CalDivengence(const vector< vector<OpenMesh::Vec3f> > &vf)
{
	//calculate the divergence for each vertex using the the vector field vf defined in each face.
	// the divergence result should divide the bcdarea 
	if (vf.size() != ObjTriMesh.n_faces()) {
		cout << "The size of vector field does not equal to the number of triangular faces..." << endl;
		return false;
	}
	if (m_divergence.size() != ObjTriMesh.n_vertices()) {
		m_divergence.clear(); m_divergence.resize(ObjTriMesh.n_vertices());
	}
	double div_value;	map<int, OpenMesh::Vec3f> t_TPPIBG;
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		vector <double> t_divergence; t_divergence.clear();
		for (int i = 0; i < 3; ++ i) {	// for the x, y, z channels
			div_value = 0.0;
			for (MyMesh::VertexFaceIter vf_it = ObjTriMesh.vf_iter(v_it.handle()); vf_it; ++ vf_it) {
				t_TPPIBG = m_TPPIBG[vf_it.handle().idx()];
				for (map<int, OpenMesh::Vec3f>::iterator iter = t_TPPIBG.begin(); iter != t_TPPIBG.end(); iter ++) {
					if (v_it.handle().idx() == iter->first) {
						div_value -= m_TriangleArea[vf_it.handle().idx()]*OpenMesh::dot(vf[vf_it.handle().idx()][i], iter->second);
					}
				}
			}
			t_divergence.push_back(div_value/m_BCDArea[v_it.handle().idx()]);
		}
		m_divergence[v_it.handle().idx()] = t_divergence;
	}
	return true;
}

bool MeshTVRefine::CalGradient(const vector<OpenMesh::Vec3f> &vertexvec, vector< vector<OpenMesh::Vec3f> > &out_TriGradient)
{
	// Calculate the \nabla v for each triangle in x, y, z channel
	if (vertexvec.size() != ObjTriMesh.n_vertices()) {
		cout << "The size of vertex number does not equal to the number of mesh vertices..." << endl;
		return false;
	}
	if (out_TriGradient.size() != ObjTriMesh.n_faces()) {
		out_TriGradient.clear(); out_TriGradient.resize(ObjTriMesh.n_faces());
	}
	map<int, OpenMesh::Vec3f> t_TPPIBG; OpenMesh::Vec3f t_gradient;
	vector<int> idx_list; vector<OpenMesh::Vec3f> tppibg_list; 
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		t_TPPIBG = m_TPPIBG[f_it.handle().idx()];
		if (t_TPPIBG.size() != 3) {
			cout << "The size of TPPIG for face " << f_it.handle().idx() << " not equal to 3..." << endl;
			return false;
		}
		idx_list.clear();  tppibg_list.clear();
		for (map<int, OpenMesh::Vec3f>::iterator iter = t_TPPIBG.begin(); iter != t_TPPIBG.end(); iter ++) {
			idx_list.push_back(iter->first);	tppibg_list.push_back(iter->second);
		}
		vector<OpenMesh::Vec3f> vec_gradient; vec_gradient.clear();
		for (int i = 0; i < 3; ++ i) {
			t_gradient[0] = vertexvec[idx_list[0]][i]*tppibg_list[0][0] + vertexvec[idx_list[1]][i]*tppibg_list[1][0] + vertexvec[idx_list[2]][i]*tppibg_list[2][0];
			t_gradient[1] = vertexvec[idx_list[0]][i]*tppibg_list[0][1] + vertexvec[idx_list[1]][i]*tppibg_list[1][1] + vertexvec[idx_list[2]][i]*tppibg_list[2][1];
			t_gradient[2] = vertexvec[idx_list[0]][i]*tppibg_list[0][2] + vertexvec[idx_list[1]][i]*tppibg_list[1][2] + vertexvec[idx_list[2]][i]*tppibg_list[2][2];
			vec_gradient.push_back(t_gradient);
		}
		out_TriGradient[f_it.handle().idx()] = vec_gradient;
	}
}

void MeshTVRefine::LocalNormalRefinement(const double ialpha, const double ibeta, const double igamma)
{
	// refine the fitted surface normal locally using the neighbor points normal
	// min \sum_{i}(a/2 * (n_i - n_i')^2 + b/2 \sum_{j\in N(i)} (n_i - n_j)^2 )
	// need to load LoadColorFile and LoadPSNormalFile first, then refine the PS normal
	double *bx = new double[ObjTriMesh.n_vertices()]; 
	double *nx = new double[ObjTriMesh.n_vertices()]; 
	numc::RowMat<double> AlphaNormVLap; AlphaNormVLap.resize(ObjTriMesh.n_vertices(), ObjTriMesh.n_vertices()); AlphaNormVLap *= 0.0;
	numc::SparseSolver solver;
	cout << "Start the local mesh normal refinement process: ";
	// construct the coefficient matrix for local mesh normal refinement
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		double totalweight = 0.0; int totalcount = 0;
		for (MyMesh::ConstVertexVertexIter cvv_it = ObjTriMesh.cvv_iter(v_it.handle()); cvv_it; ++ cvv_it) {
			double ijweight = std::exp(-std::pow((Vertex_Intensity[v_it.handle().idx()] - Vertex_Intensity[cvv_it.handle().idx()]), 2.0)/igamma);
			totalweight += ijweight; totalcount ++;
		}
		for (MyMesh::ConstVertexVertexIter cvv_it = ObjTriMesh.cvv_iter(v_it.handle()); cvv_it; ++ cvv_it) {
			double ijweight = std::exp(-std::pow((Vertex_Intensity[v_it.handle().idx()] - Vertex_Intensity[cvv_it.handle().idx()]), 2.0)/igamma);
			AlphaNormVLap(v_it.handle().idx(), cvv_it.handle().idx()) = -ibeta*ijweight/totalweight;
		}
		AlphaNormVLap(v_it.handle().idx(), v_it.handle().idx()) = ialpha + ibeta;
	}
	solver.getMatA() = AlphaNormVLap;
	solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
	solver.init();

	vector<OpenMesh::Vec3f> TPSNormal; TPSNormal.resize(ObjTriMesh.n_vertices());
	for (int i = 0; i < 3; ++ i) {
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			bx[v_it.handle().idx()]	= ialpha*Vertex_PSNormal[v_it.handle().idx()][i];	// for the x,y,z channel
		}
		solver.solve(bx, nx);
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			TPSNormal[v_it.handle().idx()][i] = nx[v_it.handle().idx()];
		}
	}
	
	// solve the nonlinear equation to calculate the new vertex position in vx;
	solver.clear();
//#ifdef _DEBUG
	fstream of("localrefinedpsnormal.txt",std::ios::out);
	if (!of) {
		cout << "Can not open localrefinedpsnormal.txt to save data..." << endl;
		return;
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		TPSNormal[v_it.handle().idx()].normalize();
		of << TPSNormal[v_it.handle().idx()][0] << "  " << TPSNormal[v_it.handle().idx()][1] << "  " << TPSNormal[v_it.handle().idx()][2] << "  " 
			<< Vertex_PSNormal[v_it.handle().idx()][0] << "  " << Vertex_PSNormal[v_it.handle().idx()][1] << "  " << Vertex_PSNormal[v_it.handle().idx()][2] << endl;
	}
	of.close();
//#endif
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		Vertex_PSNormal[v_it.handle().idx()] = TPSNormal[v_it.handle().idx()].normalize(); 
	}
	delete nx, bx; AlphaNormVLap.clear();
	cout << "Done..." << endl; 
	return;
}

void MeshTVRefine::MeshRefineByNormal(double ialpha, double ibeta)
{
	this->UpdateTriangleArea();
	this->UpdateBCDArea();
	// refine the mesh by normal direction, using least square method
	double *b = new double[ObjTriMesh.n_vertices()*3]; 
	double *vx = new double[ObjTriMesh.n_vertices()*3]; 
	numc::RowMat<double> AlphaNormVLap; AlphaNormVLap.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::SparseSolver solver;
	// Initialization process
	vector < OpenMesh::Vec3f > oldV, curV; oldV.clear(); curV.clear();
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		oldV.push_back(ObjTriMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
		curV.push_back(ObjTriMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
	}
	double StopCond = 0.00000;  int innerL = 1; int outL = 0;
	double outTole = 5e-10;
	OpenMesh::Vec3f pointA , pointB , pointC; int ida, idb, idc;
	cout << "Start the mesh refinement by normal process: ";
	do {
		outL ++;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			oldV[v_it.handle().idx()] = curV[v_it.handle().idx()];
		}
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			b[v_it.handle().idx()] =							ialpha*ObjTriMesh.point(v_it.handle())[0];	// for the x channel
			b[v_it.handle().idx()+ObjTriMesh.n_vertices()] =	ialpha*ObjTriMesh.point(v_it.handle())[1];	// for the y channel
			b[v_it.handle().idx()+ObjTriMesh.n_vertices()*2] =	ialpha*ObjTriMesh.point(v_it.handle())[2];	// for the z channel
		}
		// Construct the nonlinear matrix A
		//cout << AlphaNormVLap.clearZero() << " non-zero elements are erased." << endl;
		AlphaNormVLap *= 0.0;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
			// fidelity term 
			AlphaNormVLap(vertex_id, vertex_id) =														ialpha; // for x channel
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) =		ialpha; // for y channel
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) =	ialpha; // for z channel
			//end of fidelity term

			// normal term
			// first for all neighbor triangle, store its vertex in counterclock wise, A(current vi), B, C
			vector< vector<OpenMesh::Vec3f> >	NeiFaceVertexList;		NeiFaceVertexList.clear();
			vector< vector<int> >				NeiFaceVertexIDList;	NeiFaceVertexIDList.clear();
			vector< double >					NeiFaceAreaList;		NeiFaceAreaList.clear();
			double								SumNeiFaceArea = 0.0000;
			for (MyMesh::VertexFaceIter vf_it = ObjTriMesh.vf_iter(v_it.handle()); vf_it; ++ vf_it) {
				vector<OpenMesh::Vec3f> FaceVertexList;		FaceVertexList.clear();
				vector<int>				FaceVertexIDList;	FaceVertexIDList.clear();
				MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(vf_it.handle());
				pointA = curV[cfvIt.handle().idx()];			ida = cfvIt.handle().idx();
				pointB = curV[(++cfvIt).handle().idx()];		idb = cfvIt.handle().idx();
				pointC = curV[(++cfvIt).handle().idx()];		idc = cfvIt.handle().idx();
				if (vertex_id == ida) {
					//pointA v_i; pointB v_i+1; pointC v_i+2
					FaceVertexList.push_back(pointA); FaceVertexList.push_back(pointB); FaceVertexList.push_back(pointC);
					FaceVertexIDList.push_back(ida); FaceVertexIDList.push_back(idb); FaceVertexIDList.push_back(idc);
				} else if (vertex_id == idb){ 
					FaceVertexList.push_back(pointB); FaceVertexList.push_back(pointC); FaceVertexList.push_back(pointA);
					FaceVertexIDList.push_back(idb); FaceVertexIDList.push_back(idc); FaceVertexIDList.push_back(ida);
				} else if (vertex_id == idc) {
					FaceVertexList.push_back(pointC); FaceVertexList.push_back(pointA); FaceVertexList.push_back(pointB);
					FaceVertexIDList.push_back(idc); FaceVertexIDList.push_back(ida); FaceVertexIDList.push_back(idb);
				} else {
					cout << "Something wrong happened in the OpenMesh iteration ..." << endl;
				}
				NeiFaceVertexList.push_back(FaceVertexList);
				NeiFaceVertexIDList.push_back(FaceVertexIDList);
				NeiFaceAreaList.push_back(OpenMesh::cross((pointC-pointB), (pointA-pointB)).norm());
				SumNeiFaceArea += OpenMesh::cross((pointC-pointB), (pointA-pointB)).norm();
			}
			// then the norm of each neighbor face is a 3*1 vector Nkj = (C-B)*(A-B), its derivative to a 3*3 jacob matrix Mki; then Mki*Nkj
			// actually Mki is a known matrix
			AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()) = 0.0;
			AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) = 0.0; // set the initial value for v_i.y and v_i.z; v_i.x is already set
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id) = 0.0;
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) = 0.0;
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
			float KNum = (float)NeiFaceVertexList.size();  //double BArea = m_BCDArea[vertex_id];
			double TwoSumNeiFaceArea = SumNeiFaceArea * 2;
			OpenMesh::Vec3f tpa, tpb, tpc;
			for (int i = 0; i < NeiFaceVertexList.size(); ++ i) {
				tpa = NeiFaceVertexList[i][0];	tpb = NeiFaceVertexList[i][1];	tpc = NeiFaceVertexList[i][2]; // in Mki
				for (int j = 0; j < NeiFaceVertexList.size(); ++ j) {
					// calculate the 36 term for the normal constraint
					pointA = NeiFaceVertexList[j][0];	pointB = NeiFaceVertexList[j][1];	pointC = NeiFaceVertexList[j][2]; // in Nkj
					// for x channel
					AlphaNormVLap(vertex_id, vertex_id) += 
						((tpb[2]-tpc[2])*(pointC[2]-pointB[2]) - (tpc[1]-tpb[1])*(pointC[1]-pointB[1]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()) += 
						(tpc[1]-tpb[1])*(pointC[0]-pointB[0])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) += 
						(tpc[2]-tpb[2])*(pointC[0]-pointB[0])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id] += 
						((tpb[2]-tpc[2])*(pointC[2]-pointB[2])*pointB[0] + (tpc[1]-tpb[1])*(pointC[0]-pointB[0])*pointB[1] - 
						(tpb[2]-tpc[2])*(pointC[0]-pointB[0])*pointB[2] - (tpc[1]-tpb[1])*(pointC[1]-pointB[1])*pointB[0] ) * ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);

					// for y channel
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id) += 
						(tpc[0]-tpb[0])*(pointC[1]-pointB[1])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) += 
						((tpb[0]-tpc[0])*(pointC[0]-pointB[0])-(tpc[2]-tpb[2])*(pointC[2]-pointB[2]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) += 
						(tpc[2]-tpb[2])*(pointC[1]-pointB[1])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id+ObjTriMesh.n_vertices()] += 
						((tpc[2]-tpb[2])*(pointC[1]-pointB[1])*pointB[2] + (tpb[0]-tpc[0])*(pointC[0]-pointB[0])*pointB[1] - 
						(tpc[2]-tpb[2])*(pointC[2]-pointB[2])*pointB[1] - (tpb[0]-tpc[0])*(pointC[1]-pointB[1])*pointB[0] ) * ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);

					// for z channel
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) += 
						(tpc[0]-tpb[0])*(pointC[2]-pointB[2])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) += 
						((tpc[1]-tpb[1])*(pointC[2]-pointB[2]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) += 
						((tpb[1]-tpc[1])*(pointC[1]-pointB[1])-(tpc[0]-tpb[0])*(pointC[0]-pointB[0]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id+ObjTriMesh.n_vertices()*2] += 
						((tpb[1]-tpc[1])*(pointC[1]-pointB[1])*pointB[2] + (tpc[0]-tpb[0])*(pointC[2]-pointB[2])*pointB[0] - 
						(tpb[1]-tpc[1])*(pointC[2]-pointB[2])*pointB[1] - (tpc[0]-tpb[0])*(pointC[0]-pointB[0])*pointB[2] ) * ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
				}

				// put the known value to b, ba*\beta*n^{in}(vi) * 1/k \sum jacobi
				b[vertex_id] += 
					ibeta*(Vertex_PSNormal[vertex_id][1]*(tpb[2]-tpc[2]) + Vertex_PSNormal[vertex_id][2]*(tpc[1]-tpb[1]))/TwoSumNeiFaceArea; // x channel
				b[vertex_id+ObjTriMesh.n_vertices()] += 
					ibeta*(Vertex_PSNormal[vertex_id][0]*(tpc[2]-tpb[2]) + Vertex_PSNormal[vertex_id][2]*(tpb[0]-tpc[0]))/TwoSumNeiFaceArea; // y channel
				b[vertex_id+ObjTriMesh.n_vertices()*2] += 
					ibeta*(Vertex_PSNormal[vertex_id][0]*(tpb[1]-tpc[1]) + Vertex_PSNormal[vertex_id][1]*(tpc[0]-tpb[0]))/TwoSumNeiFaceArea; // z channel
			}
			// end of norm term
		}
		// finish the reconstruction of the left matrix
		solver.getMatA() = AlphaNormVLap;
		solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
		solver.init();
		// solve the nonlinear equation to calculate the new vertex position in vx;
		solver.solve(b, vx);
		solver.clear();
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			curV[v_it.handle().idx()][0] = vx[v_it.handle().idx()];
			curV[v_it.handle().idx()][1] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()];
			curV[v_it.handle().idx()][2] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()*2];
		}
		MyMesh T_Mesh = ObjTriMesh;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			T_Mesh.set_point(v_it.handle(), curV[v_it.handle().idx()]);
		}
		OpenMesh::IO::Options write_options;
		write_options.set(OpenMesh::IO::Options::VertexNormal); 
		char buffer[255];
		sprintf(buffer, "TempMeshResult_%d_%.0f_%.0f.obj", outL, ialpha, ibeta);
		if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
			std::cerr << "Cannot write mesh to file " << buffer << std::endl;
		}
		StopCond = 0.0;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			StopCond += (oldV[v_it.handle().idx()]-curV[v_it.handle().idx()]).norm()*m_BCDArea[v_it.handle().idx()];
		}
		cout << "The error is: " << StopCond << "; the threshold is: " << outTole << endl;
	} while (StopCond > outTole);
}

void MeshTVRefine::MeshRefineByFaceNormal(double ialpha, double ibeta)
{
	this->UpdateBCDArea();
	this->UpdateEdgeLaplaceWeight();
	// refine the mesh by normal direction, using least square method
	double *b = new double[ObjTriMesh.n_vertices()*3]; 
	double *vx = new double[ObjTriMesh.n_vertices()*3]; 
	numc::RowMat<double> AlphaNormVLap; AlphaNormVLap.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::SparseSolver solver;
	// Initialization process
	vector < OpenMesh::Vec3f > oldV, curV; oldV.clear(); curV.clear();
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		oldV.push_back(ObjTriMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
		curV.push_back(ObjTriMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
	}
	double StopCond = 0.00000;  int innerL = 1; int outL = 0;
	double outTole = 5e-10;
	OpenMesh::Vec3f pointA, pointB, pointC, tpa, tpb, tpc; int ida, idb, idc;
	cout << "Start the mesh refinement by face normal process: ";
	do {
		outL ++;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			oldV[v_it.handle().idx()] = curV[v_it.handle().idx()];
		}
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			b[v_it.handle().idx()] =							ialpha*ObjTriMesh.point(v_it.handle())[0];	// for the x channel
			b[v_it.handle().idx()+ObjTriMesh.n_vertices()] =	ialpha*ObjTriMesh.point(v_it.handle())[1];	// for the y channel
			b[v_it.handle().idx()+ObjTriMesh.n_vertices()*2] =	ialpha*ObjTriMesh.point(v_it.handle())[2];	// for the z channel
		}
		// Construct the nonlinear matrix A
		//cout << AlphaNormVLap.clearZero() << " non-zero elements are erased." << endl;
		AlphaNormVLap *= 0.0;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
			// fidelity term 
			AlphaNormVLap(vertex_id, vertex_id) =														ialpha; // for x channel
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) =		ialpha; // for y channel
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) =	ialpha; // for z channel
			//end of fidelity term
		}

		for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
			// normal term
			MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(f_it.handle());
			pointA = curV[cfvIt.handle().idx()];			ida = cfvIt.handle().idx();		tpa = pointA;
			pointB = curV[(++cfvIt).handle().idx()];		idb = cfvIt.handle().idx();		tpb = pointB;
			pointC = curV[(++cfvIt).handle().idx()];		idc = cfvIt.handle().idx();		tpc = pointC;
			double FaceArea = OpenMesh::cross((pointC-pointB), (pointA-pointB)).norm();	double TwoSumNeiFaceArea = FaceArea * 2;
			int vertex_id = ida;
			AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()) = 0.0;
			AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) = 0.0; // set the initial value for v_i.y and v_i.z; v_i.x is already set
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id) = 0.0;
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) = 0.0;
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set

			// for x channel
			AlphaNormVLap(vertex_id, vertex_id) += 
				((tpb[2]-tpc[2])*(pointC[2]-pointB[2]) - (tpc[1]-tpb[1])*(pointC[1]-pointB[1]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()) += 
				(tpc[1]-tpb[1])*(pointC[0]-pointB[0])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			AlphaNormVLap(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) += 
				(tpc[2]-tpb[2])*(pointC[0]-pointB[0])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			b[vertex_id] += 
				((tpb[2]-tpc[2])*(pointC[2]-pointB[2])*pointB[0] + (tpc[1]-tpb[1])*(pointC[0]-pointB[0])*pointB[1] - 
				(tpb[2]-tpc[2])*(pointC[0]-pointB[0])*pointB[2] - (tpc[1]-tpb[1])*(pointC[1]-pointB[1])*pointB[0] ) * ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			b[vertex_id] += 
				ibeta*(Vertex_PSNormal[vertex_id][1]*(tpb[2]-tpc[2]) + Vertex_PSNormal[vertex_id][2]*(tpc[1]-tpb[1]))/TwoSumNeiFaceArea; // x channel

			// for y channel
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id) += 
				(tpc[0]-tpb[0])*(pointC[1]-pointB[1])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) += 
				((tpb[0]-tpc[0])*(pointC[0]-pointB[0])-(tpc[2]-tpb[2])*(pointC[2]-pointB[2]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) += 
				(tpc[2]-tpb[2])*(pointC[1]-pointB[1])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			b[vertex_id+ObjTriMesh.n_vertices()] += 
				((tpc[2]-tpb[2])*(pointC[1]-pointB[1])*pointB[2] + (tpb[0]-tpc[0])*(pointC[0]-pointB[0])*pointB[1] - 
				(tpc[2]-tpb[2])*(pointC[2]-pointB[2])*pointB[1] - (tpb[0]-tpc[0])*(pointC[1]-pointB[1])*pointB[0] ) * ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			b[vertex_id+ObjTriMesh.n_vertices()] += 
				ibeta*(Vertex_PSNormal[vertex_id][0]*(tpc[2]-tpb[2]) + Vertex_PSNormal[vertex_id][2]*(tpb[0]-tpc[0]))/TwoSumNeiFaceArea; // y channel

			// for z channel
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) += 
				(tpc[0]-tpb[0])*(pointC[2]-pointB[2])*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) += 
				((tpc[1]-tpb[1])*(pointC[2]-pointB[2]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) += 
				((tpb[1]-tpc[1])*(pointC[1]-pointB[1])-(tpc[0]-tpb[0])*(pointC[0]-pointB[0]))*ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			b[vertex_id+ObjTriMesh.n_vertices()*2] += 
				((tpb[1]-tpc[1])*(pointC[1]-pointB[1])*pointB[2] + (tpc[0]-tpb[0])*(pointC[2]-pointB[2])*pointB[0] - 
				(tpb[1]-tpc[1])*(pointC[2]-pointB[2])*pointB[1] - (tpc[0]-tpb[0])*(pointC[0]-pointB[0])*pointB[2] ) * ibeta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
			b[vertex_id+ObjTriMesh.n_vertices()*2] += 
				ibeta*(Vertex_PSNormal[vertex_id][0]*(tpb[1]-tpc[1]) + Vertex_PSNormal[vertex_id][1]*(tpc[0]-tpb[0]))/TwoSumNeiFaceArea; // z channel
		}
		// finish the reconstruction of the left matrix
		solver.getMatA() = AlphaNormVLap;
		solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
		solver.init();
		// solve the nonlinear equation to calculate the new vertex position in vx;
		solver.solve(b, vx);
		solver.clear();
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			curV[v_it.handle().idx()][0] = vx[v_it.handle().idx()];
			curV[v_it.handle().idx()][1] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()];
			curV[v_it.handle().idx()][2] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()*2];
		}
		MyMesh T_Mesh = ObjTriMesh;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			T_Mesh.set_point(v_it.handle(), curV[v_it.handle().idx()]);
		}
		OpenMesh::IO::Options write_options;
		write_options.set(OpenMesh::IO::Options::VertexNormal); 
		char buffer[255];
		sprintf(buffer, "TempMeshResult_%d_%.0f_%.0f.obj", outL, ialpha, ibeta);
		if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
			std::cerr << "Cannot write mesh to file " << buffer << std::endl;
		}
		StopCond = 0.0;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			StopCond += (oldV[v_it.handle().idx()]-curV[v_it.handle().idx()]).norm()*m_BCDArea[v_it.handle().idx()];
		}
		cout << "The error is: " << StopCond << "; the threshold is: " << outTole << endl;
	} while (StopCond > outTole);
	return;
}

void MeshTVRefine::LocalNormalTVRefinement(const double ialpha, const double igamma, const double isigma)
{
	// refine the fitted surface normal locally using the neighbor points normal
	// min \sum_{i}(a/2 * (n_i - n_i')^2 + b/2 \sum_{j\in N(i)} (n_i - n_j)^2 )
	// need to load LoadColorFile and LoadPSNormalFile first, then refine the PS normal
	vector < vector<OpenMesh::Vec3f> > lambda; lambda.clear();			// x,y,z channels in each triangular face
	vector < vector<OpenMesh::Vec3f> > npfield; npfield.clear();			// x,y,z channels in each triangular face
	vector < OpenMesh::Vec3f > oldN, curN; oldN.clear(); curN.clear();
	OpenMesh::Vec3f wx, wy, wz;
	double StopCond = 0.00000;  int innerL = 1; int outL = 0;
	double outTole = 1e-3;
	// Initialization process
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		oldN.push_back(Vertex_PSNormal[v_it.handle().idx()]);
		curN.push_back(Vertex_PSNormal[v_it.handle().idx()]);
	}
	lambda.resize(ObjTriMesh.n_faces()); npfield.resize(ObjTriMesh.n_faces());
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		lambda[f_it.handle().idx()].clear(); npfield[f_it.handle().idx()].clear();
		for (int i = 0; i < 3; ++i) {	// for the x, y, z channels
			lambda[f_it.handle().idx()].push_back(OpenMesh::Vec3f(0.0,0.0,0.0));
			npfield[f_it.handle().idx()].push_back(OpenMesh::Vec3f(0.0,0.0,0.0));
		}
	}

	double *nb = new double[ObjTriMesh.n_vertices()*3];
	double *nx = new double[ObjTriMesh.n_vertices()*3];
	numc::RowMat<double> AlphaNormVLap; AlphaNormVLap.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::SparseSolver solver;
	cout << "Start the mesh normal local refinement process: " << endl;
	do { // Start the iteration
		outL ++;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			oldN[v_it.handle().idx()] = curN[v_it.handle().idx()];
		}
		for (int iter_step = 0; iter_step < innerL; ++ iter_step) { // The inner loop
			cout << "(" << outL << ", " << iter_step << "): ";
			cout << "Solve the nonlinear v-subproblem ";
			for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
				for (int i = 0; i < 3; ++ i) {	// for the x, y, z channels
					npfield[f_it.handle().idx()][i] = lambda[f_it.handle().idx()][i] + npfield[f_it.handle().idx()][i]*igamma;
				}
			}
			// Calculate the divergence of pfield, then set it to b as the right part
			CalDivengence(npfield);
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				nb[v_it.handle().idx()] = 
					m_BCDArea[v_it.handle().idx()]*(ialpha*Vertex_PSNormal[v_it.handle().idx()][0] - m_divergence[v_it.handle().idx()][0]);	// for the x channel
				nb[v_it.handle().idx()+ObjTriMesh.n_vertices()] = 
					m_BCDArea[v_it.handle().idx()]*(ialpha*Vertex_PSNormal[v_it.handle().idx()][1] - m_divergence[v_it.handle().idx()][1]);	// for the y channel
				nb[v_it.handle().idx()+ObjTriMesh.n_vertices()*2] = 
					m_BCDArea[v_it.handle().idx()]*(ialpha*Vertex_PSNormal[v_it.handle().idx()][2] - m_divergence[v_it.handle().idx()][2]);	// for the z channel
			}
#ifdef _DEBUG
			fstream of("nb.txt",std::ios::out);
			if (!of) {
				cout << "Can not open b.txt to save data..." << endl;
				return;
			}
			for (int i = 0; i < ObjTriMesh.n_vertices()*3; ++ i) {
				of << nb[i] << endl;
			}
			of.close();
#endif
			// Construct the nonlinear matrix A
			//cout << AlphaNormVLap.clearZero() << " non-zero elements are erased." << endl;
			AlphaNormVLap *= 0.0;
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
				// fidelity term and laplacian term: the laplace coefficient of vi and vj is std::exp(-std::pow((I(vi)-I(vj)), 2.0)/isigma)
				for (MyMesh::ConstVertexVertexIter cvv_it = ObjTriMesh.cvv_iter(v_it.handle()); cvv_it; ++ cvv_it) {
					double ijweight = std::exp(-std::pow((Vertex_Intensity[v_it.handle().idx()] - Vertex_Intensity[cvv_it.handle().idx()]), 2.0)/isigma);
					AlphaNormVLap(vertex_id, cvv_it.handle().idx()) = -igamma * ijweight; 
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), cvv_it.handle().idx()+ObjTriMesh.n_vertices()) = -igamma * ijweight; 
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, cvv_it.handle().idx()+ObjTriMesh.n_vertices()*2) = -igamma * ijweight; 
					diag_value -= ijweight;
				}
				AlphaNormVLap(vertex_id, vertex_id) = ialpha*m_BCDArea[vertex_id] - igamma*diag_value; // for x channel
				AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) = ialpha*m_BCDArea[vertex_id] - igamma*diag_value; // for y channel
				AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) = ialpha*m_BCDArea[vertex_id] - igamma*diag_value; // for z channel
				// end of fidelity term and laplacian term
			}
			// finish the reconstruction of the left matrix
			solver.getMatA() = AlphaNormVLap;
			solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
			solver.init();
			// solve the nonlinear equation to calculate the new vertex position in vx;
			solver.solve(nb, nx);
			solver.clear();
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				curN[v_it.handle().idx()][0] = nx[v_it.handle().idx()];
				curN[v_it.handle().idx()][1] = nx[v_it.handle().idx()+ObjTriMesh.n_vertices()];
				curN[v_it.handle().idx()][2] = nx[v_it.handle().idx()+ObjTriMesh.n_vertices()*2];
				curN[v_it.handle().idx()].normalize();
			}
			// test for alm tv coefficient
			double udiff = 0.0;
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				udiff += (curN[v_it.handle().idx()] - Vertex_PSNormal[v_it.handle().idx()]).norm();
			}
			cout << "The udiff is: " << udiff << endl; // end

			cout << "Solve the p sub problem ";
			double pdiff = 0.0;
			CalGradient(curN, m_TriGradient);
			for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
				wx = m_TriGradient[f_it.handle().idx()][0] - lambda[f_it.handle().idx()][0]/igamma;
				wy = m_TriGradient[f_it.handle().idx()][1] - lambda[f_it.handle().idx()][1]/igamma;
				wz = m_TriGradient[f_it.handle().idx()][2] - lambda[f_it.handle().idx()][2]/igamma;
				vector<OpenMesh::Vec3f> old_ppfield = npfield[f_it.handle().idx()];
				if (std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm()) <= 1/igamma) {
					npfield[f_it.handle().idx()][0] = npfield[f_it.handle().idx()][1] = npfield[f_it.handle().idx()][2] = OpenMesh::Vec3f(0.0,0.0,0.0);
				} else {
					float temp_norm = std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm());
					npfield[f_it.handle().idx()][0] = wx*(1 - 1/(igamma*temp_norm));
					npfield[f_it.handle().idx()][1] = wy*(1 - 1/(igamma*temp_norm));
					npfield[f_it.handle().idx()][2] = wz*(1 - 1/(igamma*temp_norm));
				}
				pdiff += (npfield[f_it.handle().idx()][0] - old_ppfield[0]).norm() + 
						 (npfield[f_it.handle().idx()][1] - old_ppfield[1]).norm() +
						 (npfield[f_it.handle().idx()][2] - old_ppfield[2]).norm();
			} cout << "The pdiff is: " << pdiff << endl;
			cout << "Done..." << endl;
		}
		//cout << endl << "Update Lambda." << endl;
		// update lambda
		for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
			lambda[f_it.handle().idx()][0] += (npfield[f_it.handle().idx()][0]-m_TriGradient[f_it.handle().idx()][0])*igamma;
			lambda[f_it.handle().idx()][1] += (npfield[f_it.handle().idx()][1]-m_TriGradient[f_it.handle().idx()][1])*igamma;
			lambda[f_it.handle().idx()][2] += (npfield[f_it.handle().idx()][2]-m_TriGradient[f_it.handle().idx()][2])*igamma;
		}
		StopCond = 0.00000;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			StopCond += (oldN[v_it.handle().idx()]-curN[v_it.handle().idx()]).sqrnorm()*m_BCDArea[v_it.handle().idx()];
		}
		cout << "The error is: " << StopCond << "; the threshold is: " << outTole << endl;
	} while (StopCond > outTole);
	solver.clear();

	fstream of("localtvrefinedpsnormal.txt",std::ios::out);
	if (!of) {
		cout << "Can not open localrefinedpsnormal.txt to save data..." << endl;
		return;
	}
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		of << curN[v_it.handle().idx()][0] << "  " << curN[v_it.handle().idx()][1] << "  " << curN[v_it.handle().idx()][2] << "  " 
			<< Vertex_PSNormal[v_it.handle().idx()][0] << "  " << Vertex_PSNormal[v_it.handle().idx()][1] << "  " << Vertex_PSNormal[v_it.handle().idx()][2] << endl;
		Vertex_PSNormal[v_it.handle().idx()] = curN[v_it.handle().idx()]; 
	}
	cout << "Done..." << endl;
	return;
}

double MeshTVRefine::CalculateVEnergy(MyMesh &CurMesh, const std::vector<OpenMesh::Vec3f> InputV)
{
	if (CurMesh.n_vertices() != InputV.size()) {
		cout << "The size of CurV and InputV is not equal..." << endl;
		return 0.0;
	}
	double result = 0.00;
	for (MyMesh::VertexIter v_it = CurMesh.vertices_begin(); v_it != CurMesh.vertices_end(); ++ v_it) {
		result += (CurMesh.point(v_it.handle()) - InputV[v_it.handle().idx()]).norm();
	}
	return result;
}

double MeshTVRefine::CalculateNEnergy(MyMesh &CurMesh, const std::vector<OpenMesh::Vec3f> InputN)
{
	if (CurMesh.n_vertices() != InputN.size()) {
		cout << "The size of CurN and InputN is not equal..." << endl;
		return 0.0;
	}
	CurMesh.request_face_normals();
	CurMesh.update_face_normals();
	CurMesh.request_vertex_normals();
	// let the mesh update the normals
	CurMesh.update_vertex_normals();
	double result = 0.00;
	for (MyMesh::VertexIter v_it = CurMesh.vertices_begin(); v_it != CurMesh.vertices_end(); ++ v_it) {
		result += (CurMesh.normal(v_it.handle()) - InputN[v_it.handle().idx()]).norm();
	}
	return result;
}

double MeshTVRefine::CalculateTVEnergy(MyMesh &CurMesh)
{
	vector< map<int, OpenMesh::Vec3f> > temp_TPPIBG;
	this->BuildTPPIBG(CurMesh, temp_TPPIBG);
	vector < OpenMesh::Vec3f > CurV; CurV.clear();
	for (MyMesh::VertexIter v_it = CurMesh.vertices_begin(); v_it != CurMesh.vertices_end(); ++ v_it) {
		CurV.push_back(CurMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
	}
	this->CalGradient(CurV, m_TriGradient);
	double result = 0.0;
	for (int idx = 0; idx < m_TriGradient.size(); ++ idx) {
		for (int ii = 0; ii < m_TriGradient[idx].size(); ++ ii) {
			result += m_TriGradient[idx][ii].norm();
		}
	}
	return result;
}

void MeshTVRefine::ALMTVMeshRefine(double m_alpha, double m_beta, double m_gamma, int iter_step = 0, double m_eta = 1.0)
{
	this->UpdateEdgeLaplaceWeight();
	this->UpdateTriangleArea();
	this->UpdateBCDArea();
	this->BuildTPPIBG(ObjTriMesh, m_TPPIBG);
	vector < vector<OpenMesh::Vec3f> > lambda; lambda.clear();			// x,y,z channels in each triangular face
	vector < vector<OpenMesh::Vec3f> > pfield; pfield.clear();			// x,y,z channels in each triangular face
	vector < OpenMesh::Vec3f > inputV, oldV, curV; inputV.clear(); oldV.clear(); curV.clear();
	OpenMesh::Vec3f pointA , pointB , pointC; int ida, idb, idc;
	OpenMesh::Vec3f vi , vi1 , vi2; int idi, idi1, idi2; 
	OpenMesh::Vec3f wx, wy, wz;		char buffer[255];
	double StopCond = 0.00000;  int innerL = 1; int outL = 0;
	double outTole = 5e-11;
	// Initialization process
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		oldV.push_back(ObjTriMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
		curV.push_back(ObjTriMesh.point(v_it.handle()));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
		inputV.push_back(ObjTriMesh.point(v_it.handle()));
	}
	lambda.resize(ObjTriMesh.n_faces()); pfield.resize(ObjTriMesh.n_faces());
	//this->CalGradient(curV, m_TriGradient);
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		lambda[f_it.handle().idx()].clear(); pfield[f_it.handle().idx()].clear();
		for (int i = 0; i < 3; ++i) {	// for the x, y, z channels
			lambda[f_it.handle().idx()].push_back(OpenMesh::Vec3f(0.0,0.0,0.0));
			pfield[f_it.handle().idx()].push_back(OpenMesh::Vec3f(0.0,0.0,0.0));//(m_TriGradient[f_it.handle().idx()][i]);
		}
	}

	double *b = new double[ObjTriMesh.n_vertices()*3];
	double *vx = new double[ObjTriMesh.n_vertices()*3];
	numc::RowMat<double> AlphaNormVLap; AlphaNormVLap.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::RowMat<double> NormVMatrix; NormVMatrix.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::SparseSolver solver;
	cout << "Start the mesh refinement process: " << endl;
	do { // Start the iteration
		outL ++;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			oldV[v_it.handle().idx()] = curV[v_it.handle().idx()];
		}
		for (int iter = 0; iter < innerL; ++ iter) { // The inner loop
			cout << "(" << outL << ", " << iter << "): " << endl;
			cout << "Solve the u sub problem " ;
			for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
				for (int i = 0; i < 3; ++ i) {	// for the x, y, z channels
					pfield[f_it.handle().idx()][i] = lambda[f_it.handle().idx()][i] + pfield[f_it.handle().idx()][i]*m_gamma;
				}
			}
			// Calculate the divergence of pfield, then set it to b as the right part
			CalDivengence(pfield);
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				b[v_it.handle().idx()] = 
					m_BCDArea[v_it.handle().idx()]*(m_alpha*inputV[v_it.handle().idx()][0] - m_divergence[v_it.handle().idx()][0]);	// for the x channel
				b[v_it.handle().idx()+ObjTriMesh.n_vertices()] = 
					m_BCDArea[v_it.handle().idx()]*(m_alpha*inputV[v_it.handle().idx()][1] - m_divergence[v_it.handle().idx()][1]);	// for the y channel
				b[v_it.handle().idx()+ObjTriMesh.n_vertices()*2] = 
					m_BCDArea[v_it.handle().idx()]*(m_alpha*inputV[v_it.handle().idx()][2] - m_divergence[v_it.handle().idx()][2]);	// for the z channel
			}
#ifdef _DEBUG
			fstream of("b.txt",std::ios::out);
			if (!of) {
				cout << "Can not open b.txt to save data..." << endl;
				return;
			}
			for (int i = 0; i < ObjTriMesh.n_vertices()*3; ++ i) {
				of << b[i] << endl;
			}
			of.close();
#endif
			for (int ii = 0; ii < 1; ii++) {
				// Construct the nonlinear matrix A
				//cout << AlphaNormVLap.clearZero() << " non-zero elements are erased." << endl;
				AlphaNormVLap *= 0.0;
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
					// fidelity term and laplacian term: the laplace coefficient of vi and vj is 1/2 * (cot(\alpha_{ij})+cot(\beta_{ij}))
					for (MyMesh::VertexOHalfedgeIter vohe_it = ObjTriMesh.voh_iter(v_it.handle()); vohe_it; ++vohe_it) {
						int vv_id = ObjTriMesh.to_vertex_handle(vohe_it.handle()).idx();
						int edge_id = ObjTriMesh.edge_handle(vohe_it.handle()).idx(); //=ObjTriMesh.edge_handle(ObjTriMesh.find_halfedge(v_it.handle(), vv_it.handle())).idx();
						AlphaNormVLap(vertex_id, vv_id) = -m_gamma * m_edgelaplaceweight[edge_id]; 
						AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vv_id+ObjTriMesh.n_vertices()) = -m_gamma * m_edgelaplaceweight[edge_id]; 
						AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vv_id+ObjTriMesh.n_vertices()*2) = -m_gamma * m_edgelaplaceweight[edge_id]; 
						diag_value -= m_edgelaplaceweight[edge_id];
					}
					AlphaNormVLap(vertex_id, vertex_id) = m_alpha*m_BCDArea[vertex_id] - m_gamma*diag_value; // for x channel
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) = m_alpha*m_BCDArea[vertex_id] - m_gamma*diag_value; // for y channel
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) = m_alpha*m_BCDArea[vertex_id] - m_gamma*diag_value; // for z channel
					//end of fidelity term and laplacian term
				}
				//AlphaNormVLap.SaveToFile("Matrix_00.txt");

				// normal term, need to multiply bcdarea for each vertex
				NormVMatrix *= 0.0;
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
					// first for all neighbor triangle, store its vertex in counterclock wise, A(current vi), B, C
					vector< vector<OpenMesh::Vec3f> >	NeiFaceVertexList;		NeiFaceVertexList.clear();
					vector< vector<int> >				NeiFaceVertexIDList;	NeiFaceVertexIDList.clear();
					vector< double >					NeiFaceAreaList;		NeiFaceAreaList.clear();
					double								SumNeiFaceArea = 0.0000;
					OpenMesh::Vec3f Temp_sum_norm = OpenMesh::Vec3f(0.0,0.0,0.0);
					for (MyMesh::VertexFaceIter vf_it = ObjTriMesh.vf_iter(v_it.handle()); vf_it; ++ vf_it) {
						vector<OpenMesh::Vec3f> FaceVertexList;		FaceVertexList.clear();
						vector<int>				FaceVertexIDList;	FaceVertexIDList.clear();
						MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(vf_it.handle());
						pointA = curV[cfvIt.handle().idx()];			ida = cfvIt.handle().idx();
						pointB = curV[(++cfvIt).handle().idx()];		idb = cfvIt.handle().idx();
						pointC = curV[(++cfvIt).handle().idx()];		idc = cfvIt.handle().idx();
						if (vertex_id == ida) {
							//pointA v_i; pointB v_i+1; pointC v_i+2
							FaceVertexList.push_back(pointA); FaceVertexList.push_back(pointB); FaceVertexList.push_back(pointC);
							FaceVertexIDList.push_back(ida); FaceVertexIDList.push_back(idb); FaceVertexIDList.push_back(idc);
						} else if (vertex_id == idb){ 
							FaceVertexList.push_back(pointB); FaceVertexList.push_back(pointC); FaceVertexList.push_back(pointA);
							FaceVertexIDList.push_back(idb); FaceVertexIDList.push_back(idc); FaceVertexIDList.push_back(ida);
						} else if (vertex_id == idc) {
							FaceVertexList.push_back(pointC); FaceVertexList.push_back(pointA); FaceVertexList.push_back(pointB);
							FaceVertexIDList.push_back(idc); FaceVertexIDList.push_back(ida); FaceVertexIDList.push_back(idb);
						} else {
							cout << "Something wrong happened in the OpenMesh iteration ..." << endl;
						}
						NeiFaceVertexList.push_back(FaceVertexList);
						NeiFaceVertexIDList.push_back(FaceVertexIDList);
						NeiFaceAreaList.push_back(OpenMesh::cross((pointB-pointA), (pointC-pointA)).norm());
						Temp_sum_norm += OpenMesh::cross((pointB-pointA), (pointC-pointA));
						//SumNeiFaceArea += OpenMesh::cross((pointB-pointA), (pointC-pointA)).norm();
					}
					SumNeiFaceArea = max((double)(Temp_sum_norm.norm()), 1e-10);
					// then the norm of each neighbor face is a 3*1 vector Nkj = (C-B)*(A-B), its derivative to a 3*3 jacob matrix Mki; then Mki*Nkj
					// actually Mki is a known matrix
					NormVMatrix(vertex_id, vertex_id) = 0.0;
					NormVMatrix(vertex_id, vertex_id+ObjTriMesh.n_vertices()) = 0.0;
					NormVMatrix(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) = 0.0; // set the initial value for v_i.y and v_i.z; v_i.x is already set
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id) = 0.0;
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) = 0.0;
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) = 0.0;
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) = 0.0;

					float KNum = (float)NeiFaceVertexList.size();  double BArea = m_BCDArea[vertex_id];
					double TwoSumNeiFaceArea = SumNeiFaceArea; // * 2
					OpenMesh::Vec3f tpa, tpb, tpc;
					for (int i = 0; i < NeiFaceVertexList.size(); ++ i) {
						tpa = NeiFaceVertexList[i][0];	tpb = NeiFaceVertexList[i][1];	tpc = NeiFaceVertexList[i][2]; // in Mki
						for (int j = 0; j < NeiFaceVertexList.size(); ++ j) {
							//cout << endl << i << ", " << j << ": ";
							// calculate the 36 term for the normal constraint
							pointA = NeiFaceVertexList[j][0];	pointB = NeiFaceVertexList[j][1];	pointC = NeiFaceVertexList[j][2]; // in Nkj
							// for x channel
							NormVMatrix(vertex_id, vertex_id) += 
								((tpb[2]-tpc[2])*(pointC[2]-pointB[2]) - (tpc[1]-tpb[1])*(pointC[1]-pointB[1]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpb[2]-tpc[2])*(pointC[2]-pointB[2]) - (tpc[1]-tpb[1])*(pointC[1]-pointB[1]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							NormVMatrix(vertex_id, vertex_id+ObjTriMesh.n_vertices()) += 
								(tpc[1]-tpb[1])*(pointC[0]-pointB[0])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << (tpc[1]-tpb[1])*(pointC[0]-pointB[0])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							NormVMatrix(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) += 
								(tpc[2]-tpb[2])*(pointC[0]-pointB[0])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << (tpc[2]-tpb[2])*(pointC[0]-pointB[0])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							b[vertex_id] += 
								((tpb[2]-tpc[2])*(pointC[2]-pointB[2])*pointB[0] + (tpc[1]-tpb[1])*(pointC[0]-pointB[0])*pointB[1] - 
								(tpb[2]-tpc[2])*(pointC[0]-pointB[0])*pointB[2] - (tpc[1]-tpb[1])*(pointC[1]-pointB[1])*pointB[0] ) * BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpb[2]-tpc[2])*(pointC[2]-pointB[2])*pointB[0] + (tpc[1]-tpb[1])*(pointC[0]-pointB[0])*pointB[1] - 
							//	(tpb[2]-tpc[2])*(pointC[0]-pointB[0])*pointB[2] - (tpc[1]-tpb[1])*(pointC[1]-pointB[1])*pointB[0] ) * BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) <<";  ";

							// for y channel
							NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id) += 
								(tpc[0]-tpb[0])*(pointC[1]-pointB[1])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << (tpc[0]-tpb[0])*(pointC[1]-pointB[1])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) += 
								((tpb[0]-tpc[0])*(pointC[0]-pointB[0])-(tpc[2]-tpb[2])*(pointC[2]-pointB[2]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpb[0]-tpc[0])*(pointC[0]-pointB[0])-(tpc[2]-tpb[2])*(pointC[2]-pointB[2]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) += 
								(tpc[2]-tpb[2])*(pointC[1]-pointB[1])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << (tpc[2]-tpb[2])*(pointC[1]-pointB[1])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							b[vertex_id+ObjTriMesh.n_vertices()] += 
								((tpc[2]-tpb[2])*(pointC[1]-pointB[1])*pointB[2] + (tpb[0]-tpc[0])*(pointC[0]-pointB[0])*pointB[1] - 
								(tpc[2]-tpb[2])*(pointC[2]-pointB[2])*pointB[1] - (tpb[0]-tpc[0])*(pointC[1]-pointB[1])*pointB[0] ) * BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpc[2]-tpb[2])*(pointC[1]-pointB[1])*pointB[2] + (tpb[0]-tpc[0])*(pointC[0]-pointB[0])*pointB[1] - 
							//	(tpc[2]-tpb[2])*(pointC[2]-pointB[2])*pointB[1] - (tpb[0]-tpc[0])*(pointC[1]-pointB[1])*pointB[0] ) * BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea)<<";  ";

							// for z channel
							NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) += 
								(tpc[0]-tpb[0])*(pointC[2]-pointB[2])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << (tpc[0]-tpb[0])*(pointC[2]-pointB[2])*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) += 
								((tpc[1]-tpb[1])*(pointC[2]-pointB[2]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpc[1]-tpb[1])*(pointC[2]-pointB[2]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) += 
								((tpb[1]-tpc[1])*(pointC[1]-pointB[1])-(tpc[0]-tpb[0])*(pointC[0]-pointB[0]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpb[1]-tpc[1])*(pointC[1]-pointB[1])-(tpc[0]-tpb[0])*(pointC[0]-pointB[0]))*BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << "  ";
							b[vertex_id+ObjTriMesh.n_vertices()*2] += 
								((tpb[1]-tpc[1])*(pointC[1]-pointB[1])*pointB[2] + (tpc[0]-tpb[0])*(pointC[2]-pointB[2])*pointB[0] - 
								(tpb[1]-tpc[1])*(pointC[2]-pointB[2])*pointB[1] - (tpc[0]-tpb[0])*(pointC[0]-pointB[0])*pointB[2] ) * BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
							//cout << ((tpb[1]-tpc[1])*(pointC[1]-pointB[1])*pointB[2] + (tpc[0]-tpb[0])*(pointC[2]-pointB[2])*pointB[0] - 
							//	(tpb[1]-tpc[1])*(pointC[2]-pointB[2])*pointB[1] - (tpc[0]-tpb[0])*(pointC[0]-pointB[0])*pointB[2] ) * BArea*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea) << endl;
						}
						//NormVMatrix.SaveToFile("NormalMatrix_01.txt");

						// put the known value to b, ba*\beta*n^{in}(vi) * 1/k \sum jacobi
						b[vertex_id] += 
							BArea*m_beta*(Vertex_PSNormal[vertex_id][1]*(tpb[2]-tpc[2]) + Vertex_PSNormal[vertex_id][2]*(tpc[1]-tpb[1]))/TwoSumNeiFaceArea; // x channel
						//cout << BArea*m_beta*(Vertex_PSNormal[vertex_id][1]*(tpb[2]-tpc[2]) + Vertex_PSNormal[vertex_id][2]*(tpc[1]-tpb[1]))/TwoSumNeiFaceArea << "  ";
						b[vertex_id+ObjTriMesh.n_vertices()] += 
							BArea*m_beta*(Vertex_PSNormal[vertex_id][0]*(tpc[2]-tpb[2]) + Vertex_PSNormal[vertex_id][2]*(tpb[0]-tpc[0]))/TwoSumNeiFaceArea; // y channel
						//cout << BArea*m_beta*(Vertex_PSNormal[vertex_id][0]*(tpc[2]-tpb[2]) + Vertex_PSNormal[vertex_id][2]*(tpb[0]-tpc[0]))/TwoSumNeiFaceArea<< "  ";
						b[vertex_id+ObjTriMesh.n_vertices()*2] += 
							BArea*m_beta*(Vertex_PSNormal[vertex_id][0]*(tpb[1]-tpc[1]) + Vertex_PSNormal[vertex_id][1]*(tpc[0]-tpb[0]))/TwoSumNeiFaceArea; // z channel
						//cout << BArea*m_beta*(Vertex_PSNormal[vertex_id][0]*(tpb[1]-tpc[1]) + Vertex_PSNormal[vertex_id][1]*(tpc[0]-tpb[0]))/TwoSumNeiFaceArea << endl;
					}

					//NormVMatrix.SaveToFile("NormalMatrix_02.txt");
					// end of norm term
				}
				//NormVMatrix.SaveToFile("NormalMatrix.txt");
				// finish the reconstruction of the left matrix
				AlphaNormVLap += NormVMatrix;
				//sprintf(buffer, "%s_%2d_%d.txt","Matrix", outL,ii);
				//AlphaNormVLap.SaveToFile(buffer);
				solver.getMatA() = AlphaNormVLap;
				solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
				solver.init();
				// solve the nonlinear equation to calculate the new vertex position in vx;
				solver.solve(b, vx);
				solver.clear();
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					curV[v_it.handle().idx()][0] = vx[v_it.handle().idx()];
					curV[v_it.handle().idx()][1] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()];
					curV[v_it.handle().idx()][2] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()*2];
				}
			}

			// test for alm tv coefficient
			double udiff = 0.0;
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				udiff += (curV[v_it.handle().idx()] - ObjTriMesh.point(v_it.handle())).sqrnorm();
			}
			cout << "Udiff is: " << udiff*m_alpha << endl; // end
			if (outL % 1 == 0) {
				MyMesh T_Mesh = ObjTriMesh;
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					T_Mesh.set_point(v_it.handle(), curV[v_it.handle().idx()]);
				}
				double v_energy = CalculateVEnergy(T_Mesh, inputV); 
				double n_energy = CalculateNEnergy(T_Mesh, this->Vertex_PSNormal);
				double tv_energy = CalculateTVEnergy(T_Mesh);
				cout << "The current energy is V:" << v_energy << ", N:" << n_energy << ", TV: " << tv_energy << "; Sum: " 
					<< m_alpha * v_energy + m_beta * n_energy + m_eta * tv_energy << endl;
				OpenMesh::IO::Options write_options;
				write_options.set(OpenMesh::IO::Options::VertexNormal); 
				char buffer[255];
				sprintf(buffer, "Results\\TempMeshResult_%d_%d_%d_%.0f_%.0f_%.3f.obj", iter_step, outL, iter, m_alpha, m_beta, m_gamma);
				if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
					std::cerr << "Cannot write mesh to file " << buffer << std::endl;
				}
			}

			cout << "Solve the p sub problem ...";
			CalGradient(curV,m_TriGradient); double pdiff = 0.0, avg_wxy = 0.0; int i1 = 0, i2 = 0;
			for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
				wx = m_TriGradient[f_it.handle().idx()][0] - lambda[f_it.handle().idx()][0]/m_gamma;
				wy = m_TriGradient[f_it.handle().idx()][1] - lambda[f_it.handle().idx()][1]/m_gamma;
				wz = m_TriGradient[f_it.handle().idx()][2] - lambda[f_it.handle().idx()][2]/m_gamma;
				avg_wxy += std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm());
				if (std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm()) <= m_eta/m_gamma) {
					pfield[f_it.handle().idx()][0] = pfield[f_it.handle().idx()][1] = pfield[f_it.handle().idx()][2] = OpenMesh::Vec3f(0.0,0.0,0.0);
					i1++;
				} else {
					float temp_norm = std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm());
					pfield[f_it.handle().idx()][0] = wx*(1 - m_eta/(m_gamma*temp_norm));
					pfield[f_it.handle().idx()][1] = wy*(1 - m_eta/(m_gamma*temp_norm));
					pfield[f_it.handle().idx()][2] = wz*(1 - m_eta/(m_gamma*temp_norm)); 
					i2++;
				}
				pdiff += (pfield[f_it.handle().idx()][0]).norm() + (pfield[f_it.handle().idx()][1]).norm() + (pfield[f_it.handle().idx()][2]).norm();
			} cout << "Pdiff is: " << pdiff << "  " << i1 << "/" << i2 << " under " << avg_wxy/ObjTriMesh.n_faces() << endl;
		} // end of innerL loop
		cout << endl << "Update Lambda." << endl;
		// update lambda
		for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
			lambda[f_it.handle().idx()][0] += (pfield[f_it.handle().idx()][0]-m_TriGradient[f_it.handle().idx()][0])*m_gamma;
			lambda[f_it.handle().idx()][1] += (pfield[f_it.handle().idx()][1]-m_TriGradient[f_it.handle().idx()][1])*m_gamma;
			lambda[f_it.handle().idx()][2] += (pfield[f_it.handle().idx()][2]-m_TriGradient[f_it.handle().idx()][2])*m_gamma;
		}
		StopCond = 0.00000;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			StopCond += (oldV[v_it.handle().idx()]-curV[v_it.handle().idx()]).sqrnorm()*m_BCDArea[v_it.handle().idx()];
		}
		cout << "The error is: " << StopCond << "; the threshold is: " << outTole << endl;
	} while (StopCond > outTole);
	solver.clear(); delete b, vx; AlphaNormVLap.clear();
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		ObjTriMesh.set_point(v_it.handle(), curV[v_it.handle().idx()]);
	}
	sprintf(buffer, "FinalMeshResult_%d_%.0f_%.0f_%.3f.obj", iter_step, m_alpha, m_beta, m_gamma);
	OpenMesh::IO::Options write_options;
	write_options.set(OpenMesh::IO::Options::VertexNormal); 
	if ( !OpenMesh::IO::write_mesh(ObjTriMesh, string(buffer), write_options) ) {
		std::cerr << "Cannot write mesh to file" << buffer << std::endl;
	}
	return;
}


void MeshTVRefine::ALMTVMeshRefineByFaceNormal(double m_alpha, double m_beta, double m_gamma, int iter_step = 0)
{
	this->UpdateEdgeLaplaceWeight();
	this->UpdateTriangleArea();
	this->UpdateBCDArea();
	this->BuildTPPIBG(ObjTriMesh, m_TPPIBG);
	vector < vector<OpenMesh::Vec3f> > lambda; lambda.clear();			// x,y,z channels in each triangular face
	vector < vector<OpenMesh::Vec3f> > pfield; pfield.clear();			// x,y,z channels in each triangular face
	vector < OpenMesh::Vec3f > inputV, oldV, curV; inputV.clear(); oldV.clear(); curV.clear();
	OpenMesh::Vec3f pointA , pointB , pointC, tpa, tpb, tpc; int ida, idb, idc;
	OpenMesh::Vec3f vi , vi1 , vi2; int idi, idi1, idi2; 
	OpenMesh::Vec3f wx, wy, wz;		char buffer[255];
	double StopCond = 0.00000;  int innerL = 1; int outL = 0;
	double outTole = 5e-10;
	// Initialization process
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		oldV.push_back(OpenMesh::Vec3f(0.0,0.0,0.0));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
		curV.push_back(OpenMesh::Vec3f(0.0,0.0,0.0));//(ObjTriMesh.point(v_it.handle()));//(OpenMesh::Vec3f(0.0,0.0,0.0));//
		inputV.push_back(ObjTriMesh.point(v_it.handle()));
	}
	lambda.resize(ObjTriMesh.n_faces()); pfield.resize(ObjTriMesh.n_faces());
	for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
		lambda[f_it.handle().idx()].clear(); pfield[f_it.handle().idx()].clear();
		for (int i = 0; i < 3; ++i) {	// for the x, y, z channels
			lambda[f_it.handle().idx()].push_back(OpenMesh::Vec3f(0.0,0.0,0.0));
			pfield[f_it.handle().idx()].push_back(OpenMesh::Vec3f(0.0,0.0,0.0));
		}
	}

	double *b = new double[ObjTriMesh.n_vertices()*3];
	double *vx = new double[ObjTriMesh.n_vertices()*3];
	numc::RowMat<double> AlphaNormVLap; AlphaNormVLap.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::RowMat<double> NormVMatrix; NormVMatrix.resize(ObjTriMesh.n_vertices()*3, ObjTriMesh.n_vertices()*3); 
	numc::SparseSolver solver;
	cout << "Start the mesh refinement process: " << endl;
	do { // Start the iteration
		outL ++;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			oldV[v_it.handle().idx()] = curV[v_it.handle().idx()];
		}
		for (int iter = 0; iter < innerL; ++ iter) { // The inner loop
			cout << "(" << outL << ", " << iter << "): " << endl;
			cout << "Solve the u sub problem " ;
			for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
				for (int i = 0; i < 3; ++ i) {	// for the x, y, z channels
					pfield[f_it.handle().idx()][i] = lambda[f_it.handle().idx()][i] + pfield[f_it.handle().idx()][i]*m_gamma;
				}
			}
			// Calculate the divergence of pfield, then set it to b as the right part
			CalDivengence(pfield);
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				b[v_it.handle().idx()] = 
					m_BCDArea[v_it.handle().idx()]*(m_alpha*inputV[v_it.handle().idx()][0] - m_divergence[v_it.handle().idx()][0]);	// for the x channel
				b[v_it.handle().idx()+ObjTriMesh.n_vertices()] = 
					m_BCDArea[v_it.handle().idx()]*(m_alpha*inputV[v_it.handle().idx()][1] - m_divergence[v_it.handle().idx()][1]);	// for the y channel
				b[v_it.handle().idx()+ObjTriMesh.n_vertices()*2] = 
					m_BCDArea[v_it.handle().idx()]*(m_alpha*inputV[v_it.handle().idx()][2] - m_divergence[v_it.handle().idx()][2]);	// for the z channel
			}
#ifdef _DEBUG
			fstream of("b.txt",std::ios::out);
			if (!of) {
				cout << "Can not open b.txt to save data..." << endl;
				return;
			}
			for (int i = 0; i < ObjTriMesh.n_vertices()*3; ++ i) {
				of << b[i] << endl;
			}
			of.close();
#endif
			for (int ii = 0; ii < 1; ii++) {
				// Construct the nonlinear matrix A
				//cout << AlphaNormVLap.clearZero() << " non-zero elements are erased." << endl;
				AlphaNormVLap *= 0.0;
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					const int vertex_id = v_it.handle().idx(); double diag_value = 0.0;
					// fidelity term and laplacian term: the laplace coefficient of vi and vj is 1/2 * (cot(\alpha_{ij})+cot(\beta_{ij}))
					for (MyMesh::VertexOHalfedgeIter vohe_it = ObjTriMesh.voh_iter(v_it.handle()); vohe_it; ++vohe_it) {
						int vv_id = ObjTriMesh.to_vertex_handle(vohe_it.handle()).idx();
						int edge_id = ObjTriMesh.edge_handle(vohe_it.handle()).idx(); //=ObjTriMesh.edge_handle(ObjTriMesh.find_halfedge(v_it.handle(), vv_it.handle())).idx();
						AlphaNormVLap(vertex_id, vv_id) = -m_gamma * m_edgelaplaceweight[edge_id]; 
						AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vv_id+ObjTriMesh.n_vertices()) = -m_gamma * m_edgelaplaceweight[edge_id]; 
						AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vv_id+ObjTriMesh.n_vertices()*2) = -m_gamma * m_edgelaplaceweight[edge_id]; 
						diag_value -= m_edgelaplaceweight[edge_id];
					}
					AlphaNormVLap(vertex_id, vertex_id) = m_alpha*m_BCDArea[vertex_id] - m_gamma*diag_value; // for x channel
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) = m_alpha*m_BCDArea[vertex_id] - m_gamma*diag_value; // for y channel
					AlphaNormVLap(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) = m_alpha*m_BCDArea[vertex_id] - m_gamma*diag_value; // for z channel
					//end of fidelity term and laplacian term
				}
				//AlphaNormVLap.SaveToFile("Matrix_00.txt");

				// normal term, need to multiply bcdarea for each vertex
				NormVMatrix *= 0.0;
				for (int ii = 0; ii < ObjTriMesh.n_vertices(); ++ ii) {
					NormVMatrix(ii, ii) = 0.0;
					NormVMatrix(ii, ii+ObjTriMesh.n_vertices()) = 0.0;
					NormVMatrix(ii, ii+ObjTriMesh.n_vertices()*2) = 0.0; // set the initial value for v_i.y and v_i.z; v_i.x is already set
					NormVMatrix(ii+ObjTriMesh.n_vertices(), ii) = 0.0;
					NormVMatrix(ii+ObjTriMesh.n_vertices(), ii+ObjTriMesh.n_vertices()) = 0.0;
					NormVMatrix(ii+ObjTriMesh.n_vertices(), ii+ObjTriMesh.n_vertices()*2) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
					NormVMatrix(ii+ObjTriMesh.n_vertices()*2, ii) = 0.0;
					NormVMatrix(ii+ObjTriMesh.n_vertices()*2, ii+ObjTriMesh.n_vertices()) = 0.0; //set initial value for v_i.y and v_i.z; v_i.x is already set
					NormVMatrix(ii+ObjTriMesh.n_vertices()*2, ii+ObjTriMesh.n_vertices()*2) = 0.0;
				}
				for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
					// normal term
					MyMesh::ConstFaceVertexIter cfvIt = ObjTriMesh.cfv_iter(f_it.handle());
					pointA = curV[cfvIt.handle().idx()];			ida = cfvIt.handle().idx();		tpa = pointA;
					pointB = curV[(++cfvIt).handle().idx()];		idb = cfvIt.handle().idx();		tpb = pointB;
					pointC = curV[(++cfvIt).handle().idx()];		idc = cfvIt.handle().idx();		tpc = pointC;
					double FaceArea = OpenMesh::cross((pointC-pointB), (pointA-pointB)).norm();	double TwoSumNeiFaceArea = FaceArea * 2;
					int vertex_id = ida;
					// for x channel
					NormVMatrix(vertex_id, vertex_id) += 
						((tpb[2]-tpc[2])*(pointC[2]-pointB[2]) - (tpc[1]-tpb[1])*(pointC[1]-pointB[1]))*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					NormVMatrix(vertex_id, vertex_id+ObjTriMesh.n_vertices()) += 
						(tpc[1]-tpb[1])*(pointC[0]-pointB[0])*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					NormVMatrix(vertex_id, vertex_id+ObjTriMesh.n_vertices()*2) += 
						(tpc[2]-tpb[2])*(pointC[0]-pointB[0])*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id] += 
						((tpb[2]-tpc[2])*(pointC[2]-pointB[2])*pointB[0] + (tpc[1]-tpb[1])*(pointC[0]-pointB[0])*pointB[1] - 
						(tpb[2]-tpc[2])*(pointC[0]-pointB[0])*pointB[2] - (tpc[1]-tpb[1])*(pointC[1]-pointB[1])*pointB[0] ) * m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id] += 
						m_beta*(Vertex_PSNormal[vertex_id][1]*(tpb[2]-tpc[2]) + Vertex_PSNormal[vertex_id][2]*(tpc[1]-tpb[1]))/TwoSumNeiFaceArea; // x channel

					// for y channel
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id) += 
						(tpc[0]-tpb[0])*(pointC[1]-pointB[1])*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()) += 
						((tpb[0]-tpc[0])*(pointC[0]-pointB[0])-(tpc[2]-tpb[2])*(pointC[2]-pointB[2]))*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices(), vertex_id+ObjTriMesh.n_vertices()*2) += 
						(tpc[2]-tpb[2])*(pointC[1]-pointB[1])*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id+ObjTriMesh.n_vertices()] += 
						((tpc[2]-tpb[2])*(pointC[1]-pointB[1])*pointB[2] + (tpb[0]-tpc[0])*(pointC[0]-pointB[0])*pointB[1] - 
						(tpc[2]-tpb[2])*(pointC[2]-pointB[2])*pointB[1] - (tpb[0]-tpc[0])*(pointC[1]-pointB[1])*pointB[0] ) * m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id+ObjTriMesh.n_vertices()] += 
						m_beta*(Vertex_PSNormal[vertex_id][0]*(tpc[2]-tpb[2]) + Vertex_PSNormal[vertex_id][2]*(tpb[0]-tpc[0]))/TwoSumNeiFaceArea; // y channel

					// for z channel
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id) += 
						(tpc[0]-tpb[0])*(pointC[2]-pointB[2])*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()) += 
						((tpc[1]-tpb[1])*(pointC[2]-pointB[2]))*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					NormVMatrix(vertex_id+ObjTriMesh.n_vertices()*2, vertex_id+ObjTriMesh.n_vertices()*2) += 
						((tpb[1]-tpc[1])*(pointC[1]-pointB[1])-(tpc[0]-tpb[0])*(pointC[0]-pointB[0]))*m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id+ObjTriMesh.n_vertices()*2] += 
						((tpb[1]-tpc[1])*(pointC[1]-pointB[1])*pointB[2] + (tpc[0]-tpb[0])*(pointC[2]-pointB[2])*pointB[0] - 
						(tpb[1]-tpc[1])*(pointC[2]-pointB[2])*pointB[1] - (tpc[0]-tpb[0])*(pointC[0]-pointB[0])*pointB[2] ) * m_beta/(TwoSumNeiFaceArea*TwoSumNeiFaceArea);
					b[vertex_id+ObjTriMesh.n_vertices()*2] += 
						m_beta*(Vertex_PSNormal[vertex_id][0]*(tpb[1]-tpc[1]) + Vertex_PSNormal[vertex_id][1]*(tpc[0]-tpb[0]))/TwoSumNeiFaceArea; // z channel
				}
				// finish the reconstruction of the left matrix
				AlphaNormVLap += NormVMatrix;
				sprintf(buffer, "%s_%2d_%d.txt","Matrix", outL,ii);
				AlphaNormVLap.SaveToFile(buffer);
				solver.getMatA() = AlphaNormVLap;
				solver.getMatA().mMtype = numc::CSRMatrix<double>::RealUnSymm;
				solver.init();
				// solve the nonlinear equation to calculate the new vertex position in vx;
				solver.solve(b, vx);
				solver.clear();
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					curV[v_it.handle().idx()][0] = vx[v_it.handle().idx()];
					curV[v_it.handle().idx()][1] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()];
					curV[v_it.handle().idx()][2] = vx[v_it.handle().idx()+ObjTriMesh.n_vertices()*2];
				}
			}

			// test for alm tv coefficient
			double udiff = 0.0;
			for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
				udiff += (curV[v_it.handle().idx()] - ObjTriMesh.point(v_it.handle())).sqrnorm();
			}
			cout << "Udiff is: " << udiff*m_alpha << endl; // end
			if (outL % 1 == 0) {
				MyMesh T_Mesh = ObjTriMesh;
				for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
					T_Mesh.set_point(v_it.handle(), curV[v_it.handle().idx()]);
				}
				OpenMesh::IO::Options write_options;
				write_options.set(OpenMesh::IO::Options::VertexNormal); 
				char buffer[255];
				sprintf(buffer, "TempMeshResult_%d_%d_%d_%.0f_%.0f_%.3f.obj", iter_step, outL, iter, m_alpha, m_beta, m_gamma);
				if ( !OpenMesh::IO::write_mesh(T_Mesh, string(buffer), write_options) ) {
					std::cerr << "Cannot write mesh to file " << buffer << std::endl;
				}
			}

			cout << "Solve the p sub problem ...";
			CalGradient(curV, m_TriGradient); double pdiff = 0.0, avg_wxy = 0.0; int i1 = 0, i2 = 0;
			for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
				wx = m_TriGradient[f_it.handle().idx()][0] - lambda[f_it.handle().idx()][0]/m_gamma;
				wy = m_TriGradient[f_it.handle().idx()][1] - lambda[f_it.handle().idx()][1]/m_gamma;
				wz = m_TriGradient[f_it.handle().idx()][2] - lambda[f_it.handle().idx()][2]/m_gamma;
				avg_wxy += std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm());
				if (std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm()) <= 1/m_gamma) {
					pfield[f_it.handle().idx()][0] = pfield[f_it.handle().idx()][1] = pfield[f_it.handle().idx()][2] = OpenMesh::Vec3f(0.0,0.0,0.0);
					i1++;
				} else {
					float temp_norm = std::sqrt(wx.sqrnorm()+wy.sqrnorm()+wz.sqrnorm());
					pfield[f_it.handle().idx()][0] = wx*(1 - 1/(m_gamma*temp_norm));
					pfield[f_it.handle().idx()][1] = wy*(1 - 1/(m_gamma*temp_norm));
					pfield[f_it.handle().idx()][2] = wz*(1 - 1/(m_gamma*temp_norm)); 
					i2++;
				}
				pdiff += (pfield[f_it.handle().idx()][0]).norm() + 
					(pfield[f_it.handle().idx()][1]).norm() +
					(pfield[f_it.handle().idx()][2]).norm();
			} cout << "Pdiff is: " << pdiff << "  " << i1 << "/" << i2 << " under " << avg_wxy/ObjTriMesh.n_faces() << endl;
		} // end of innerL loop
		cout << endl << "Update Lambda." << endl;
		// update lambda
		for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
			lambda[f_it.handle().idx()][0] += (pfield[f_it.handle().idx()][0]-m_TriGradient[f_it.handle().idx()][0])*m_gamma;
			lambda[f_it.handle().idx()][1] += (pfield[f_it.handle().idx()][1]-m_TriGradient[f_it.handle().idx()][1])*m_gamma;
			lambda[f_it.handle().idx()][2] += (pfield[f_it.handle().idx()][2]-m_TriGradient[f_it.handle().idx()][2])*m_gamma;
		}
		StopCond = 0.00000;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
			StopCond += (oldV[v_it.handle().idx()]-curV[v_it.handle().idx()]).sqrnorm()*m_BCDArea[v_it.handle().idx()];
		}
		cout << "The error is: " << StopCond << "; the threshold is: " << outTole << endl;
	} while (StopCond > outTole);
	solver.clear(); delete b, vx; AlphaNormVLap.clear();
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		ObjTriMesh.set_point(v_it.handle(), curV[v_it.handle().idx()]);
	}
	sprintf(buffer, "FinalMeshResult_%d_%.0f_%.0f_%.3f.obj", iter_step, m_alpha, m_beta, m_gamma);
	OpenMesh::IO::Options write_options;
	write_options.set(OpenMesh::IO::Options::VertexNormal); 
	if ( !OpenMesh::IO::write_mesh(ObjTriMesh, string(buffer), write_options) ) {
		std::cerr << "Cannot write mesh to file" << buffer << std::endl;
	}
	return;
}