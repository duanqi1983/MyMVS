#ifndef __PCLUTILITY_H__
#define __PCLUTILITY_H__

//struct MyPointType // My point type with xyz, normal, rgb and view point
//{
//	PCL_ADD_POINT4D;    // This adds the members x,y,z which can also be accessed using the point (which is float[4])
//	PCL_ADD_NORMAL4D;   // This adds the member normal[3] which can also be accessed using the point (which is float[4])
//	vector<int> views; // list of views
//	vector<int> px; // 2D coordinate X
//	vector<int> py; // 2D coordinate Y
//	double fx;
//	vector<Visible_Info> vec_visible_info;
//	double *vec_visible_color;
//	//EIGEN_MAKE_ALIGNED_OPERATOR_NEW   // make sure our new allocators are aligned
//} ; 
//inline std::ostream& operator << (std::ostream& os, const MyPointType& p)
//{
//	os << "(" << p.x << "," << p.y << "," << p.z << " - " << 
//		p.normal[0] << "," << p.normal[1] << "," << p.normal[2] << " - " ;//<< 
//		//p.rgb << p.intensity << p.curvature << p.vp_x << "," << p.vp_y << "," << p.vp_z<< ")";
//	return (os);
//}

extern pcl::PointCloud<pcl::PointNormal> my_pcl_cloud;
extern pcl::PolygonMesh my_pcl_mesh;
std::vector<int> parts;
std::vector<int> states;

void ConvertTracks2PCL(vector<track> &_tracks)
{
	my_pcl_cloud.points.resize(_tracks.size());
	for (size_t i = 0; i < my_pcl_cloud.points.size (); ++i)
	{
		my_pcl_cloud.points[i].x = _tracks[i].X[0];
		my_pcl_cloud.points[i].y = _tracks[i].X[1];
		my_pcl_cloud.points[i].z = _tracks[i].X[2];

		my_pcl_cloud.points[i].normal_x = _tracks[i].Norm[0];
		my_pcl_cloud.points[i].normal_y = _tracks[i].Norm[1];
		my_pcl_cloud.points[i].normal_z = _tracks[i].Norm[2];

		//my_pcl_cloud.points[i].r = _tracks[i].material_label * 25;
		//my_pcl_cloud.points[i].g = _tracks[i].material_label * 25;
		//my_pcl_cloud.points[i].b = _tracks[i].material_label * 25;
		//my_pcl_cloud.points[i].rgb = _tracks[i].material_label * 25;
	}
}

//template <typename PointClass>
void CloudTriangulation(pcl::PointCloud<pcl::PointNormal> &cloud_in)
{
	// Create search tree*
	pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
	tree2->setInputCloud (cloud_in.makeShared());

	//pcl::OrganizedFastMesh<pcl::PointNormal> ofm;
	//ofm.setInputCloud(cloud_in.makeShared());
	//ofm.performReconstruction(my_pcl_mesh);

	pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
	// Set the maximum distance between connected points (maximum edge length)
	gp3.setSearchRadius (5000);

	// Set typical values for the parameters
	gp3.setMu (4);
	gp3.setMaximumNearestNeighbors (1000);
	gp3.setMaximumSurfaceAngle(M_PI/2); // 60 degrees
	gp3.setMinimumAngle(M_PI/180); // 10 degrees
	gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
	gp3.setNormalConsistency(false);

	gp3.setInputCloud(cloud_in.makeShared());
	gp3.setSearchMethod (tree2);
	gp3.reconstruct (my_pcl_mesh);
	parts = gp3.getPartIDs();
	states= gp3.getPointStates();
}


template <typename PointClass>
void ViewCloud(const char* WindowTitle, pcl::PointCloud<PointClass> &cloud_show) {
#if !defined(_APPLE_)
	std::cerr << "Cloud Information: " << std::endl << cloud_show << std::endl;
	pcl::visualization::PointCloudGeometryHandlerXYZ<PointClass> handler (cloud_show.makeShared());
	//pcl::visualization::PointCloudColorHandlerRGBField<PointClass> handler(cloud_show.makeShared());
	pcl::visualization::PCLVisualizer viewer (WindowTitle);
	viewer.setBackgroundColor (0.5,0.5,0);
	viewer.addCoordinateSystem (1.0,0.2,0.3,0.5);
	viewer.addPointCloud (cloud_show.makeShared(),handler,"cloud_rgb");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud_rgb"); //points size 
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, "cloud_rgb"); //color (red)
	viewer.addText(WindowTitle,10,10,1.0,1.0,1.0,"cloud_rgb");
	viewer.initCameraParameters();
	viewer.spin();
	// Remove visualization data
	viewer.removePointCloud ("cloud_rgb");
#endif
}

void ViewPolygonMesh(const char* WindowTitle, pcl::PolygonMesh mesh_show) {
#if !defined(_APPLE_)
	//std::cerr << "Cloud Information: " << std::endl << cloud_show << std::endl;
	//pcl::visualization::PointCloudGeometryHandlerXYZ<PointClass> handler (cloud_show.makeShared());
	//pcl::visualization::PointCloudColorHandlerRGBField<PointClass> handler(cloud_show.makeShared());
	pcl::visualization::PCLVisualizer viewer (WindowTitle);
	viewer.setBackgroundColor (0.5,0.5,0);
	viewer.addCoordinateSystem (1.0,0.2,0.3,0.5);
	//viewer.addPointCloud (cloud_show.makeShared(),handler,"cloud_rgb");
	//viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud_rgb"); //points size 
	//viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, "cloud_rgb"); //color (red)
	viewer.addText(WindowTitle,10,10,1.0,1.0,1.0,"cloud_mesh");
	viewer.initCameraParameters();
	viewer.addPolygonMesh(mesh_show);
	viewer.spin();
	// Remove visualization data
	viewer.removeAllPointClouds();
#endif
}





























#endif //__PCLUTILITY_H__