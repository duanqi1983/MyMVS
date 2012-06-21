#ifndef _PARAMETERS_DEF
#define _PARAMETERS_DEF


#include "cv.h"
#include "highgui.h"
#include "cxcore.h"

#define USING_NAGC

#include <nag.h>
#include <nage02.h>
#include <nag_stdlib.h>
#include <nag_example_file_io.h>
#pragma comment(lib,"CLW6I09DA_nag.lib")
#pragma comment(lib,"nagc_nag_MD.lib")
#pragma comment(lib,"libdecimal.lib")
#pragma comment(lib,"libiomp5md.lib")
#pragma comment(lib,"libirc.lib")
#pragma comment(lib,"libmmd.lib")

//#include "MeshTVRefine.h"

//#include <pcl/point_types.h>
//#include <pcl/point_cloud.h>
//#include <pcl/PointIndices.h>
//
//#include <pcl/console/time.h>
//#include <pcl/console/parse.h>
//#include <pcl/console/print.h>
//
//#include <pcl/io/pcd_io.h>
//#include <pcl/io/ply_io.h>
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <pcl/visualization/point_cloud_handlers.h>
//
//#include <pcl/registration/icp.h>
//#include <pcl/registration/icp_nl.h>
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/registration/registration.h>
//#include <pcl/registration/correspondence_estimation.h>
//
//#include <pcl/kdtree/kdtree_flann.h>
//#include <pcl/kdtree/kdtree.h>
//
//#include <pcl/octree/octree.h>
//
//#include <pcl/surface/reconstruction.h>
//#include <pcl/surface/mls.h>
//#include <pcl/surface/gp3.h>
//#include <pcl/surface/organized_fast_mesh.h>
//
//#include <pcl/common/pca.h>
//#include <pcl/common/common_headers.h>
//#include <pcl/common/transforms.h>
//
//#include <pcl/filters/filter.h>
//#include <pcl/filters/passthrough.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/statistical_outlier_removal.h>
//
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/features/pfh.h>
//#include <pcl/features/feature.h>
//
//#include <pcl/keypoints/keypoint.h>
//#include <pcl/keypoints/narf_keypoint.h>
//#include <pcl/keypoints/sift_keypoint.h>

//#include <cholmod.h>
//#include <TriMesh.h>
//#include <TriMesh_algo.h>

//#include <CGAL/basic.h>
//#include <CGAL/Cartesian.h>
//#include <CGAL/Segment_2.h>
//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Polyhedron_3.h>

bool MiddleBuryData = true;
bool LiJianguoData = false;

//string inputFileName = "camera_estimated.txt"; string path="../../data/Monster";
//string inputFileName = "Monster_par.txt"; string path="../data/Monster_4";
//string inputFileName = "camera_estimated.txt"; string path="../../data/IronLion";
//string inputFileName = "camera_estimated.txt"; string path="../../data/SculptFace";
//string inputFileName = "SculptFace_par.txt"; string path="../../data/SculptFace_4";

string inputFileName;
string path;
//string inputFileName = "dinoSR_par.txt"; string path="../data/dinoSparseRing";
//string inputFileName = "dinoSR_par.txt"; string path="../../data/dinoSparseRingx4";
//string inputFileName = "templeSR_par.txt"; string path="../data/templeSparseRing";
//string inputFileName = "templeSR_par.txt"; string path="../../data/templeSparseRingx4";
//string inputFileName = "dinoR_par.txt"; string path="../../data/dinoRing";
//string inputFileName = "templeR_par.txt"; string path="../data/templeRing";


//string inputFileName = "bowl_par.txt"; string path="../data/bowl";
//string inputFileName = "venus_par.txt"; string path="../data/venus";
//string inputFileName = "pot_par.txt"; string path="../data/pot";

string scaleX = "";
int scaleXd = 1;

int patchGridRes = 4;

//string outputOBJ ="dinoSparseRing.obj";

int daisyParas[] = {15,3,8,4};
int daisyParasNormals[] = {15,3,8,4};

int daisyLength = (daisyParas[1]*daisyParas[2]+1)*daisyParas[3];
int daisyNormalsLength = (daisyParasNormals[1]*daisyParasNormals[2]+1)*daisyParasNormals[3];

const int patchSize=5;
const int patchSize2=patchSize/2;

#define theta_lower 5.0
#define theta_upper 45.0

#define distance_lower 0.05
#define distance_upper 2

#define top_K 10
#define sampling_step 2
#define NCC .8
#define normalNCC .8
#define beta 3
#define gamma 2
#define angleSearchSpace 40
#define angleSearchStep 5

#define angleThreshold .5

#define isObject true

#define deltaX 1e-9

#define buildTrackDepth 3
#define DQPI     3.141592653589793

//#define BRDFFITTING

#endif