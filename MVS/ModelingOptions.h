#ifndef __MODELINGOPTIONS_H__
#define __MODELINGOPTIONS_H__

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <limits>
#include <fstream>
#include <cassert>
#include <ctime>
#include <cmath>
#include <map>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/function.hpp>
#include <boost/random.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include <boost/make_shared.hpp>
#include <boost/thread.hpp>
#include <boost/throw_exception.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

typedef boost::numeric::ublas::bounded_matrix<double, 4, 4>Matrix_4D;
typedef boost::numeric::ublas::bounded_matrix<float, 4, 4>	Matrix_4F;
typedef boost::numeric::ublas::bounded_vector<double, 4>	Vector_4D;
typedef boost::numeric::ublas::bounded_vector<float, 4>	Vector_4F;
typedef boost::numeric::ublas::bounded_matrix<double, 3, 3>Matrix_3D;
typedef boost::numeric::ublas::bounded_matrix<float, 3, 3>	Matrix_3F;
typedef boost::numeric::ublas::bounded_vector<double, 3>	Vector_3D;
typedef boost::numeric::ublas::bounded_vector<float, 3>	Vector_3F;

typedef Vector_3D MyColor;

void CompleteDirName(string &dirname)
{
	if (dirname[dirname.length()-1] != '\\') {
		dirname = dirname + "\\";
	}
}

class ModelingOptions
{
public:
	ModelingOptions() {
		DirName				= ".\\";
		FileTypes			= "jpg";
		MaskDirName			= ".\\masks";
		MaskFileTypes		= "bmp";
		list_file			= "list.txt";
		bundle_dir_name		= "bundle\\";
		bundle_file			= "bundle.out";
		pmvs_dir_name		= "pmvs\\";
		pmvs_file			= "pmvs_result";
		middlebury_dir_name = "MiddleBury\\";
		middlebury_par_file = "MiddleBury_par.txt";
		brdf_db_name		= DirName + "brdfdb";

		range_value			= 0.5;
		fitting_choice		= 0;

		fidParam			= 0;
		m_beta				= 0;
		pcd_eta				= 0;
		fcd_eta				= 0;
		fnd_eta				= 0;
		varsigma			= 20;
		ALM_TVNorm			= false;
		ALM_TVU				= false;
		penParam			= 0.01;
		regParam			= 1.0;
		pcParam				= 0.0;
		lapParam			= 0;

		Is_Matched			= true;
		Is_Bundler			= true;
		Using_PMVS			= true;
		Using_BRDF			= true;
		Is_PoissonRec		= false;
		UseFaceArea			= false;

		UseFittedColor		= false;
		UseMatlabSolver		= true;
		ScaleDelta			= true;
		AnisotropicLaplace	= true;
	};
	bool ReadOptionFile(const char* optionfilename) {
		//Initialization from the input file name
		ifstream ifstr;
		ifstr.open(optionfilename);
		if (!ifstr.is_open()) {
			return false;
		}
		while (1) {
			string name;
			ifstr >> name;
			if (ifstr.eof()) {
				break;
			}
			if (name[0] == '#') {
				char buffer[1024];
				ifstr.putback('#');
				ifstr.getline(buffer, 1024);
				continue;
			}
			if (name.substr(name.find_last_of("-")+1) == "imagedir" || 
				name.substr(name.find_last_of("-")+1) == "idir") {
					ifstr >> DirName; 
					CompleteDirName(DirName);
					brdf_db_name	= DirName + "brdfdb";
					CompleteDirName(brdf_db_name);
					MaskDirName		= DirName + "mask";
					CompleteDirName(MaskDirName);
			}
			if (name.substr(name.find_last_of("-")+1) == "filetype" ||
				name.substr(name.find_last_of("-")+1) == "ft") {
					ifstr >> FileTypes;
			}
			if (name.substr(name.find_last_of("-")+1) == "maskdir" || 
				name.substr(name.find_last_of("-")+1) == "mdir") {
					ifstr >> MaskDirName; 
					CompleteDirName(MaskDirName);
			}
			if (name.substr(name.find_last_of("-")+1) == "maskfiletype" ||
				name.substr(name.find_last_of("-")+1) == "mft") {
					ifstr >> MaskFileTypes;
			}
			if (name.substr(name.find_last_of("-")+1) == "method" ||
				name.substr(name.find_last_of("-")+1) == "me") {
					string upmethod; ifstr >> upmethod;
					transform(upmethod.begin(), upmethod.end(),upmethod.begin(), toupper);
					if (upmethod == "PMVS") {
						Using_PMVS = true; Using_BRDF = false;
					}
					if (upmethod == "BRDF") {
						Using_PMVS = false; Using_BRDF = true;
					}
			}
			if (name.substr(name.find_last_of("-")+1) == "brdfdb" ||
				name.substr(name.find_last_of("-")+1) == "bdb") {
					ifstr >> brdf_db_name;
					CompleteDirName(brdf_db_name);
			}

			if (name.substr(name.find_last_of("-")+1) == "RAV" ||
				name.substr(name.find_last_of("-")+1) == "rav") {
					ifstr >> range_value;
			}
			if (name.substr(name.find_last_of("-")+1) == "FC" ||
				name.substr(name.find_last_of("-")+1) == "fc") {
					ifstr >> fitting_choice;
			}

			if (name.substr(name.find_last_of("-")+1) == "FID" ||
				name.substr(name.find_last_of("-")+1) == "fid") {
					ifstr >> fidParam;
			}
			if (name.substr(name.find_last_of("-")+1) == "BETA" ||
				name.substr(name.find_last_of("-")+1) == "beta") {
					ifstr >> m_beta;
			}
			if (name.substr(name.find_last_of("-")+1) == "PC" ||
				name.substr(name.find_last_of("-")+1) == "pc") {
					ifstr >> pcParam; 
			}
			if (name.substr(name.find_last_of("-")+1) == "FCD" ||
				name.substr(name.find_last_of("-")+1) == "fcd") {
					ifstr >> fcd_eta; 
			}
			if (name.substr(name.find_last_of("-")+1) == "FND" ||
				name.substr(name.find_last_of("-")+1) == "fnd") {
					ifstr >> fnd_eta; 
			}
			if (name.substr(name.find_last_of("-")+1) == "PCD" ||
				name.substr(name.find_last_of("-")+1) == "pcd") {
					ifstr >> pcd_eta; 
			}
			if (name.substr(name.find_last_of("-")+1) == "VAR" ||
				name.substr(name.find_last_of("-")+1) == "var") {
					ifstr >> varsigma;
			}
			if (name.substr(name.find_last_of("-")+1) == "TVN" ||
				name.substr(name.find_last_of("-")+1) == "tvn") {
					ifstr >> penParam; ALM_TVNorm = true;
			}
			if (name.substr(name.find_last_of("-")+1) == "TVU" ||
				name.substr(name.find_last_of("-")+1) == "tvu") {
					ifstr >> penParam; ALM_TVU = true;
			}
			if (name.substr(name.find_last_of("-")+1) == "REG" ||
				name.substr(name.find_last_of("-")+1) == "reg") {
					ifstr >> regParam; 
			}
			if (name.substr(name.find_last_of("-")+1) == "LAP" ||
				name.substr(name.find_last_of("-")+1) == "lap") {
					ifstr >> lapParam; 
			}
			if (name.substr(name.find_last_of("-")+1) == "FA" ||
				name.substr(name.find_last_of("-")+1) == "fa") {
					int temp; ifstr >> temp; UseFaceArea = temp>0 ? true:false; 
			}

			if (name.substr(name.find_last_of("-")+1) == "FITCOLOR" ||
				name.substr(name.find_last_of("-")+1) == "fitcolor") {
					int temp; ifstr >> temp; UseFittedColor = temp>0 ? true:false; 
			}
			if (name.substr(name.find_last_of("-")+1) == "MAT" ||
				name.substr(name.find_last_of("-")+1) == "mat") {
					int temp; ifstr >> temp; UseMatlabSolver = temp>0 ? true:false; 
			}
			if (name.substr(name.find_last_of("-")+1) == "SCD" ||
				name.substr(name.find_last_of("-")+1) == "scd") {
					int temp; ifstr >> temp; ScaleDelta = temp>0 ? true:false; 
			}
			if (name.substr(name.find_last_of("-")+1) == "AL" ||
				name.substr(name.find_last_of("-")+1) == "al") {
					int temp; ifstr >> temp; AnisotropicLaplace = temp>0 ? true:false; 
			}
		}
		ifstr.close();
		return true;
	}
	string	DirName;
	string	FileTypes;
	string	MaskDirName;
	string  MaskFileTypes;
	string	list_file;
	string	bundle_dir_name;
	string	bundle_file;
	string	pmvs_dir_name;
	string	pmvs_file;
	string  middlebury_dir_name;
	string  middlebury_par_file;
	string	brdf_db_name;

	double  range_value;
	int		fitting_choice;
	string	meshname;

	double	fidParam;
	double	m_beta;
	double	pcd_eta;
	double	fcd_eta;
	double	fnd_eta;
	double	varsigma;
	bool	ALM_TVNorm;
	bool	ALM_TVU;
	double	penParam;
	double	regParam;
	double	pcParam;
	double  lapParam;

	bool	Is_Matched;
	bool	Is_Bundler;
	bool	Is_PoissonRec;
	bool	Using_PMVS;
	bool	Using_BRDF;
	bool	UseFaceArea;

	bool	UseFittedColor;
	bool	UseMatlabSolver;
	bool	ScaleDelta;
	bool	AnisotropicLaplace;
};

extern ModelingOptions MOptions;

void ShowUsage()
{
	std::cout
		<<"BRDF Based 3D Model Reconstruction Usage:\n"
		<<"-h -help							: Parameter information\n"
		<<"-optionfile (-of) <string>		: Option file with options inside, default: option.txt\n"
		<<"-imagedir (-idir) <strings>		: The directory name of input image sources, default: current directory \n"
		<<"-filetype (-ft) <string>			: File type to be processed in the database, default: jpg\n"
		<<"-maskdir (-mdir) <strings>		: The directory name of input mask images, default: current directory\mask \n"
		<<"-maskfiletype (-mft) <string>	: File type to be processed in the database, default: bmp\n"
		<<"-method (-me) <string>			: Reconstruction method, PMVS or BRDF, default: BRDF\n"
		<<"-brdfdb (-bdb) <string>			: The directory of brdf database, default: 'imagedir'\brdfdb \n"
		//<<"-o <string>       : Where to save SIFT features\n"
		//<<"-f <float>        : Filter width factor; Width will be 2*factor+1 (default : 4.0)\n"
		//<<"-w  <float>       : Orientation sample window factor (default: 2.0)\n"
		//<<"-dw <float>  *    : Descriptor grid size factor (default : 3.0)\n"
		//<<"-fo <int>    *    : First octave to detect DOG keypoints(default : 0)\n"
		//<<"-no <int>         : Maximum number of Octaves (default : no limit)\n"
		//<<"-d <int>          : Number of DOG levels in an octave (default : 3)\n"
		//<<"-t <float>        : DOG threshold (default : 0.02/3)\n"
		//<<"-e <float>        : Edge Threshold (default : 10.0)\n"
		//<<"-m  <int=2>       : Multi Feature Orientations (default : 1)\n"
		//<<"-m2p              : 2 Orientations packed as one float\n"
		//<<"-s  <int=1>       : Sub-Pixel, Sub-Scale Localization, Multi-Refinement(num)\n"
		//<<"-lcpu -lc <int>   : CPU/GPU mixed Feature List Generation (defaut : 6)\n"
		//<<"                    Use GPU first, and use CPU when reduction size <= pow(2,num)\n"
		//<<"                    When <num> is missing or equals -1, no GPU will be used\n"
		//<<"-noprep           : Upload raw data to GPU (default: RGB->LUM and down-sample on CPU)\n"
		//<<"-sd               : Skip descriptor computation if specified\n"
		//<<"-unn    *         : Write unnormalized descriptor if specified\n"
		//<<"-b      *         : Write binary sift file if specified\n"
		//<<"-fs <int>         : Block Size for freature storage <default : 4>\n"
		//<<"-cuda <int=0>     : Use CUDA SiftGPU, and specifiy the device index\n"
		//<<"-tight            : Automatically resize pyramid to fit new images tightly\n"
		//<<"-p  <W>x<H>       : Inititialize the pyramids to contain image of WxH (eg -p 1024x768)\n"
		//<<"-tc[1|2|3] <int> *: Threshold for limiting the overall number of features (3 methods)\n"
		//<<"-v <int>          : Level of timing details. Same as calling Setverbose() function\n"
		//<<"-loweo            : (0, 0) at center of top-left pixel (defaut: corner)\n"
		//<<"-maxd <int> *     : Max working dimension (default : 2560 (unpacked) / 3200 (packed))\n"
		//<<"-nomc             : Disabling auto-downsamping that try to fit GPU memory cap\n"
		//<<"-exit             : Exit program after processing the input image\n"
		//<<"-unpack           : Use the old unpacked implementation\n"
		//<<"-di               : Use dynamic array indexing if available (defualt : no)\n"
		//<<"                    It could make computation faster on cards like GTX 280\n"
		//<<"-ofix     *       : use 0 as feature orientations.\n"
		//<<"-ofix-not *       : disable -ofix.\n"
		//<<"-winpos <X>x<Y> * : Screen coordinate used in Win32 to select monitor/GPU.\n"
		//<<"-display <string>*: Display name used in Linux/Mac to select monitor/GPU.\n"
		<<"\n"
		<<"NOTE: parameters marked with * can be changed after initialization\n"
		<<"\n";
}

void ParseParam(int argc, char* argv[], ModelingOptions &options)
{
	if (argc < 2 && !options.ReadOptionFile("option.txt")) {
		ShowUsage();
		return;
	}
	for (int i = 1; i < argc; ++i) {
		string arg = string(argv[i]);
		if(arg.length() <1 || arg[0] != '-')continue;
		string opt = arg.substr(arg.find_first_of("-")+1);
		string param = string(argv[i+1]);
		if (opt == "optionfile" || opt == "of") {
			options.ReadOptionFile(param.c_str());
		}
		if (opt == "imagedir" || opt == "idir") {
			options.DirName			= param;	i++;
			CompleteDirName(options.DirName);
			options.brdf_db_name	= options.DirName + "brdfdb";
			CompleteDirName(options.brdf_db_name);
			options.MaskDirName		= options.DirName + "masks";
		}
		if (opt == "filetype" || opt == "ft") {
			options.FileTypes		= param;	i++;
		}
		if (opt == "maskdir" || opt == "mdir") {
			options.MaskDirName		= param;	i++;
			CompleteDirName(options.MaskDirName);
		}
		if (opt == "maskfiletype" || opt == "mft") {
			options.MaskFileTypes	= param;	i++;
		}
		if (opt == "method" || opt == "me") {
			string upmethod = param;
			transform(upmethod.begin(), upmethod.end(),upmethod.begin(), toupper);
			if (upmethod == "PMVS") {
				options.Using_PMVS = true; options.Using_BRDF = false;
			}
			if (upmethod == "BRDF") {
				options.Using_PMVS = false; options.Using_BRDF = true;
			}	
			i++;
		}
		if (opt == "brdfdb" || opt == "bdb") {
			options.brdf_db_name	= param;	i++;
			CompleteDirName(options.brdf_db_name);
		}
		if (opt == "RAV" || opt == "rav") {
			options.range_value	= atof(param.c_str());	i++;
		}
		if (opt == "FC" || opt == "fc") {
			options.fitting_choice	= atoi(param.c_str());	i++;
		}

		if (opt == "FID" || opt == "fid") {
			options.fidParam	= atof(param.c_str());	i++;
		}
		if (opt == "BETA" || opt == "beta") {
			options.m_beta	= atof(param.c_str());	i++;
		}
		if (opt == "PC" || opt == "pc") {
			options.pcParam	= atof(param.c_str());	i++;
		}
		if (opt == "FCD" || opt == "fcd") {
			options.fcd_eta	= atof(param.c_str());	i++;
		}
		if (opt == "FND" || opt == "fnd") {
			options.fnd_eta	= atof(param.c_str());	i++;
		}
		if (opt == "PCD" || opt == "pcd") {
			options.pcd_eta	= atof(param.c_str());	i++;
		}
		if (opt == "VAR" || opt == "var") {
			options.varsigma	= atof(param.c_str());	i++;
		}
		if (opt == "TVU" || opt == "tvu") {
			options.ALM_TVU = true;
			options.penParam	= atof(param.c_str());	i++;
		}
		if (opt == "TVN" || opt == "tvn") {
			options.ALM_TVNorm = true;
			options.penParam	= atof(param.c_str());	i++;
		}
		if (opt == "REG" || opt == "reg") {
			options.regParam	= atof(param.c_str());	i++;
		}
		if (opt == "LAP" || opt == "lap") {
			options.lapParam	= atof(param.c_str());	i++;
		}
		if (opt == "FA" || opt == "fa") {
			int temp = atoi(param.c_str());
			options.UseFaceArea	= temp>0?true:false;	i++;
		}
		if (opt == "FITCOLOR" || opt == "fitcolor") {
			int temp = atoi(param.c_str());
			options.UseFittedColor	= temp>0?true:false;	i++;
		}
		if (opt == "MAT" || opt == "mat") {
			int temp = atoi(param.c_str());
			options.UseMatlabSolver	= temp>0?true:false;	i++;
		}
		if (opt == "SCD" || opt == "scd") {
			int temp = atoi(param.c_str());
			options.ScaleDelta	= temp>0?true:false;	i++;
		}
		if (opt == "AL" || opt == "al") {
			int temp = atoi(param.c_str());
			options.AnisotropicLaplace	= temp>0?true:false;	i++;
		}

		if (opt == "h" || opt == "help") {
			ShowUsage();
		}
	}
}

void MyCreateDirectory(string dir_name)
{
	boost::filesystem::path dir_path(dir_name);
	if (!boost::filesystem::exists(dir_path) || !is_directory(dir_path)) {
		boost::filesystem::create_directory(dir_path);
	}
}

void GetDirFileList(string DirName, vector<string> &filelist, string filetype = "jpg")
{
	string upfiletype = filetype; transform(upfiletype.begin(), upfiletype.end(),upfiletype.begin(),toupper);
	string lowfiletype = filetype; transform(lowfiletype.begin(), lowfiletype.end(),lowfiletype.begin(),tolower);
	boost::filesystem::path dir_path(DirName);
	if (!boost::filesystem::exists(dir_path)) {
		cout << DirName << " is not exist ..." << endl;
		return;
	}
	if (!is_directory(dir_path)) {
		cout << DirName << " is not a directory ..." << endl;
		return;
	}
	boost::filesystem::directory_iterator end_iter;
	for (boost::filesystem::directory_iterator iter(dir_path); iter != end_iter; ++iter) {
		string fn = iter->path().string();
		if (filetype.length() < 1 || 
			(filetype.length() > 1 && fn.substr(fn.find_last_of(".")+1) == upfiletype) ||
			(filetype.length() > 1 && fn.substr(fn.find_last_of(".")+1) == lowfiletype) ) {
				//filelist.push_back(fn.substr(fn.find_last_of("\\")+1));
				filelist.push_back(fn);
		}
	}
}


void CopyDirectory(string Src_dir, string Dst_dir)
{
	boost::filesystem::path src(Src_dir); 
	boost::filesystem::path dst(Dst_dir);
	if (! boost::filesystem::exists(dst)) {  
		boost::filesystem::create_directories(dst);  
	} 
	for (boost::filesystem::directory_iterator it(src); it != boost::filesystem::directory_iterator(); ++it) {
		const boost::filesystem::path newSrc = src / it->path().string();  
		const boost::filesystem::path newDst = dst / it->path().string();  
		if (boost::filesystem::is_directory(newSrc)) {  
			CopyDirectory(newSrc.string(), newDst.string());  
		}  
		else if (boost::filesystem::is_regular_file(newSrc)) {  
			boost::filesystem::copy_file(newSrc, newDst, boost::filesystem::copy_option::overwrite_if_exists);  
		}  
		else {
			fprintf(stderr, "Error: unrecognized file - %s", newSrc.string().c_str());  
		} 
	}
}

void MyCopyFile(string src_filename, string dst_filename)  
{  
	boost::filesystem::path src(src_filename); 
	boost::filesystem::path dst(dst_filename);
	if (! boost::filesystem::exists(src)) {  
		cout << src_filename << " is not existed" << endl; 
		return;  
	}  
	if (boost::filesystem::is_regular_file(src)) {
		boost::filesystem::copy_file(src, dst, boost::filesystem::copy_option::overwrite_if_exists); 
	} 
} 

void MyDeleteFile(string filename)
{
	boost::filesystem::path src(filename);
	if (!boost::filesystem::exists(src)) {
		cout << filename << " is not existed" << endl; 
		return; 
	}
	boost::filesystem::remove(src);
}

void MyMoveFile(string src_filename, string dst_filename)
{
	MyCopyFile(src_filename, dst_filename);
	MyDeleteFile(src_filename);
}

bool FileExisted(const char* name)
{
	FILE * myfile;
	myfile = fopen (name,"r");

	if (myfile == NULL) 
	{
		//cout << name << " is not existed...." << endl;
		return false;
	}
	fclose(myfile);
	return true;
}

#endif //__MODELINGOPTIONS_H__