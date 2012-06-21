#ifndef _VIEWPOINT_HPP
#define _VIEWPOINT_HPP

#include "rectify.h"

class viewPoint
{
public:
	int height;
	int width;

	int ViewID;
	string imageName;
	string maskName;
	IplImage* image;  IplImage* m_Image;
	IplImage* xImage;
	uchar *ucharImage;
	uchar *ucharXImage;
	bool **sil;
	double K[9],Kinv[9]; // intrinsic matrix
	double R[9]; // rotation matrix
	double T[3]; // translation vector
	double C[3]; // camera center
	double P[12]; // projection matrix
	double M[9]; // M = K x R
	double RT[12]; // [R T]
	double viewingVector[3]; //  transpose (M[3]) This vector start from camera position to optical center
	double InverseViewVector[3];//This direction is actually the view direction, from optical center to camera position
	int stereoPairs[2];
	daisy *desc_original;
	bool verticalStereo[2];
	int **originalMatchingResultX[2];
	int **originalMatchingResultY[2];
	bool **collected;

	vector<int> visibletrackid;

	cv::Mat LightMatrix; 
	double LightDirection[3];
	double LightLuminance[3];

	viewPoint()
	{
		image = xImage = m_Image = NULL;
		ucharImage = ucharXImage = NULL;
		sil = NULL;desc_original = NULL;
		originalMatchingResultX[0] = originalMatchingResultX[1] = NULL; 
		originalMatchingResultY[0] = originalMatchingResultY[1] = NULL;
		collected = NULL;
	}

	~viewPoint()
	{
		if (sil != NULL)									deallocate(sil,height);
		if (xImage != NULL && xImage != image)				cvReleaseImage(&xImage);
		if (image != NULL)									cvReleaseImage(&image);
		if (m_Image != NULL)								cvReleaseImage(&m_Image);
		if (ucharXImage != NULL&&ucharXImage!=ucharImage)	delete ucharXImage;
		if (ucharImage != NULL)								delete ucharImage;
		if (desc_original != NULL)							delete desc_original;
		if (collected != NULL)								deallocate(collected,height);
		if (originalMatchingResultX[0] != NULL)				deallocate(originalMatchingResultX[0],height);
		if (originalMatchingResultY[0] != NULL)				deallocate(originalMatchingResultY[0],height);
		if (originalMatchingResultX[1] != NULL)				deallocate(originalMatchingResultX[1],height);
		if (originalMatchingResultY[1] != NULL)				deallocate(originalMatchingResultY[1],height);
	}

	void calP()
	{
		for (int i=0;i<3;i++)
		{
			for (int j=0;j<3;j++) 
				RT[i*4+j] = R[i*3+j];

			RT[i*4+3] = T[i];
		}
		matrixMul(K,3,3,RT,3,4,P);
	}


	void calC()
	{

		// P : calibration matrix
		// P = KR [ I | -C ]
		// P = [KR | KT]
		// KR [ I | -C] = [KR | KT]
		// [KR | KR-C] = [KR | KT]
		// KR(-C) = KT
		// R(-C) = T
		// C = - inv(R) * T
		// Multiple View Geometry in computer vision 
		// p 156

		cv::Mat _R = cv::Mat(3,3,CV_64F, R);
		cv::Mat _T = cv::Mat(3,1,CV_64F, T);

		cv::Mat _C = - _R.inv() * _T;
		C[0] = _C.at<double>(0,0);
		C[1] = _C.at<double>(1,0);
		C[2] = _C.at<double>(2,0);


	}

	void calM()
	{	
		// M = K x R
		// P = [M KT]
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++) M[i*3+j] = P[i*4+j];
	}

	void calKinv()
	{
		cv::Mat _Kinv = cv::Mat(3,3,CV_64F,K).inv();
		for (int i=0;i<9;i++) Kinv[i] = _Kinv.at<double>(i/3,i%3);
	}

	void calViewingVector()
	{		
		cv::Mat _M = cv::Mat(3,3,CV_64F,M);
		double detM = determinant(_M);
		//printf("determiant M = %lf\n", detM);
		detM /= fabs(detM);
		for (int i=0;i<3;i++) {
			viewingVector[i] = detM*M[6+i];
			InverseViewVector[i] = -viewingVector[i];
		}
		matrixNorm(viewingVector,3); matrixNorm(InverseViewVector, 3);
	}

	viewPoint(int id, string filename, double _K[], double _R[], double _T[])
	{
		ViewID = id;
		imageName = filename;
		matrixCopy(_K,K,9);
		matrixCopy(_R,R,9);
		matrixCopy(_T,T,3);

		image = cvLoadImage(imageName.c_str(),0); m_Image = cvLoadImage(imageName.c_str(),1);
		height = image->height;
		width = image->width;
		//printf("Reading image %s... Height %d x Width %d\n", imageName.c_str(), height, width);		
		calP();
		calM();

		calC();
		//printf("center = %10.3lf %10.3lf %10.3lf\n",C[0],C[1],C[2]);

		calViewingVector();
		//printf("viewing vector = %10.3lf %10.3lf %10.3lf\n",viewingVector[0],viewingVector[1],viewingVector[2]);

		calKinv();
		stereoPairs[0] = stereoPairs[1] = -1;		
		verticalStereo[0] = verticalStereo[1] = false;
		ucharImage = new uchar[height * width];
		for (int y=0;y<height;y++)
		{
			for (int x=0;x<width;x++)
			{
				ucharImage[y*width+x] = (uchar)cvGet2D(image,y,x).val[0];
				//if (ucharImage[y*width+x]==0) printf("0"); else printf("1");
			}
			//printf("\n");
		}

		if (scaleX.length()>0)
		{
			xImage = cvLoadImage((path+(string)"/"+scaleX+imageName.substr(imageName.find_last_of("\\")+1)).c_str(),0);
			ucharXImage = new uchar[xImage->height*xImage->width];
			for (int y=0;y<xImage->height;y++)
			{
				for (int x=0;x<xImage->width;x++)
				{
					ucharXImage[y*xImage->width+x] = (uchar)cvGet2D(xImage,y,x).val[0];
				}
			}
		} else 
		{
			xImage = image;
			ucharXImage = ucharImage;
		}
	}

	viewPoint(int id, string filename, double _K[], double _R[], double _T[], double _P[])
	{
		ViewID = id;
		imageName = filename;
		matrixCopy(_K,K,9);
		matrixCopy(_R,R,9);
		matrixCopy(_T,T,3);
		matrixCopy(_P,P,12);

		image = cvLoadImage(imageName.c_str(),0); m_Image = cvLoadImage(imageName.c_str(),1);
		height = image->height;
		width = image->width;
		//printf("Reading image %s... Height %d x Width %d\n", imageName.c_str(), height, width);		

		calM();

		calC();
		//printf("center = %10.3lf %10.3lf %10.3lf\n",C[0],C[1],C[2]);

		calViewingVector();
		//printf("viewing vector = %10.3lf %10.3lf %10.3lf\n",viewingVector[0],viewingVector[1],viewingVector[2]);

		calKinv();

		stereoPairs[0] = stereoPairs[1] = -1;
		verticalStereo[0] = verticalStereo[1] = false;
		ucharImage = new uchar[height * width];
		for (int y=0;y<height;y++)
			for (int x=0;x<width;x++)
			{
				ucharImage[y*width+x] = (uchar)cvGet2D(image,y,x).val[0];
				//printf("%d ", im[y*w+x]);
			}
	}

	void readSilhouette(const char *name)
	{		
		maskName = string(name);
		if (name!=NULL && FileExisted(name))
		{
			//printf("Read mask file: %s\n",name);
			cv::Mat silImage = cv::Mat(cvLoadImage(name,0));
			double back_color = silImage.at<unsigned char>(1,1);
			sil = new bool* [silImage.rows];
			for (int x=0;x<silImage.rows;x++)
			{
				sil[x] = new bool[silImage.cols];
				for (int y=0;y<silImage.cols;y++)
				{
					//printf("%d ",silImage.at<unsigned char>(x,y));
					if (abs(silImage.at<unsigned char>(x,y)-back_color)< 100) 
						sil[x][y] = false; 
					else 
						sil[x][y] = true;
				}
			}
		} else
		{

			sil = allocate<bool>(height,width);
			for (int x=0;x<height;x++) 
				for (int y=0;y<width;y++) sil[x][y] = true;

		}
	}


	void rectifyImages(int indexPair, viewPoint* stereoViewPoint, 
		CvMat* &mx1, CvMat* &my1, CvMat* &mx2, CvMat* &my2, CvMat* &img1r, CvMat* &img2r)
	{

		CvSize imageSize = cvSize(width, height);

		double K1[9], K2[9];
		for (int i=0;i<9;i++) K1[i] = K[i];
		for (int i=0;i<9;i++) K2[i] = stereoViewPoint->K[i];
		CvMat _K1 = cvMat(3,3, CV_64F, K1);
		CvMat _K2 = cvMat(3,3, CV_64F, K2);	

		double R1[9],R2[9];
		for (int i=0;i<9;i++) R1[i] = R[i];
		for (int i=0;i<9;i++) R2[i] = stereoViewPoint->R[i];
		CvMat _R1 = cvMat(3, 3, CV_64F,  R1);
		CvMat _R2 = cvMat(3, 3, CV_64F,  R2);

		double T1[3],T2[3];
		for (int i=0;i<3;i++) T1[i] = T[i];
		for (int i=0;i<3;i++) T2[i] = stereoViewPoint->T[i];

		CvMat _T1 = cvMat(3, 1, CV_64F,  T1);
		CvMat _T2 = cvMat(3, 1, CV_64F,  T2);

		double _K_l_rect[9] = {0}, _K_r_rect[9] = {0}, _R_l_rect[9] = {0}, _R_r_rect[9] = {0}, _T_New[3] = {0};

		CvMat K_l_rect = cvMat(3, 3, CV_64F, _K_l_rect);
		CvMat K_r_rect = cvMat(3, 3, CV_64F, _K_r_rect);
		CvMat R_l_rect = cvMat(3, 3, CV_64F, _R_l_rect);
		CvMat R_r_rect = cvMat(3, 3, CV_64F, _R_r_rect);
		CvMat T_New = cvMat(3, 1, CV_64F, _T_New);

		int stereo_type;	//horizontal or vertical epipolar line

		icvRectify(&_K1, &_K2, &_R1, &_R2, &_T1, &_T2, 
			NULL, NULL, imageSize, 
			&K_l_rect, &K_r_rect, &R_l_rect, &R_r_rect, &T_New, stereo_type);

		bool isVerticalStereo = (stereo_type == VER_STEREO);

		mx1 = cvCreateMat( imageSize.height,imageSize.width, CV_32F );
		my1 = cvCreateMat( imageSize.height,imageSize.width, CV_32F );
		mx2 = cvCreateMat( imageSize.height,imageSize.width, CV_32F );
		my2 = cvCreateMat( imageSize.height,imageSize.width, CV_32F );
		img1r = cvCreateMat( imageSize.height, imageSize.width, CV_8U );
		img2r = cvCreateMat( imageSize.height, imageSize.width, CV_8U );

		double distort[4] = {0};
		CvMat D1 = cvMat(4,1,CV_64F, distort);

		double _R_rected_to_org[9];
		CvMat R_rected_to_org = cvMat(3, 3, CV_64F, _R_rected_to_org);
		cvTranspose(&R_l_rect, &R_rected_to_org);

		double baseline;
		double newf = cvmGet(&K_l_rect, 0, 0);		
		if( fabs(cvmGet(&T_New, 0, 0)) > fabs(cvmGet(&T_New, 1, 0)) )
		{
			baseline = cvmGet(&T_New, 0, 0);
			cvmSet(&T_New, 1, 0, 0);
			cvmSet(&T_New, 2, 0, 0);
		}
		else
		{
			baseline = cvmGet(&T_New, 1, 0);			
			cvmSet(&T_New, 0, 0, 0);
			cvmSet(&T_New, 2, 0, 0);
		}

		//Precompute maps for cvRemap()
		icvInitUndistortRectifyMap(&_K1, &D1, &R_l_rect, &K_l_rect, mx1, my1);
		icvInitUndistortRectifyMap(&_K2, &D1, &R_r_rect, &K_r_rect, mx2, my2);
		// after rectify, the rotation matrix is identity matrix, only translation existed


		IplImage* img1 = image;
		IplImage* img2 = stereoViewPoint->image;


		if( img1 && img2 )
		{
			//CvMat part;
			// dst(x,y) <= src(mapx(x,y),mapy(x,y))
			cvRemap( img1, img1r, mx1, my1 );
			cvRemap( img2, img2r, mx2, my2 );

			//cvSaveImage("a.png", img1r);
			//cvSaveImage("b.png", img2r);
			/*

			CvMat* pair;
			if( !isVerticalStereo )
			pair = cvCreateMat( imageSize.height, imageSize.width*2, CV_8UC1 );
			else
			pair = cvCreateMat( imageSize.height*2, imageSize.width, CV_8UC1 );
			//Setup for finding stereo corrrespondences

			//cvNamedWindow( "rectified", 1 );
			if( !isVerticalStereo )
			{
			cvGetCols( pair, &part, 0, imageSize.width );
			cvConvert( img1r, &part);
			cvGetCols( pair, &part, imageSize.width, imageSize.width*2 );
			cvConvert( img2r, &part);
			for(int j = 0; j < imageSize.height; j += 16 )
			cvLine( pair, cvPoint(0,j), cvPoint(imageSize.width*2,j), CV_RGB(0,0,0));
			}
			else
			{
			cvGetRows( pair, &part, 0, imageSize.height );
			cvConvert( img1r, &part);
			cvGetRows( pair, &part, imageSize.height, imageSize.height*2 );
			cvConvert( img2r, &part);
			for(int j = 0; j < imageSize.width; j += 16 )
			cvLine( pair, cvPoint(j,0), cvPoint(j,imageSize.height*2), CV_RGB(0,0,0));
			}
			*/

			//cvShowImage( "rectified", pair );
			//cvWaitKey();
			//cvSaveImage("c.png", pair);
		}

		// free memory
		//cvReleaseImage( &img1 );
		//cvReleaseImage( &img2 );

		//cvReleaseMat( &mx1 );
		//cvReleaseMat( &my1 );
		//cvReleaseMat( &mx2 );
		//cvReleaseMat( &my2 );
		//cvReleaseMat( &img1r );
		//cvReleaseMat( &img2r );
		verticalStereo[indexPair] = isVerticalStereo;
	}



	void stereoMatching(int indexPair, daisy* desc_original_rectified, daisy *desc_paired_rectified, int** stereoMatchingResult)
	{

		for (int r=0;r<height;r++)
			for (int c=0;c<width;c++) stereoMatchingResult[r][c] = -1;

		if (verticalStereo[indexPair]) 
		{
			// vertically rectified
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
			for (int j=0;j<width;j++)
			{	

				int **matches = allocate<int>(height,top_K);				
				double **dp = allocate<double>(height,top_K);
				int **dp_res = allocate<int>(height,top_K);
				double **c = allocate<double>(height,top_K);				


				// kdtree method
				float* thor0 = NULL;
				float* thor1 = NULL;
				float *dataset = allocate<float>(height*daisyLength);
				for (int x=0;x<height;x++) 
				{					
					desc_paired_rectified->get_descriptor(x,j,thor0);
					for (int y=0;y<daisyLength;y++) dataset[x*daisyLength+y] = thor0[y];
				}
				cv::Mat features = cv::Mat(height, daisyLength, CV_32F, dataset);
				cv::flann::Index index(features, cv::flann::KDTreeIndexParams());	

				vector<float> query;
				vector<int> indices;
				vector<float> dists;
				query.resize(daisyLength);
				indices.resize(top_K);
				dists.resize(top_K);

				for (int i=0;i<height;i+=sampling_step)
				{					
					desc_original_rectified->get_descriptor(i,j,thor0);
					for (int k=0;k<daisyLength;k++) query[k] = thor0[k];
					index.knnSearch(query, indices, dists, top_K, cv::flann::SearchParams());
					for (int k1=0;k1<top_K;k1++)
					{
						for(int k2=k1+1;k2<top_K;k2++)
							if (indices[k1]>indices[k2])
							{
								int tmp = indices[k1];
								indices[k1] = indices[k2];
								indices[k2] = tmp;
							}
					}

					for (int k=0;k<top_K;k++) 
					{
						matches[i][k] = indices[k];						
						desc_paired_rectified->get_descriptor(indices[k],j,thor1);
						c[i][k] = cross_correlation(thor0, thor1, daisyLength);
						//printf("%f ", dists[k]);
					}				
					//printf("\n");

				}
				deallocate(dataset);
				//*/


				/*




				// linear method
				// find top_K candidates for each pixel in scan-line

				double **cc = allocate<double>(height,height);				
				for (int i=0;i<height;i+=sampling_step) 
				{
				vector<double> value;

				float* thor0 = NULL;
				desc_original_rectified->get_descriptor(i,j,thor0);
				for ( int ii = 0; ii < height ; ii ++ )
				{
				float* thor1 = NULL;
				desc_paired_rectified->get_descriptor(ii,j,thor1);
				cc[i][ii] = matrixL1(thor0, thor1, daisyLength);
				value.push_back(cc[i][ii]);
				}
				std::sort(value.begin(),value.end());
				double K_value = value[top_K];
				//printf("K_value = %lf	", K_value);
				int cntK = 0;
				for (int ii=0; ii<height; ii++)
				if (cc[i][ii]<=K_value && cntK < top_K )
				{
				matches[i][cntK] = ii;
				//printf("%d ",ii);
				float* thor1 = NULL;
				desc_paired_rectified->get_descriptor(ii,j,thor1);
				c[i][cntK] = cross_correlation(thor0, thor1, daisyLength);
				//if (c[i][cntK]<0.8) c[i][cntK] = 0.1;
				//c[i][ii] = 100 - c[i][ii];
				cntK++;
				}// else c[i][ii] = 0;
				//printf("\n");
				}

				//*/


				// DP part
				double tmp;
				for (int u=0; u<height; u+=sampling_step)
					for (int v=0; v<top_K; v++) 
					{
						dp[u][v] = c[u][v],  dp_res[u][v] = -1;
						if (u>=sampling_step)
						{
							for (int preV = 0; preV<top_K; preV++)
								if (matches[u][v]>matches[u-sampling_step][preV])
								{
									tmp = dp[u-sampling_step][preV] + c[u][v];
									if (dp[u][v]<tmp) 
									{
										// [u-sampling_step][preV] -> [u][v]
										dp[u][v] = tmp;
										dp_res[u][v] = preV;
									}
								} else if (dp[u][v]<dp[u-sampling_step][preV]) 
								{
									// not choose [u][v]
									dp[u][v] = dp[u-sampling_step][preV];
									dp_res[u][v] = preV;
								}
						}

						if (v>0 && dp[u][v]<dp[u][v-1])  
						{
							// not choose [u][v]
							dp[u][v] = dp[u][v-1];
							dp_res[u][v] = -2;
						}
					}				


					//printf("DP = %lf		", dp[h-1][top_K-1]);

					// retrieve dp result				
					int u = 0;
					while (u + sampling_step < height) u += sampling_step;
					int v = top_K-1;
					while (true)
					{
						int preV = dp_res[u][v];
						if (preV == -1) 
						{
							// only accept if cross-correlation >= cross_correlation_threshold
							if (c[u][v]>= NCC)
								stereoMatchingResult[u][j] = matches[u][v];
							break;
						}
						else if (preV == -2) v--;
						else 
						{
							if (matches[u][v]>matches[u-sampling_step][preV]) 
							{
								// only accept if cross-correlation >= cross_correlation_threshold
								if (c[u][v]>= NCC)
									stereoMatchingResult[u][j] = matches[u][v];						
							}
							u -= sampling_step;
							v = preV;						
						}
					}
					deallocate(matches,height);
					deallocate(c,height);
					deallocate(dp,height);
					deallocate(dp_res,height);
			}	

		}
		else 
		{

			// horizontally rectified
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
			for (int i=0;i<height;i++)
			{

				int **matches = allocate<int>(width,top_K);
				double **c = allocate<double>(width,top_K);
				double **dp = allocate<double>(width,top_K);
				int **dp_res = allocate<int>(width,top_K);


				///*
				// kdtree method
				float* thor0 = NULL;
				float* thor1 = NULL;
				float *dataset = allocate<float>(width*daisyLength);
				for (int x=0;x<width;x++) 
				{					
					desc_paired_rectified->get_descriptor(i,x,thor0);
					for (int y=0;y<daisyLength;y++) dataset[x*daisyLength+y] = thor0[y];
				}
				cv::Mat features = cv::Mat(width, daisyLength, CV_32F, dataset);
				cv::flann::Index index(features, cv::flann::KDTreeIndexParams());	
				vector<float> query;
				vector<int> indices;
				vector<float> dists;
				query.resize(daisyLength);
				indices.resize(top_K);
				dists.resize(top_K);

				for (int j=0;j<width;j+=sampling_step)
				{					
					desc_original_rectified->get_descriptor(i,j,thor0);
					for (int k=0;k<daisyLength;k++) query[k] = thor0[k];
					index.knnSearch(query, indices, dists, top_K, cv::flann::SearchParams());
					for (int k1=0;k1<top_K;k1++)
						for(int k2=k1+1;k2<top_K;k2++)
							if (indices[k1]>indices[k2])
							{
								int tmp = indices[k1];
								indices[k1] = indices[k2];
								indices[k2] = tmp;
							}

							for (int k=0;k<top_K;k++) 
							{
								matches[j][k] = indices[k];						
								desc_paired_rectified->get_descriptor(i,indices[k],thor1);
								c[j][k] = cross_correlation(thor0, thor1, daisyLength);
							}				

				}
				deallocate(dataset);
				//*/

				/*

				//linear method
				// find top_K candidates for each pixel in scan-line
				for (int j=0;j<width;j+=sampling_step) 
				{

				vector<double> value;

				float* thor0 = NULL;
				desc_original_rectified->get_descriptor(i,j,thor0);
				for ( int jj = 0; jj < width ; jj ++ )
				{
				float* thor1 = NULL;
				desc_paired_rectified->get_descriptor(i,jj,thor1);
				c[j][jj] = matrixL1(thor0, thor1, daisyLength);
				value.push_back(c[j][jj]);
				}
				std::sort(value.begin(),value.end());
				double K_value = value[top_K];
				//printf("K_value = %lf	", K_value);
				int cntK = 0;
				for (int jj=0;jj<width;jj++)
				if (c[j][jj]<=K_value && cntK < top_K )
				{
				matches[j][cntK++] = jj;
				//printf("%d ",ii);
				float* thor1 = NULL;
				desc_paired_rectified->get_descriptor(i,jj,thor1);
				c[j][jj] = cross_correlation(thor0, thor1, daisyLength);
				//c[i][ii] = 100 - c[i][ii];
				} else c[j][jj] = 0;
				//printf("\n");
				}
				*/


				// DP part
				double tmp;
				for (int u=0;u<width;u+=sampling_step)
					for (int v=0;v<top_K;v++) 
					{
						dp[u][v] = c[u][v],  dp_res[u][v] = -1;
						if (u>=sampling_step)
							for (int preV = 0; preV<top_K; preV++)
								if (matches[u][v]>matches[u-sampling_step][preV])
								{
									tmp = dp[u-sampling_step][preV] + c[u][v];
									if (dp[u][v]<tmp) dp[u][v] = tmp,  dp_res[u][v] = preV;
								} 
								else if (dp[u][v]<dp[u-sampling_step][preV]) dp[u][v] = dp[u-sampling_step][preV], dp_res[u][v] = preV;
								if (v>0 && dp[u][v]<dp[u][v-1])  dp[u][v] = dp[u][v-1], dp_res[u][v] = -2;
					}				


					//printf("DP = %lf		", dp[h-1][top_K-1]);

					// retreive dp result				
					int u = 0;
					while (u + sampling_step < width) u += sampling_step;
					int v = top_K-1;
					while (1)
					{
						int preV = dp_res[u][v];
						if (preV == -1) 
						{
							// only accept if cross-correlation >= cross_correlation_threshold
							if (c[u][v]>= NCC)
								stereoMatchingResult[i][u] = matches[u][v];
							break;
						}
						else if (preV == -2) v--;
						else if (matches[u][v]>matches[u-sampling_step][preV]) 
						{
							// only accept if cross-correlation >= cross_correlation_threshold
							if (c[u][v]>= NCC)
								stereoMatchingResult[i][u] = matches[u][v];
							u -= sampling_step;
							v = preV;
						}
						else v = preV, u-=sampling_step;
					}
					deallocate(matches,width);
					deallocate(c,width);
					deallocate(dp,width);
					deallocate(dp_res,width);
			}

		}		




	}

	void originalImageMatching(int indexPair, CvMat *mx1_map, CvMat *my1_map, CvMat *mx2_map, CvMat *my2_map,
		int **stereoMatchingResult, daisy *originalDaisy, daisy *rectifiedDaisy)
	{
		originalMatchingResultX[indexPair] = allocate<int>(height,width);
		originalMatchingResultY[indexPair] = allocate<int>(height,width);		
		double **c = allocate<double>(height,width);

		for (int i=0;i<height;i++)
			for (int j=0;j<width;j++)
			{
				originalMatchingResultX[indexPair][i][j] = -1;
				originalMatchingResultY[indexPair][i][j] = -1;

			}

			int cntError = 0;
			int fixedError = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
			for (int i=0;i<height;i++)
			{			
				for (int j=0;j<width;j++)
				{
					float *thor0, *thor1;				
					// original image point
					int i_original = (int)(cvmGet(my1_map,i,j)+0.5);
					int j_original = (int)(cvmGet(mx1_map,i,j)+0.5);
					if (matrixBound(i_original,0,height) && matrixBound(j_original,0,width))
					{
						if (stereoMatchingResult[i][j]>=0)
						{						
							int ii,jj;
							if (verticalStereo[indexPair])
							{
								ii = stereoMatchingResult[i][j];
								jj = j;
							} else
							{
								ii = i;
								jj = stereoMatchingResult[i][j];
							}
							// rectified image point
							int ii_original = (int)(cvmGet(my2_map,ii,jj)+0.5);
							int jj_original = (int)(cvmGet(mx2_map,ii,jj)+0.5);
							if (matrixBound(ii_original,0,height) && matrixBound(jj_original,0,width))
							{							
								originalDaisy->get_descriptor(i,j,thor0);								
								rectifiedDaisy->get_descriptor(ii,jj,thor1);
								//double _c = matrixL1(thor0,thor1,daisyLength);
								double _c = cross_correlation(thor0,thor1,daisyLength);
								//if (_c>cross_correlation_threshold)
								{

									if (originalMatchingResultX[indexPair][i_original][j_original]>=0)
									{

#pragma omp critical
										{
											cntError++;									
										}
										if (_c>c[i_original][j_original])
#pragma omp critical
										{
											//printf("%lf -> %lf \n", c[i_original][j_original],_c);
											fixedError ++;
											c[i_original][j_original] = _c;
											originalMatchingResultX[indexPair][i_original][j_original] = ii_original;
											originalMatchingResultY[indexPair][i_original][j_original] = jj_original;
										}



									} else
#pragma omp critical
									{
										originalMatchingResultX[indexPair][i_original][j_original] = ii_original;
										originalMatchingResultY[indexPair][i_original][j_original] = jj_original;	
										c[i_original][j_original] = _c;
									}
								}
							}
						}
					}
				}
			}
			printf("Total duplicated pixels = %d    -   changed = %d\n ", cntError, fixedError);
			deallocate(c, height);
	}

	bool SaveToFile(string filename)
	{
		//save the data into file
		std::fstream fout(filename.c_str(), std::ios::out);
		if (!fout) {
			cout << "Can not open " << filename << " to save data..." << endl;
			return false;
		}
		cout << "Start to save " << filename;
		// first save parameters
		fout << ViewID << endl;
		fout << imageName << endl;
		fout << maskName << endl;
		for (int j = 0; j < 9; ++j) {
			fout << K[j] << " ";
		}fout << endl;
		for (int j = 0; j < 9; ++j) {
			fout << R[j] << " ";
		}fout << endl;
		for (int j = 0; j < 3; ++j) {
			fout << T[j] << " ";
		}fout << endl;

		// save the matching result of mvs algorithm
		//fout << stereoPairs[0] << " " << stereoPairs[1] << endl;
		//fout << verticalStereo[0] << " " << verticalStereo[1] << endl;
		//for (int id = 0; id < 2; id ++) {
		//	for (int i=0;i<height;i++) {
		//		for (int j=0;j<width;j++) {
		//			//cout << id << ",," << i <<"," << j;
		//			fout << originalMatchingResultX[id][i][j] << " ";
		//		}
		//		fout << endl;
		//	}fout << endl;
		//	for (int i=0;i<height;i++) {
		//		for (int j=0;j<width;j++) {
		//			//cout << id << ",," << i <<"," << j;
		//			fout << originalMatchingResultY[id][i][j] << " ";
		//		}
		//		fout << endl;
		//	}fout << endl;
		//}
		fout.close();
		cout << "...end"<<endl;
		return true;
	}

	bool LoadFromFile(string filename)
	{
		//save the data into file
		std::fstream fin(filename.c_str(), std::ios::in);
		if (!fin) {
			cout << "Can not open " << filename << " to load data..." << endl;
			return false;
		}
		// first save parameters

		fin >> ViewID;
		char buff[256]; fin.getline(buff, 256);
		fin.getline(buff, 256); imageName = string(buff);
		fin.getline(buff, 256); maskName = string(buff);

		double _K[9],_R[9],_T[3];
		for (int j = 0; j < 9; ++j) {
			fin >> _K[j];
		}
		for (int j = 0; j < 9; ++j) {
			fin >> _R[j];
		}
		for (int j = 0; j < 3; ++j) {
			fin >> _T[j];
		}

		matrixCopy(_K,K,9);
		matrixCopy(_R,R,9);
		matrixCopy(_T,T,3);

		image = cvLoadImage(imageName.c_str(),0); m_Image = cvLoadImage(imageName.c_str(),1);
		height = image->height;
		width = image->width;
		printf("Reading image %s... Height %d x Width %d\n", imageName.c_str(), height, width);		
		calP();
		calM();

		calC();
		printf("center = %10.3lf %10.3lf %10.3lf\n",C[0],C[1],C[2]);

		calViewingVector();
		printf("viewing vector = %10.3lf %10.3lf %10.3lf\n",viewingVector[0],viewingVector[1],viewingVector[2]);

		calKinv();
		ucharImage = new uchar[height * width];
		for (int y=0;y<height;y++)
		{
			for (int x=0;x<width;x++)
			{
				ucharImage[y*width+x] = (uchar)cvGet2D(image,y,x).val[0];
				//if (ucharImage[y*width+x]==0) printf("0"); else printf("1");
			}
			//printf("\n");
		}

		if (scaleX.length()>0)
		{
			xImage = cvLoadImage((path+(string)"/"+scaleX+imageName.substr(imageName.find_last_of("\\")+1)).c_str(),0);
			ucharXImage = new uchar[xImage->height*xImage->width];
			for (int y=0;y<xImage->height;y++)
			{
				for (int x=0;x<xImage->width;x++)
				{
					ucharXImage[y*xImage->width+x] = (uchar)cvGet2D(xImage,y,x).val[0];
				}
			}
		} else 
		{
			xImage = image;
			ucharXImage = ucharImage;
		}

		readSilhouette(buff);

		// save the matching result of mvs algorithm
		//fin >> stereoPairs[0] >> stereoPairs[1];
		//fin >> verticalStereo[0] >> verticalStereo[1];
		//for (int id = 0; id < 2; id ++) {
		//	originalMatchingResultX[id] = allocate<int>(height,width);
		//	originalMatchingResultY[id] = allocate<int>(height,width);	
		//	for (int i=0;i<height;i++) {
		//		for (int j=0;j<width;j++) {
		//			fin >> originalMatchingResultX[id][i][j];
		//		}
		//	}
		//	for (int i=0;i<height;i++) {
		//		for (int j=0;j<width;j++) {
		//			fin >> originalMatchingResultY[id][i][j];
		//		}
		//	}
		//}
		fin.close();
		return true;
	}


};


#endif