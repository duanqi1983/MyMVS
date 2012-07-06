extern int nViews;
extern viewPoint** listViewPoints;
extern bool **paired;
extern double* camera_dis;
extern double *thetas;
extern double medianD;
extern float **** daisyMEM;

double _max_angle = -1.0;
double _max_dist = -1.0;
int max_theta = 180; int min_theta = 0; int theta_step = 5;
int total_angle_num = 0;

double LuminanceColor[3];
double Gamma = 200;
double LightPosition[3];
extern char buffer[255];

extern Engine *m_ep;
extern MyMesh ObjTriMesh;
vector<double> Vertex_Color_List;
vector< vector<double> > Vertex_PSNormal_list;
int FittingOrder = 12;
extern double SplineFittingDecreasePara;
//vector< list< pair<double, pair<int, double> > > > VertexIntensityRecords;//view_angle, view_index, intensity

void Mat_Print(string text, cv::Mat m)
{
	fprintf(stderr, "\n%s", text.c_str());
	fprintf(stderr, "\nrow = %d; col = %d\n",m.rows, m.cols);
	for (int i = 0; i < m.rows; i ++) {
		for (int j = 0; j < m.cols; j ++) {
			fprintf(stderr, "%lf  ", m.at<double>(i,j));
		}
		fprintf(stderr, "\n");
	}
}

template <class T>
bool SaveArray (T* a, int size, const char* filename) {
#ifdef _DEBUG
	fstream fout(filename, ios::out);
	if (!fout) {
		cout << "Can not open " << filename << " to save data..." << endl;
		return false;
	}
	for (int i = 0; i < size; i ++) {
		fout << a[i] << endl;
	}
	fout.close();
#endif
	return true;
}

void chooseStereoPairs()
{

	// initialize
	int *pairs = new int [nViews * 2];
	double *metric_theta_d = new double[nViews *2];
	for (int i=0;i<nViews;i++) pairs[i*2] = pairs[i*2+1] = -1;

	// find all distances
	camera_dis = new double[nViews*nViews];	
	printf("\nDistances between cameras: \n");
	for (int i=0;i<nViews;i++) 
	{
		for (int j=0;j<nViews;j++) 
		{
			camera_dis[i*nViews+j] = matrixDist(listViewPoints[i]->C,listViewPoints[j]->C,3);
			printf("%10.3lf",camera_dis[i*nViews+j]);
		}
		printf("\n");
	}

	// find all thetas
	thetas = new double[nViews*nViews];
	printf("\nAngles between cameras: \n");
	for (int i=0;i<nViews;i++)
	{
		for (int j=0;j<nViews;j++) 
		{
			if (i!=j) thetas[i*nViews+j] = calAngle(listViewPoints[i]->viewingVector,listViewPoints[j]->viewingVector);
			else thetas[i*nViews+j] = 0;
			printf("%10.3lf",thetas[i*nViews+j]);
		} 
		printf("\n");
	}

	// find median distance
	vector<double> listD;
	listD.clear();
	for (int i=0;i<nViews;i++)
		for (int j=i+1;j<nViews;j++) listD.push_back(camera_dis[i*nViews+j]);
	sort(listD.begin(),listD.end());
	medianD = listD[listD.size()/2];

	printf("\nMedian distance = %lf\n",medianD);

	// select pair
	for (int i=0;i<nViews;i++)
	{
		for (int j=0;j<nViews;j++)
		{
			if (i!=j)
				if (theta_lower<thetas[i*nViews+j] && thetas[i*nViews+j]<theta_upper)
					if (distance_lower*medianD < camera_dis[i*nViews+j] && camera_dis[i*nViews+j] < distance_upper*medianD)
					{
						double metricValue = thetas[i*nViews+j] * camera_dis[i*nViews+j];
						if (pairs[i*2]==-1)
						{
							pairs[i*2] = j;
							metric_theta_d[i*2] = metricValue;
						} else
							if (metricValue < metric_theta_d[i*2])
							{
								pairs[i*2+1] = pairs[i*2];
								metric_theta_d[i*2+1] = metric_theta_d[i*2];

								pairs[i*2] = j;
								metric_theta_d[i*2] = metricValue;
							} else
								if (pairs[i*2+1] ==-1 || metricValue < metric_theta_d[i*2+1])
								{
									pairs[i*2+1] = j;
									metric_theta_d[i*2+1] = metricValue;
								}
					}
		}
	}

	paired = new bool*[nViews];
	for (int i=0;i<nViews;i++) paired[i] = new bool[nViews];

	for (int i=0;i<nViews;i++) 
	{
		listViewPoints[i]->stereoPairs[0] = pairs[i*2];
		listViewPoints[i]->stereoPairs[1] = pairs[i*2+1];
		for(int j=0;j<nViews;j++) paired[i][j] = false;
		if (pairs[i*2]>=0) 
		{
			paired[i][pairs[i*2]] = true;
			paired[pairs[i*2]][i] = true;
		}
		if (pairs[i*2+1]>=0) 
		{
			paired[i][pairs[i*2+1]] = true;
			paired[pairs[i*2+1]][i] = true;
		}
	}

	printf("\nStereo pairs:\n");
	for (int i=0;i<nViews;i++) printf("%d : %d : %d\n",
		i,listViewPoints[i]->stereoPairs[0],listViewPoints[i]->stereoPairs[1]);

	deallocate(pairs);
	deallocate(metric_theta_d);
}


void load_daisy_image(CvMat *src, uchar* &im, int &h, int &w)
{
	h = src->height;
	w = src->width;
	im = new uchar[h * w];
	for (int y=0;y<h;y++)
		for (int x=0;x<w;x++)
		{
			im[y*w+x] = (uchar)cvGet2D(src,y,x).val[0];
			//printf("%d ", im[y*w+x]);
		}

}

void load_daisy_image(IplImage *src, uchar* &im, int &h, int &w)
{
	h = src->height;
	w = src->width;
	im = new uchar[h * w];
	for (int y=0;y<h;y++)
		for (int x=0;x<w;x++)
		{
			im[y*w+x] = (uchar)cvGet2D(src,y,x).val[0];
			//printf("%d ", im[y*w+x]);
		}

}


daisy* initDaisy(CvMat *src, double rad, int radq, int thq, int histq)
{


	uchar *im = NULL;
	int h,w;
	load_daisy_image(src,im,h,w);		

	daisy *desc = new daisy();
	desc->verbose(0);
	desc->set_image(im,h,w);
	desc->set_parameters(rad,radq,thq,histq);
	//desc->set_parameters(daisyParas[0],daisyParas[1],daisyParas[2],daisyParas[3]);
	desc->initialize_single_descriptor_mode();
	desc->compute_descriptors();
	desc->normalize_descriptors();
	return desc;
}


daisy* initDaisy(IplImage *src, double rad, int radq, int thq, int histq)
{

	uchar *im = NULL;
	int h,w;
	load_daisy_image(src,im,h,w);		

	daisy *desc = new daisy();
	desc->verbose(0);
	desc->set_image(im,h,w);
	desc->set_parameters(rad,radq,thq,histq);
	//desc->set_parameters(daisyParas[0],daisyParas[1],daisyParas[2],daisyParas[3]);
	desc->initialize_single_descriptor_mode();
	desc->compute_descriptors();
	desc->normalize_descriptors();
	return desc;
}



void stereoMatching()
{	
	for (int i=0;i<nViews;i++)
		for (int id=0; id<2; id++)
			if (listViewPoints[i]->stereoPairs[id]>=0)
			{				
				CvMat *mx1_map, *my1_map, *mx2_map, *my2_map;
				CvMat *image_original_rectified, *image_paired_rectified;

				printf("\nRectifying image %d and %d...\n", i,listViewPoints[i]->stereoPairs[id]);
				listViewPoints[i]->rectifyImages(id,listViewPoints[listViewPoints[i]->stereoPairs[id]], 
					mx1_map, my1_map, mx2_map, my2_map,image_original_rectified,image_paired_rectified);

				printf("Calculate DAISY features in rectifed Images: %d %d...\n",i,listViewPoints[i]->stereoPairs[id]);
				//double timer_start = (double)getTickCount();
				daisy *desc_original_rectified = initDaisy(image_original_rectified,daisyParas[0],daisyParas[1],daisyParas[2],daisyParas[3]);
				//float *thor;
				//for (int dx=0; dx<listViewPoints[i]->height; dx++)
				//	for (int dy=0; dy<listViewPoints[i]->width; dy++) desc_original_rectified->get_descriptor(dx,dy,thor);
				//printf("\nTime = %lfs\n",((double)getTickCount()-timer_start)/getTickFrequency());

				//timer_start = (double)getTickCount();
				//desc_original_rectified = initDaisyMode2(listViewPoints[i]->image,daisyParas[0],daisyParas[1],daisyParas[2],daisyParas[3]);
				//thor = new float[daisyLength];
				//for (int dx=0; dx<listViewPoints[i]->height; dx++)
				//	for (int dy=0; dy<listViewPoints[i]->width; dy++) desc_original_rectified->get_descriptor(dx,dy,0,thor);
				//printf("\nTime = %lfs\n",((double)getTickCount()-timer_start)/getTickFrequency());
				//getch();
				//exit(1);
				daisy *desc_paired_rectified = initDaisy(image_paired_rectified,daisyParas[0],daisyParas[1],daisyParas[2],daisyParas[3]);


				char buffer [100];
				string rectifedImageName;

				rectifedImageName = path + "_output/rectified" + string(itoa(i,buffer,10)) + "_" +string(itoa(listViewPoints[i]->stereoPairs[id],buffer,10))+"_original.png";
				cvSaveImage(rectifedImageName.c_str(),image_original_rectified);

				rectifedImageName = path + "_output/rectified" + string(itoa(i,buffer,10)) + "_" +string(itoa(listViewPoints[i]->stereoPairs[id],buffer,10))+"_paired.png";
				cvSaveImage(rectifedImageName.c_str(),image_paired_rectified);


				printf("Stereo matching %d : %d     -    ",i,listViewPoints[i]->stereoPairs[id]);
				printf(listViewPoints[i]->verticalStereo[id] ? "Vertically rectified\n" : "Horizontally rectified\n");

				int **stereoMatchingResult = allocate<int>(listViewPoints[i]->height,listViewPoints[i]->width);
				listViewPoints[i]->stereoMatching(id, desc_original_rectified, desc_paired_rectified, stereoMatchingResult);



				listViewPoints[i]->originalImageMatching(id, mx1_map, my1_map, mx2_map, my2_map,stereoMatchingResult,
					desc_original_rectified, desc_paired_rectified);

				// print out matches
				CvMat* tmpImage = cvCreateMat( listViewPoints[i]->height, listViewPoints[i]->width, CV_8U );
				CvScalar s0;
				CvScalar s1;
				s1.val[0] = 0;			
				s0.val[0] = 255;

				for (int u=0;u<listViewPoints[i]->height;u++)
					for (int v=0;v<listViewPoints[i]->width;v++)
						if (stereoMatchingResult[u][v]>=0) 
							cvSet2D(tmpImage,u,v,s1); else cvSet2D(tmpImage,u,v,s0);

				//char buffer[100];
				string tmpName;
				tmpName = path + "_output/foundMatch" + string(itoa(i,buffer,10))+"_"+string(itoa(listViewPoints[i]->stereoPairs[id],buffer,10))+".png";
				cvSaveImage(tmpName.c_str(),tmpImage);
				cvReleaseMat(&tmpImage);



				delete desc_original_rectified;
				delete desc_paired_rectified;
				cvReleaseMat(&image_original_rectified);
				cvReleaseMat(&image_paired_rectified);
				cvReleaseMat(&mx1_map);
				cvReleaseMat(&my1_map);
				cvReleaseMat(&mx2_map);
				cvReleaseMat(&my2_map);
				deallocate(stereoMatchingResult,listViewPoints[i]->height);

				//testMatchingResult(i,id);
			}

}

void buildTrack(track *trk, int viewId, int x, int y, bool *collectedView, int depth)
{	
	if (depth>=buildTrackDepth) return;
	for (int k=0;k<2;k++)
	{
		int nextView = listViewPoints[viewId]->stereoPairs[k];
		if (nextView>=0 )
			if (!collectedView[nextView])
			{	
				int xx = listViewPoints[viewId]->originalMatchingResultX[k][x][y];
				int yy = listViewPoints[viewId]->originalMatchingResultY[k][x][y];
				if (xx>=0)		
					if (listViewPoints[nextView]->sil[xx][yy])
						if (!listViewPoints[nextView]->collected[xx][yy])
						{
							trk->addView(nextView,xx,yy);
							collectedView[nextView] = true;	
							buildTrack(trk,nextView,xx,yy,collectedView, depth+1);						
						}
			}
	}
}

void buildTracks()
{
	printf("\nBuilding tracks...\n");

	// initialize flags
	for (int i=0;i<nViews;i++)
	{
		listViewPoints[i]->collected = allocate<bool>(listViewPoints[i]->height,listViewPoints[i]->width);
		for (int x=0;x<listViewPoints[i]->height;x++)
			for (int y=0;y<listViewPoints[i]->width;y++) listViewPoints[i]->collected[x][y] = false;
	}

	// build track	     
	int totalTracks = 0;
	for (int i=0;i<nViews;i++) 
	{
		int cntTracks = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static)
#endif		
		for (int x=0;x<listViewPoints[i]->height;x++)
		{
			bool *collectedView = new bool[nViews];
			for (int y=0;y<listViewPoints[i]->width;y++)
				if (listViewPoints[i]->sil[x][y])
					if (!listViewPoints[i]->collected[x][y])
					{					
						track trk;
						trk.addView(i,x,y);

						for (int l=0;l<nViews; l++) collectedView[l] = false;
						collectedView[i] = true;
						buildTrack(&trk, i, x, y, collectedView, 1);	

						bool reliable = true;

						if (trk.size>=beta) 
						{	
							while (trk.size>=beta)
							{
								// triangulate track
								triangulateTrack(&trk);
								// optimize by LMA
								optimizeByLMA(&trk);
								// remove mis-matched pixels
								if (trk.reliableViews < trk.size && trk.reliableViews>=beta ) trk.removeMismatchedPixels(); else break;
							}
							if (trk.reliableViews>=beta)
							{
								// check with silhouette
								for (int j=0;j<nViews;j++)
								{
									int px[2];
									matrixProject(trk.X,listViewPoints[j]->P,px,1);
									// filter by silhouette
									if (matrixBound(px[0],0,listViewPoints[j]->height) &&
										matrixBound(px[1],0,listViewPoints[j]->width))
									{
										if (!listViewPoints[j]->sil[px[0]][px[1]])
										{
											reliable = false;
											break;
										}
									} else if (isObject) 
									{
										reliable = false;
										break;
									}
								}
								// check by angle between viewing rays
								if (reliable)
								{
									for (int viewI=0; viewI<trk.size; viewI++)
									{
										for (int viewJ = viewI+1; viewJ<trk.size; viewJ++)
											if (paired[viewI][viewJ])
											{
												double v1[3],v2[3];
												v1[0] = listViewPoints[trk.views[viewI]]->C[0] - trk.X[0];
												v1[1] = listViewPoints[trk.views[viewI]]->C[1] - trk.X[1];
												v1[2] = listViewPoints[trk.views[viewI]]->C[2] - trk.X[2];
												matrixNorm(v1,3);

												v2[0] = listViewPoints[trk.views[viewJ]]->C[0] - trk.X[0];
												v2[1] = listViewPoints[trk.views[viewJ]]->C[1] - trk.X[1];
												v2[2] = listViewPoints[trk.views[viewJ]]->C[2] - trk.X[2];
												matrixNorm(v2,3);

												double inbetweenAngle = calAngle(v1,v2);
												if (inbetweenAngle<theta_lower || inbetweenAngle> theta_upper)
												{
													reliable = false;
													break;
												}
											}
											if (!reliable) break;
									}
								}

								if (reliable)
								{
									for (int j=0;j<trk.size;j++) 
										if (trk.reliable[j])
#pragma omp critical
										{
											listViewPoints[trk.views[j]]->collected[trk.px[j]][trk.py[j]] = true;		
										}
#pragma omp critical
										{
											tracks.push_back(trk); 
											cntTracks++;
										}

								}						
							} else reliable = false;

						}  else reliable = false;
					}
					deallocate(collectedView);
		}

		printf("View %d : collected %d\n", i, cntTracks);
		totalTracks += cntTracks;

	}

	// output
	printf("Total of built tracks: %d\n", totalTracks);

	int nPixels = 0;
	int maxLengthTrack = 0;
	for (int i=0; i<tracks.size(); i++)
	{
		nPixels += tracks[i].size;
		if (maxLengthTrack<tracks[i].size)
		{
			printf("\n\n");
			for (int j=0;j<tracks[i].size;j++) printf("%d %d %d\n", tracks[i].views[j], tracks[i].px[j], tracks[i].py[j]);
		}
		maxLengthTrack = std::max(maxLengthTrack, tracks[i].size);
	}
	printf("Total pixels = %d    -     max Length Track = %d\n",nPixels, maxLengthTrack);

	//// clear memory
	//for (int i=0; i<nViews; i++) 
	//{
	//	if (listViewPoints[i]->stereoPairs[0]>=0)
	//	{
	//		deallocate(listViewPoints[i]->originalMatchingResultX[0],listViewPoints[i]->height);
	//		deallocate(listViewPoints[i]->originalMatchingResultY[0],listViewPoints[i]->height);
	//	}
	//	if (listViewPoints[i]->stereoPairs[1]>=0)
	//	{		
	//		deallocate(listViewPoints[i]->originalMatchingResultX[1],listViewPoints[i]->height);		
	//		deallocate(listViewPoints[i]->originalMatchingResultY[1],listViewPoints[i]->height);
	//	}
	//}

}



daisy* initDaisyMode2(IplImage *src, double rad, int radq, int thq, int histq)
{
	uchar *im = NULL;
	int h,w;
	load_daisy_image(src,im,h,w);		

	daisy *desc = new daisy();
	desc->verbose(0);
	desc->set_image(im,h,w);
	desc->set_parameters(rad,radq,thq,histq);
	//desc->set_parameters(daisyParas[0],daisyParas[1],daisyParas[2],daisyParas[3]);
	//desc->set_normalization(NRM_FULL);
	desc->initialize_single_descriptor_mode();

	//desc->compute_descriptors();
	//desc->normalize_descriptors();

	return desc;
}

void calNormals()
{	

	// save track pixels to images


	
	for (int i=0;i<nViews;i++)
	{
		CvMat* tmpImage = cvCreateMat( listViewPoints[i]->height, listViewPoints[i]->width, CV_8U );
		CvScalar s0;
		CvScalar s1;
		s1.val[0] = 0;
		s1.val[1] = 0;
		s1.val[2] = 0;				
		s0.val[0] = 255;
		s0.val[1] = 255;
		s0.val[2] = 255;

		for (int u=0;u<listViewPoints[i]->height;u++)
			for (int v=0;v<listViewPoints[i]->width;v++)
				if (listViewPoints[i]->collected[u][v]) 
					cvSet2D(tmpImage,u,v,s1); else cvSet2D(tmpImage,u,v,s0);
		char buffer[100];
		string tmpName;
		tmpName = path + "TempResult/track" + string(itoa(i,buffer,10))+".png";
		cvSaveImage(tmpName.c_str(),tmpImage);
		cvReleaseMat(&tmpImage);
	}
	


	//
	cntF = 0;
	printf("\nPre-calculating Rij , Tij...");

	Rij = new double**[nViews];
	Tij = new double**[nViews];
	for (int i=0; i<nViews; i++)
	{
		Rij[i] = new double*[nViews];
		Tij[i] = new double*[nViews];
		for (int j=0;j<nViews;j++)
			if (i==j) 
			{
				Rij[i][j] = NULL;
				Tij[i][j] = NULL;
			}
			else
			{	
				Rij[i][j] = new double[9];
				Tij[i][j] = new double[3];

				cv::Mat Ri = cv::Mat(3,3,CV_64F,listViewPoints[i]->R);
				cv::Mat Rj = cv::Mat(3,3,CV_64F,listViewPoints[j]->R);
				cv::Mat _Rij = Rj * Ri.inv();
				for (int k=0;k<9;k++) Rij[i][j][k] = _Rij.at<double>(k/3,k%3);

				cv::Mat Ti = cv::Mat(3,1,CV_64F,listViewPoints[i]->T);
				cv::Mat Tj = cv::Mat(3,1,CV_64F,listViewPoints[j]->T);
				cv::Mat _Tij = Tj - _Rij * Ti;
				for (int k=0;k<3;k++) Tij[i][j][k] = _Tij.at<double>(k,0);				
			}
	}
	
	printf("\nOptimizing normals...\n");
	
	//printf("\nGet DAISY MODE 2 in view ");
	//daisyMEM = new float***[nViews];
	//int cnt = 0;
	for (int i=0;i<nViews;i++) 
	{
		//printf("%d ",i);
		
		/*
		listViewPoints[i]->desc_original = initDaisy(listViewPoints[i]->image,R_tangent,daisyParas[1],daisyParas[2],daisyParas[3]);
		daisyMEM[i] = new float **[listViewPoints[i]->height];
		for (int x=0;x<listViewPoints[i]->height;x++)
		{
			daisyMEM[i][x] = new float*[listViewPoints[i]->width];
			for (int y=0;y<listViewPoints[i]->width;y++)
			{
				daisyMEM[i][x][y] = NULL;
				if (listViewPoints[i]->collected[x][y])
				{
					cnt++;
					
					daisyMEM[i][x][y] = new float[daisyLength];
					float* thor = NULL;
					listViewPoints[i]->desc_original->get_descriptor(x,y,thor);
					for (int j=0;j<daisyLength;j++) daisyMEM[i][x][y][j] = thor[j];
					
				}
			}
		}
		delete listViewPoints[i]->desc_original;
		*/
		//cvSaveImage("0001.png",listViewPoints[i]->image);


		// apply Gaussian filter 5x5
		//cvSmooth(listViewPoints[i]->image, listViewPoints[i]->image, CV_GAUSSIAN, 51);

		//char buffer[100];
		//string tmpName;
		//tmpName = path + "_output/GaussianFilter" + string(itoa(i,buffer,10))+".png";
		//cvSaveImage(tmpName.c_str(),listViewPoints[i]->image);

		// init Daisy Mode 2
		//listViewPoints[i]->desc_original = initDaisyMode2(listViewPoints[i]->image,daisyParasNormals[0],daisyParasNormals[1],daisyParasNormals[2],daisyParasNormals[3]);

		/*
		float *thor0 = NULL;
		listViewPoints[i]->desc_original->get_descriptor(2,89,thor0);
		for (int t=0; t<daisyNormalsLength; t++)
		{
			printf("%10.9lf ", thor0[t]);
			if ((t+1)%daisyParasNormals[3]==0) printf("\n");
		}
		*/
	}
	//printf("\nDaisy cache : %d descriptors\n",cnt);

	/*
	// test

	float *thor0 = new float[daisyLength];
	float *thor1 = new float[daisyLength];

	listViewPoints[1]->desc_original->get_descriptor(2,79,0,thor0);
	for (int z=0;z<daisyLength;z++) 
	{
		printf("%f ", thor0[z]);
		if (z%daisyParas[3]==3) printf("\n");
	}
	printf("\n");	
	listViewPoints[2]->desc_original->get_descriptor(6.47,73.19,0,thor1);
	for (int z=0;z<daisyLength;z++) 
	{
		printf("%f ", thor1[z]);
		if (z%daisyParas[3]==3) printf("\n");
	}
	printf("\n");
	listViewPoints[0]->desc_original->get_descriptor(6.38,84.69,0,thor1);
	for (int z=0;z<daisyLength;z++) 
	{
		printf("%f ", thor1[z]);
		if (z%daisyParas[3]==3) printf("\n");
	}
	printf("\n");
	listViewPoints[1]->desc_original->get_descriptor(2,79,0,thor1);
	for (int z=0;z<daisyLength;z++) 
	{
		printf("%f ", thor1[z]);
		if (z%daisyParas[3]==3) printf("\n");
	}
	printf("\n");
	printf("cc = %lf\n\n\n", cross_correlation(thor0, thor1, daisyLength));
	printf("\n");


	*/

	//---------------------------------------------------------------------
	

	for (int i=0; i<nViews; i++) 
		deallocate(listViewPoints[i]->collected,listViewPoints[i]->height);

	
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic ,1)
#endif
	for (int i=0; i<tracks.size(); i++) 
		lbfgsFunction(i);

	printf("Number of calculating F function = %d\n", cntF);


	/*
	// delete DAISY memory
	for (int i=0;i<nViews;i++) 
		if (daisyMEM[i]!=NULL)
		{
			for (int x=0;x<listViewPoints[i]->height;x++) 
				if (daisyMEM[i][x]!=NULL)
				{
					for (int y=0;y<listViewPoints[i]->width;y++)
						if (daisyMEM[i][x][y]!=NULL) delete[] daisyMEM[i][x][y];
					delete[] daisyMEM[i][x];
				}
			delete listViewPoints[i]->desc_original;
			delete[] daisyMEM[i];
		}

		*/

	

	for (int i=0;i<nViews;i++)
	{
		for (int j=0;j<nViews;j++)
		{
			delete[] Rij[i][j];
			delete[] Tij[i][j];
		}
		delete[] Rij[i];
		delete[] Tij[i];

	}
	delete[] Rij;
	delete[] Tij;
	
}



void verifyTracks()
{
	// to find error
	printf("\n\nChecking tracks to find errors...\n");
	int cnt = 0;



	for (int j=0;j < tracks.size(); j++)
		if (tracks[j].valid)
	{
		bool ok = true;
		for (int k =0 ; k<tracks[j].size; k++)
		{
			for (int l = 0; l<tracks[j].size; l++)
			{
				if (k!=l)
					if (tracks[j].views[k] == tracks[j].views[l])
					{
						ok = false;
						break;
					}
			}	
			if (tracks[j].reliable[k])
			{
				double x[2];
				matrixProject(tracks[j].X, listViewPoints[tracks[j].views[k]]->P, x, 1);
				if (fabs(x[0] - tracks[j].px[k]) >gamma-deltaX || fabs(x[1] - tracks[j].py[k]) >gamma-deltaX )
				{
					ok = false;
					break;
				}
			}

		}
		if (ok)
		for (int k=0;k<nViews;k++)
		{
			int x[2];
			matrixProject(tracks[j].X,listViewPoints[k]->P,x,1);
			if (matrixBound(x[0],0,listViewPoints[k]->height) && matrixBound(x[1],0,listViewPoints[k]->width))
			{
				if (!listViewPoints[k]->sil[x[0]][x[1]])
				{
					ok =false;
					break;
				}
			} else
			if (isObject)
			{
				ok = false;
				break;
			}

		}

		double n[3];
		to3DNormal(tracks[j].X[3],tracks[j].X[4],&n[0],&n[1],&n[2]);

		if (ok)
		for (int k=0;k<tracks[j].size;k++)
		{
			int viewId = tracks[j].views[k];

			double XicamOi[3];
			matrixSubtract(listViewPoints[viewId]->C,tracks[j].X,XicamOi,3);
			matrixNorm(XicamOi,3);
			double inbetweenAngle = dotProduct(n,XicamOi,3);
			if (inbetweenAngle < angleThreshold + deltaX) 
			{
				ok = false;
				break;
			}
		}

		//if (ok) ok = crossCorrelationInHomography(j)>cross_correlation_threshold_homography;

		if (!ok)
		{
			tracks[j].valid = false;
			cnt++;
			/*
			printf("Error at track %d\n", j);
			for (int k = 0 ; k< tracks[j].size; k++)
			{
				printf("view %d     x = %d     y = %d\n", tracks[j].views[k], tracks[j].px[k], tracks[j].py[k]);
			}
			printf("\n");
			*/
		}
	} else cnt++;

	
	//vector<double> fxList;
	//for (int i=0;i<tracks.size();i++)
	//	if (tracks[i].valid) fxList.push_back(tracks[i].fx);
	//sort(fxList.begin(),fxList.end());
	//double thresholdFx = fxList[(int)(fxList.size()*.75)];
	//for (int i=0;i<tracks.size();i++)
	//	if (tracks[i].valid)
	//		if (tracks[i].fx>thresholdFx) 
	//		{
	//			tracks[i].valid = false;
	//			cnt++;
	//		}

	//		


	double minfx=1e9;
	double maxfx = -1e9;
	for (int j=0;j<tracks.size();j++)
		if (tracks[j].valid)
		{
			minfx = std::min(minfx,tracks[j].fx);
			maxfx = std::max(maxfx,tracks[j].fx);
		}
	printf("minfx = %lf\nmafx = %lf\n",minfx,maxfx);
	for (int j=0;j<tracks.size();j++)
		if (tracks[j].valid)
			if (tracks[j].fx>-normalNCC) 
			{
				tracks[j].valid = false;
				cnt++;
			}

	minfx=1e9;
	maxfx = -1e9;
	for (int j=0;j<tracks.size();j++)
		if (tracks[j].valid)
		{
			minfx = std::min(minfx,tracks[j].fx);
			maxfx = std::max(maxfx,tracks[j].fx);
		}
	printf("minfx = %lf\nmafx = %lf\n",minfx,maxfx);

	
	printf("Number of bad vertices : %d\n",cnt);
}

bool SaveMVSResult()
{
	//save the mvs result of Hoa's result
	string TempFolderName = MOptions.DirName + "TempResult\\";
	MyCreateDirectory(TempFolderName);
	// first save the paired result
	string paired_filename = TempFolderName + "pairedfile.txt";
	std::fstream pair_out(paired_filename.c_str(),std::ios::out);
	if (!pair_out) {
		cout << "Can not open " << paired_filename << " to save data..." << endl;
		return false;
	}
	pair_out << nViews << endl;
	for (int i = 0; i < nViews; i ++) {
		for (int j = 0; j < nViews; j ++) {
			pair_out << paired[i][j] << " ";
		}
		pair_out << endl;
	}
	pair_out.close();
	// save camera_dis/medianD result
	string camdis_filename = TempFolderName + "cameradistancefile.txt";
	std::fstream camdis_out(camdis_filename.c_str(),std::ios::out);
	if (!camdis_out) {
		std::cout << "Can not open " << camdis_filename << " to save data..." << std::endl;
		return false;
	}
	camdis_out << nViews << " " << medianD << endl;
	for (int i=0;i<nViews;i++) 
	{
		for (int j=0;j<nViews;j++) 
		{
			camdis_out << camera_dis[i*nViews+j] << " ";
		}
		camdis_out << endl;
	}
	camdis_out.close();
	// save thetas result
	string thetas_filename = TempFolderName + "thetasfile.txt";
	std::fstream thetas_out(thetas_filename.c_str(),std::ios::out);
	if (!thetas_out) {
		std::cout << "Can not open " << thetas_filename << " to save data..." << std::endl;
		return false;
	}
	thetas_out << nViews << std::endl;
	for (int i=0;i<nViews;i++) 
	{
		for (int j=0;j<nViews;j++) 
		{
			thetas_out << thetas[i*nViews+j] << " ";
		}
		thetas_out << std::endl;
	}
	thetas_out.close();

	// save CameraView data
	char buf[256];
	for (int i = 0; i < nViews; i ++) {
		sprintf(buf, "%sViewPoint%.5d.txt",TempFolderName.c_str(), i);
		if (!listViewPoints[i]->SaveToFile(string(buf))) {
			return false;
		}
	}

	// save the tracks result
	string tracks_filename = TempFolderName + "tracksfile.txt";
	std::fstream tracks_out(tracks_filename.c_str(), std::ios::out);
	if (!tracks_out) {
		std::cout << "Can not open " << tracks_filename << " to save data..." << std::endl;
		return false;
	}
	tracks_out << tracks.size() << std::endl; tracks_out << setprecision(8) << fixed;
	for (int idx = 0; idx < tracks.size(); idx ++) {
		for (int i = 0; i < 5; i ++) {
			tracks_out << tracks[idx].X[i] << " ";
		}tracks_out << std::endl;

		//to3DNormal(tracks[idx].X[3],tracks[idx].X[4],&tracks[idx].Norm[0],&tracks[idx].Norm[1],&tracks[idx].Norm[2]);
		//matrixNorm(tracks[idx].Norm, 3);
		//tracks_out << tracks[idx].Norm[0] << " " << tracks[idx].Norm[1] << " " << tracks[idx].Norm[2] << endl;

		//tracks_out << tracks[idx].size << std::endl;

		//for (int i = 0; i < tracks[idx].views.size(); i ++) {
		//	tracks_out << tracks[idx].views[i] << " ";
		//}tracks_out << std::endl;
		//for (int i = 0; i < tracks[idx].px.size(); i ++) {
		//	tracks_out << tracks[idx].px[i] << " ";
		//}tracks_out << std::endl;
		//for (int i = 0; i < tracks[idx].py.size(); i ++) {
		//	tracks_out << tracks[idx].py[i] << " ";
		//}tracks_out << std::endl;
		//for (int i = 0; i < tracks[idx].reliable.size(); i ++) {
		//	tracks_out << tracks[idx].reliable[i] << " ";
		//}tracks_out << std::endl;
		//tracks_out << tracks[idx].reliableViews << " " << tracks[idx].valid << " " 
		//	<< tracks[idx].fx << std::endl;
		//tracks_out << std::endl;
	}
	tracks_out.close();
	cout << "Save MVS Result Data Sucesssfully...." << std::endl;
	return true;
}

bool LoadMVSResult()
{
	std::cout<< "Start to load MVS result from file" << std::endl;
	//load the mvs result of Hoa's result
	string TempFolderName = MOptions.DirName + "TempResult\\";
	// first load the paired result
	string paired_filename = TempFolderName + "pairedfile.txt";
	std::fstream pair_in(paired_filename.c_str(), std::ios::in);
	if (!pair_in) {
		std::cout << "Can not open " << paired_filename << " to load data..." << std::endl;
		return false;
	}
	pair_in >> nViews;
	if (paired != NULL) {
		for (int i=0;i<nViews;i++) {
			delete[] paired[i];
		}
		delete[] paired;
	}
	paired = new bool*[nViews];
	for (int i=0;i<nViews;i++) 
		paired[i] = new bool[nViews];
	for (int i = 0; i < nViews; i ++) {
		for (int j = 0; j < nViews; j ++) {
			pair_in >> paired[i][j];
		}
	}
	pair_in.close();
	// load camera_dis/medianD result
	string camdis_filename = TempFolderName + "cameradistancefile.txt";
	std::fstream camdis_in(camdis_filename.c_str(), std::ios::in);
	if (!camdis_in) {
		std::cout << "Can not open " << camdis_filename << " to load data..." << std::endl;
		return false;
	}
	camdis_in >> nViews >> medianD;
	if (camera_dis != NULL) {
		delete[] camera_dis;
	}
	camera_dis = new double[nViews*nViews];	
	for (int i=0;i<nViews;i++) 
	{
		for (int j=0;j<nViews;j++) 
		{
			camdis_in >> camera_dis[i*nViews+j];
		}
	}
	camdis_in.close();
	// load thetas result
	string thetas_filename = TempFolderName + "thetasfile.txt";
	std::fstream thetas_in(thetas_filename.c_str(), std::ios::in);
	if (!thetas_in) {
		std::cout << "Can not open " << thetas_filename << " to load data..." << std::endl;
		return false;
	}
	thetas_in >> nViews;
	if (thetas != NULL) {
		delete[] thetas;
	}
	thetas = new double[nViews*nViews];	
	for (int i=0;i<nViews;i++) 
	{
		for (int j=0;j<nViews;j++) 
		{
			thetas_in >> thetas[i*nViews+j];
		}
	}
	thetas_in.close();

	// load CameraView data
	listViewPoints =  new viewPoint*[nViews];
	char buf[256];
	for (int i = 0; i < nViews; i ++) {
		listViewPoints[i] = new viewPoint();
		sprintf(buf, "%sViewPoint%.5d.txt",TempFolderName.c_str(), i);
		if (!listViewPoints[i]->LoadFromFile(string(buf))) {
			return false;
		}
	}

	// load the tracks result
	string tracks_filename = TempFolderName + "tracksfile.txt";
	std::fstream tracks_in(tracks_filename.c_str(), std::ios::in);
	if (!tracks_in) {
		std::cout << "Can not open " << tracks_filename << " to load data..." << std::endl;
		return false;
	}
	tracks.clear();
	int track_size;
	tracks_in >> track_size;
	std::cout << "Start to read tracks";
	for (int idx = 0; idx < track_size; idx ++) {
		if (idx % 5000 == 0) {
			cout << ".";
		}
		track t_track;
		for (int i = 0; i < 5; i ++) {
			tracks_in >> t_track.X[i];
		}
		//tracks_in >> t_track.Norm[0] >> t_track.Norm[1] >> t_track.Norm[2];
		//tracks_in >> t_track.size;
		//t_track.views.resize(t_track.size);
		//t_track.px.resize(t_track.size);
		//t_track.py.resize(t_track.size);
		//t_track.reliable.resize(t_track.size);

		//for (int i = 0; i < t_track.size; i ++) {
		//	tracks_in >> t_track.views[i];
		//}
		//for (int i = 0; i < t_track.size; i ++) {
		//	tracks_in >> t_track.px[i];
		//}
		//for (int i = 0; i < t_track.size; i ++) {
		//	tracks_in >> t_track.py[i];
		//}
		//for (int i = 0; i < t_track.size; i ++) {
		//	int t_ralia; tracks_in >> t_ralia;
		//	t_track.reliable[i] = (t_ralia>0) ? true : false;
		//}
		//tracks_in >> t_track.reliableViews >> t_track.valid >> t_track.fx;
		tracks.push_back(t_track);
	}
	tracks_in.close();
	std::cout << endl << "Load MVS Result Data Successfully...." << std::endl;
	return true;
}

bool SaveLightEstimationResult()
{
	// save the estimated light information of all the CameraView
	string TempFolderName = MOptions.DirName + "TempResult\\";
	string est_light_filename = TempFolderName + "LightEstimation.txt";
	std::fstream estlight_out(est_light_filename.c_str(), std::ios::out);
	if (!estlight_out) {
		std::cout << "Can not open " << est_light_filename << " to save data..." << std::endl;
		return false;
	}
	estlight_out << nViews << " " <<listViewPoints[0]->LightMatrix.rows << " " << listViewPoints[0]->LightMatrix.cols << std::endl;
	for (int idx = 0; idx < nViews; idx ++) {
		for (int i = 0; i < listViewPoints[idx]->LightMatrix.rows; i ++) {
			for (int j = 0; j < listViewPoints[idx]->LightMatrix.cols; j ++) {
				estlight_out << listViewPoints[idx]->LightMatrix.at<double>(i,j) << " ";
			}
			estlight_out << std::endl;
		}
		for (int i = 0; i < 3; i ++) {
			estlight_out << listViewPoints[idx]->LightDirection[i] << " ";
		}estlight_out << std::endl;
		for (int i = 0; i < 3; i ++) {
			estlight_out << listViewPoints[idx]->LightLuminance[i] << " ";
		}estlight_out << std::endl;
	}
	estlight_out.close();
	return true;
}

bool LoadLightEstimationResult()
{
	// load the estimated light information of all the CameraView
	string TempFolderName = MOptions.DirName + "TempResult\\";
	string est_light_filename = TempFolderName + "LightEstimation.txt";
	std::fstream estlight_in(est_light_filename.c_str(), std::ios::in);
	if (!estlight_in) {
		std::cout << "Can not open " << est_light_filename << " to load data..." << std::endl;
		return false;
	}
	int rows, cols;
	estlight_in >> nViews >> rows >> cols;
	for (int idx = 0; idx < nViews; idx ++) {
		listViewPoints[idx]->LightMatrix = cv::Mat(rows,cols,CV_64F);
		for (int i = 0; i < listViewPoints[idx]->LightMatrix.rows; i ++) {
			for (int j = 0; j < listViewPoints[idx]->LightMatrix.cols; j ++) {
				estlight_in >> listViewPoints[idx]->LightMatrix.at<double>(i,j) ;
			}
		}
		for (int i = 0; i < 3; i ++) {
			estlight_in >> listViewPoints[idx]->LightDirection[i];
		}
		for (int i = 0; i < 3; i ++) {
			estlight_in >> listViewPoints[idx]->LightLuminance[i];
		}
	}
	estlight_in.close();
	return true;
}

/*
void EstimateLightInformation(int iter_step)
{
	std::cout << "Start the light estimation process of step " << iter_step<<std::endl;
	//LuminanceColor(0) = LuminanceColor(1) = LuminanceColor(2) = 0.0;
	int total_count = 0;
	//estimate the light direction for each CameraView, using lambertian model and RANSAC process
	//repeat the estimation process many times randomly for ransac
	int Ransac_Iteration_Num = 300;
	cv::Mat Norm3 = cv::Mat(3,3,CV_64F);
#ifdef BRDFFITTING
	cv::Mat Clr = cv::Mat(3,4,CV_64F);
	cv::Mat LambdaI = cv::Mat(3,4,CV_64F); 
#else
	cv::Mat Clr = cv::Mat(3,1,CV_64F);
	cv::Mat LambdaI = cv::Mat(3,1,CV_64F); 
	vector<cv::Mat> vec_max_LambdaI; 
#endif
	int candidate_trackid[3];
	int current_trackid[3];
	int max_vote_num = -1;
	double vote_color_thres = 10;
	double ox[2]; int r,g,b, gray; // Color in the captured image

	for (int idx = 0; idx < nViews; idx ++) {
		boost::uniform_int<> distribution(0, listViewPoints[idx]->visibletrackid.size()-1) ;
		boost::mt19937 engine(time(NULL)) ;
		boost::variate_generator<boost::mt19937, boost::uniform_int<> > myrandom (engine, distribution);
		vec_max_LambdaI.clear();
		std::cout << "Visible material of image " << idx << " is: ";
		for (int i = 0; i < listViewPoints[idx]->visiblematerialcount.size(); i ++) {
			if (listViewPoints[idx]->visiblematerialcount[i] > listViewPoints[idx]->visibletrackid.size()*0.05) {
				std::cout << i << ",";
			}
		}
		std::cout << std::endl;
		for (int kn = 0; kn < K_Num; kn ++) { // for each kinds of material, we need to estimate its LambdaI separately
			cv::Mat MaxLambdaI;
			if (listViewPoints[idx]->visiblematerialcount[kn] < listViewPoints[idx]->visibletrackid.size()*0.05) {
				vec_max_LambdaI.push_back(cv::Mat(0,0,CV_64F));
				continue;
			}
			max_vote_num = -1;
			for (int num = 0; num < Ransac_Iteration_Num; ) {
				for (int j = 0; j < 3;) { // choose three point that are visible in this CameraView
					int temp = myrandom();
					if (tracks[listViewPoints[idx]->visibletrackid[temp]].material_label == kn) {
						current_trackid[j] = listViewPoints[idx]->visibletrackid[temp];
						//record normal information into matrix
						if (abs(matrixLength(tracks[current_trackid[j]].Norm, 3) - 1.0) > 0.00001) {
							matrixNorm(tracks[current_trackid[j]].Norm, 3);
						}
						Norm3.at<double>(0,j) = tracks[current_trackid[j]].Norm[0];
						Norm3.at<double>(1,j) = tracks[current_trackid[j]].Norm[1];
						Norm3.at<double>(2,j) = tracks[current_trackid[j]].Norm[2];

						//record color information
						matrixProject(tracks[current_trackid[j]].X,listViewPoints[idx]->P,ox,1);
						Clr.at<double>(j,0) = CV_IMAGE_ELEM(listViewPoints[idx]->image, uchar, (int)(ox[0]+0.5),(int)(ox[1]+0.5));
#ifdef BRDFFITTING
						Clr.at<double>(j,1) = CV_IMAGE_ELEM(listViewPoints[idx]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+2); //r
						Clr.at<double>(j,2) = CV_IMAGE_ELEM(listViewPoints[idx]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+1); //g
						Clr.at<double>(j,3) = CV_IMAGE_ELEM(listViewPoints[idx]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+0); //b
#endif
						j ++;
					}
					
				}
				//Mat_Print("Norm3", Norm3);Mat_Print("Color", Clr);
				LambdaI = Norm3.inv() * Clr; // the result of light direction and lambda parameter.
				//Mat_Print("LambdaI", LambdaI); 
				bool valid_flag = true;
				for (int jj = 0; jj < LambdaI.cols; jj ++) {
					double temp[3];
					temp[0] = LambdaI.at<double>(0,jj);temp[1] = LambdaI.at<double>(1,jj);temp[2] = LambdaI.at<double>(2,jj);
					if (matrixLength(temp,3) > 255.0) {
						valid_flag = false;
					}
				}
				if (!valid_flag) {
					continue;
				}
				num ++;
				int current_vote = 0; 
				for (int i = 0; i < listViewPoints[idx]->visibletrackid.size(); i ++) {
					if (tracks[listViewPoints[idx]->visibletrackid[i]].material_label == kn) {
						// vote for the calculate LambdaI if it satisfy the Lambertain model.
						matrixProject(tracks[listViewPoints[idx]->visibletrackid[i]].X,listViewPoints[idx]->P,ox,1);
						b = CV_IMAGE_ELEM(listViewPoints[idx]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+0);
						g = CV_IMAGE_ELEM(listViewPoints[idx]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+1);
						r = CV_IMAGE_ELEM(listViewPoints[idx]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+2);

						gray = CV_IMAGE_ELEM(listViewPoints[idx]->image, uchar, (int)(ox[0]+0.5),(int)(ox[1]+0.5));

						if (abs(matrixLength(tracks[listViewPoints[idx]->visibletrackid[i]].Norm, 3) - 1.0) > 0.00001) {
							matrixNorm(tracks[listViewPoints[idx]->visibletrackid[i]].Norm, 3);
						}
						cv::Mat Norm1 = cv::Mat(3,1,CV_64F,tracks[listViewPoints[idx]->visibletrackid[i]].Norm);
						cv::Mat Lam_Clr = LambdaI.t() * Norm1; // Color from Lambertain model

						if (abs(Lam_Clr.at<double>(0,0) - gray) < vote_color_thres 
#ifdef BRDFFITTING
							&&abs(Lam_Clr.at<double>(1,0) - r) < vote_color_thres 
							&&abs(Lam_Clr.at<double>(2,0) - g) < vote_color_thres 
							&&abs(Lam_Clr.at<double>(3,0) - b) < vote_color_thres 
#endif
							) {
								//cout << Lam_Clr.at<double>(0,0) << " " << r << " " << Lam_Clr.at<double>(1,0) << " " << g << " "
								//	 << Lam_Clr.at<double>(2,0) << " " << b << " " << Lam_Clr.at<double>(3,0) << " " << gray << endl;
								current_vote ++;
						}
					}	
				}
				if (current_vote > max_vote_num) {
					max_vote_num = current_vote;
					matrixCopy(current_trackid, candidate_trackid,3);
					LambdaI.copyTo(MaxLambdaI);
				}
			}
			Mat_Print("Max LambdaI",MaxLambdaI);
			vec_max_LambdaI.push_back(MaxLambdaI);
		}
		// average the estimated LambdaI result of different type of material to refine light direction
		vector<double> vec_Lambda; vec_Lambda.clear(); double temp_LD[3], Aver_LD[3]; Aver_LD[0] = Aver_LD[1] = Aver_LD[2] = 0.0;
		int valid_count = 0;
		for (int ii = 0; ii < vec_max_LambdaI.size(); ii ++) {
			if (vec_max_LambdaI[ii].rows < 1) {
				vec_Lambda.push_back(0);
				continue;
			}
			valid_count ++ ;
			temp_LD[0] = vec_max_LambdaI[ii].at<double>(0,0); 
			temp_LD[1] = vec_max_LambdaI[ii].at<double>(1,0); 
			temp_LD[2] = vec_max_LambdaI[ii].at<double>(2,0);
			std::cout<<std::endl<<temp_LD[0]<<", "<<temp_LD[1]<<", "<<temp_LD[2];
			vec_Lambda.push_back(matrixLength(temp_LD, 3));		
			matrixNorm(temp_LD,3);
			Aver_LD[0] = Aver_LD[0] + temp_LD[0]; Aver_LD[1] = Aver_LD[1] + temp_LD[1]; Aver_LD[2] = Aver_LD[2] + temp_LD[2];
		}
		matrixDiv(Aver_LD, 3, valid_count);			matrixNorm(Aver_LD, 3);	
		// end of ransac process, now record the id result and the light direction result
		listViewPoints[idx]->LightMatrix = cv::Mat(3,1,CV_64F);
		listViewPoints[idx]->LightMatrix.at<double>(0,0) = 255*Aver_LD[0];
		listViewPoints[idx]->LightMatrix.at<double>(1,0) = 255*Aver_LD[1];
		listViewPoints[idx]->LightMatrix.at<double>(2,0) = 255*Aver_LD[2];
		//Mat_Print("Final LambdaI", MaxLambdaI); 
		matrixCopy(Aver_LD, listViewPoints[idx]->LightDirection, 3);

#ifdef BRDFFITTING
		double lightL[3];
		for (int jj = 0; jj < 3; jj ++) {
			lightL[0] = = listViewPoints[idx]->LightMatrix.at<double>(0,jj);
			lightL[1] = = listViewPoints[idx]->LightMatrix.at<double>(1,jj);
			lightL[2] = = listViewPoints[idx]->LightMatrix.at<double>(2,jj);
			listViewPoints[idx]->LightLuminance[jj] = matrixLength(lightL,3); 
		}
#else
		listViewPoints[idx]->LightLuminance[0] = 
			listViewPoints[idx]->LightLuminance[1] = listViewPoints[idx]->LightLuminance[2] = 255;
#endif
		for (int i = 0; i < listViewPoints[idx]->visibletrackid.size(); i ++) {
			matrixProject(tracks[listViewPoints[idx]->visibletrackid[i]].X,listViewPoints[idx]->P,ox,1);
			double t_i = CV_IMAGE_ELEM(listViewPoints[idx]->image, uchar, (int)(ox[0]+0.5),(int)(ox[1]+0.5));
			double t_cos = dotProduct(Aver_LD, tracks[listViewPoints[idx]->visibletrackid[i]].Norm, 3);
			double t_lambda = t_i/(listViewPoints[idx]->LightLuminance[0] * t_cos);
			tracks[listViewPoints[idx]->visibletrackid[i]].AddViewLambda(idx,t_lambda);
			//tracks[listViewPoints[idx]->visibletrackid[i]].AddViewLambda(idx,vec_Lambda[tracks[listViewPoints[idx]->visibletrackid[i]].material_label]/listViewPoints[idx]->LightLuminance[0]);
		}
	}
	std::cout<<std::endl<< "...done..."<<std::endl;
}
*/

bool OptimizePointLightInformation(int iter_step, bool Using_Jacobian = false)
{
	//Still need to add \omega weight into it.
	//Optimize the single point light position $l$ and lambertain parameters $\gamma$
	mxArray* PointData = mxCreateDoubleMatrix(ObjTriMesh.n_vertices(), 8, mxREAL);
	mxArray* x_res = NULL; mxArray* f_val = NULL;
	size_t m_row = mxGetM(PointData);  size_t m_col = mxGetN(PointData);
	double* p_data = mxGetPr(PointData);   char eqnbuffer[500]; int idx = 0;
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++v_it, idx ++) {
		p_data[0*m_row + idx] = ObjTriMesh.point(v_it)[0];		
		p_data[1*m_row + idx] = ObjTriMesh.point(v_it)[1];		
		p_data[2*m_row + idx] = ObjTriMesh.point(v_it)[2];		
		p_data[3*m_row + idx] = ObjTriMesh.normal(v_it)[0];		
		p_data[4*m_row + idx] = ObjTriMesh.normal(v_it)[1];		
		p_data[5*m_row + idx] = ObjTriMesh.normal(v_it)[2];		
		p_data[6*m_row + idx] = Vertex_Color_List[idx];	
		p_data[7*m_row + idx] = 1.0; // the vertex weight, currently 1 for every point	
		//cout << ObjTriMesh.point(v_it)[0] << " " << ObjTriMesh.point(v_it)[1]  << " " << ObjTriMesh.point(v_it)[2] << " "
		//	<< ObjTriMesh.normal(v_it)[0] << " " << ObjTriMesh.normal(v_it)[1] << " " << ObjTriMesh.normal(v_it)[2] << " "
		//	<< Vertex_Color_List[idx]  << endl;
	}
	//for (int idx = 0; idx < tracks.size(); idx ++) {
	//	p_data[0*m_row + idx] = tracks[idx].X[0];
	//	p_data[1*m_row + idx] = tracks[idx].X[1];
	//	p_data[2*m_row + idx] = tracks[idx].X[2];
	//	p_data[3*m_row + idx] = tracks[idx].Norm[0];
	//	p_data[4*m_row + idx] = tracks[idx].Norm[1];
	//	p_data[5*m_row + idx] = tracks[idx].Norm[2];
	//	p_data[6*m_row + idx] = tracks[idx].avg_view_intensity;
	//	//cout << tracks[idx].X[0] << " " << tracks[idx].X[1] << " " << tracks[idx].X[2] << " "
	//	//	<< tracks[idx].Norm[0] << " " << tracks[idx].Norm[1] << " " << tracks[idx].Norm[2] << " "
	//	//	<< tracks[idx].avg_view_intensity << endl;
	//}
	engPutVariable(m_ep, "PD", PointData);
	engEvalString(m_ep, "save('PointData.mat','PD');");
	//finish save point information into .mat file;

	//Generate the ElFun file for the optimization process
	std::fstream m_out("c:\\MATLAB\\ElFun.m",std::ios::out);
	if (!m_out) {
		cout << "Can not open c:\\MATLAB\\ElFun.m to save data..." << endl;
		return false;
	}
	if (Using_Jacobian) {
		m_out << "function [F,J] = ElFun(x)" << endl;
	} else {
		m_out << "function F = ElFun(x)" << endl;
	}
	m_out << "load('PointData.mat');" << endl;
	//m_out << "[m,n] = size(PD);" << endl;
	m_out << "k=1:size(PD,1);" << endl;
	m_out << "F = sqrt(PD(k,8)).*(PD(k,7) - x(4).*( (PD(k,4)).*minus(x(1),PD(k,1)) + (PD(k,5)).*minus(x(2),PD(k,2)) + (PD(k,6)).*minus(x(3), PD(k,3)) )./sqrt( minus(x(1),PD(k,1)).^2 + minus(x(2),PD(k,2)).^2 + minus(x(3),PD(k,3)).^2 ));" << endl;
	if (Using_Jacobian) {
		m_out << "if nargout > 1" << endl;
		m_out << "  for k = 1:size(PD,1)" << endl;
		m_out << "   J(k,:) = [ -PD(k,8)^(1/2)*((PD(k,4)*x(4))/((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(1/2) - (x(4)*(2*x(1) - 2*PD(k,1))*(PD(k,4)*(x(1) - PD(k,1)) + PD(k,5)*(x(2) - PD(k,2)) + PD(k,6)*(x(3) - PD(k,3))))/(2*((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(3/2))),..." << endl;
		m_out << "     -PD(k,8)^(1/2)*((PD(k,5)*x(4))/((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(1/2) - (x(4)*(2*x(2) - 2*PD(k,2))*(PD(k,4)*(x(1) - PD(k,1)) + PD(k,5)*(x(2) - PD(k,2)) + PD(k,6)*(x(3) - PD(k,3))))/(2*((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(3/2))),... " << endl;
		m_out << "     -PD(k,8)^(1/2)*((PD(k,6)*x(4))/((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(1/2) - (x(4)*(2*x(3) - 2*PD(k,3))*(PD(k,4)*(x(1) - PD(k,1)) + PD(k,5)*(x(2) - PD(k,2)) + PD(k,6)*(x(3) - PD(k,3))))/(2*((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(3/2))),... " << endl;
		m_out << "     -(PD(k,8)^(1/2)*(PD(k,4)*(x(1) - PD(k,1)) + PD(k,5)*(x(2) - PD(k,2)) + PD(k,6)*(x(3) - PD(k,3))))/((x(1) - PD(k,1))^2 + (x(2) - PD(k,2))^2 + (x(3) - PD(k,3))^2)^(1/2)];" << endl;
		m_out << "  end" << endl;
		m_out << "end" << endl;
	}
	m_out.close();

	//double avg_lx, avg_ly, avg_lz, avg_gamma;
	//avg_lx = avg_ly = avg_lz = avg_gamma = 0.0; int count = 0;
	//boost::uniform_real<> distribution(0, 100) ;
	//boost::mt19937 engine(time(NULL)) ;
	//boost::variate_generator<boost::mt19937, boost::uniform_real<> > myrandom (engine, distribution);
	//for (int ii = 0; ii < 100; ii ++) {
	//	sprintf(eqnbuffer,"x0=[%f, %f, %f, %f];", myrandom(), myrandom(), myrandom(), 
	//		Gamma+myrandom()*20/100-10);//using previous gamma as initial guess
	//	//cout << eqnbuffer << endl;
	//	engEvalString(m_ep, string(eqnbuffer).c_str()); cout << eqnbuffer << " ";
	//	//engEvalString(m_ep, "options=optimset('Display','off');");
	//	if (Using_Jacobian) {
	//		engEvalString(m_ep, "options = optimset('Jacobian','on');");
	//		engEvalString(m_ep, "xres = lsqnonlin(@ElFun,x0,[-500,-500,-500,0],[500,500,500,255], options); save('xres.mat','xres');");
	//	} else {
	//		engEvalString(m_ep, "xres = lsqnonlin(@ElFun,x0,[-500,-500,-500,0],[500,500,500,255]); save('xres.mat','xres');");
	//	}
	//	while(!FileExisted("c:\\MATLAB\\xres.mat")) 
	//	{
	//		::Sleep(10000); cout << "Wait 1000 ms for the lsqnonlin process" << endl;
	//	}
	//	engEvalString(m_ep, "load('xres.mat');");
	//	x_res = engGetVariable(m_ep, "xres");
	//	double *x_data = (double*)mxGetData(x_res);
	//	LightPosition[0] = x_data[0]; LightPosition[1] = x_data[1]; LightPosition[2] = x_data[2]; Gamma = x_data[3]; 
	//	cout << LightPosition[0] << ", " << LightPosition[1] << ", " << LightPosition[2] << ", " << Gamma << endl;
	//	avg_lx += LightPosition[0];  avg_ly += LightPosition[1];   avg_lz += LightPosition[2];   avg_gamma += Gamma; count ++;
	//	engEvalString(m_ep, "delete('xres.mat');");
	//}
	//avg_lx = avg_lx/count; avg_ly = avg_ly/count; avg_lz = avg_lz/count; avg_gamma = avg_gamma/count;

	sprintf(eqnbuffer,"x0=[%f, %f, %f, %f];", LightPosition[0], LightPosition[1], LightPosition[2], Gamma);//using previous gamma as initial guess
	cout << "The initial guess " << eqnbuffer << endl;
	engEvalString(m_ep, "delete('xres.mat');");
	engEvalString(m_ep, eqnbuffer);
	//engEvalString(m_ep, "options=optimset('Display','off');");
	if (Using_Jacobian) {
		engEvalString(m_ep, "options = optimset('Jacobian','on');");
		engEvalString(m_ep, "xres = lsqnonlin(@ElFun,x0,[-500,-500,-500,0],[500,500,500,255], options);");
	} else {
		engEvalString(m_ep, "xres = lsqnonlin(@ElFun,x0,[-500,-500,-500,0],[500,500,500,255]);");
	}
	engEvalString(m_ep, "save('xres.mat','xres');");
	while(!FileExisted("c:\\MATLAB\\xres.mat")) 
	{
		::Sleep(10000); cout << "Wait 1000 ms for the lsqnonlin process" << endl;
	}
	engEvalString(m_ep, "load('xres.mat');");
	x_res = engGetVariable(m_ep, "xres");
	double *x_data = (double*)mxGetData(x_res);
	LightPosition[0] = x_data[0]; LightPosition[1] = x_data[1]; LightPosition[2] = x_data[2]; Gamma = x_data[3]; 
	cout <<"The final result is:" << LightPosition[0] << ", " << LightPosition[1] << ", " << LightPosition[2] << ", " << Gamma << endl;
	engEvalString(m_ep, "delete('xres.mat');");

	////Here use the RANSAC ideas as the above function
	//std::cout << "Start the light estimation process of step " << iter_step<<std::endl;
	//int total_count = 0;
	//int Ransac_Iteration_Num = 300;
	//int candidate_trackid[4];
	//int current_trackid[4];
	//int max_vote_num = -1;
	//double vote_color_thres = 5;

	//int sample_number = 4;

	//double nx[100],ny[100],nz[100],px[100],py[100],pz[100],c[100];
	//double lambdaI, lx, ly, lz;

	//for (int num = 0; num < Ransac_Iteration_Num; num ++) {
	//	for (int j = 0; j < sample_number; ) {
	//		int temp = myrandom(); 
	//		if (tracks[temp].avg_view_intensity < 200) {
	//			continue;
	//		}
	//		current_trackid[j] = temp;
	//		//std::cout << "Point " << temp << ", Normal is: " << tracks[temp].Norm[0] << ", "<<tracks[temp].Norm[1] << ", " << tracks[temp].Norm[2] <<
	//		//	", Position is: " <<tracks[temp].X[0] << ", "<<tracks[temp].X[1] << ", " << tracks[temp].X[2] << std::endl;
	//		nx[j] = tracks[temp].Norm[0];	ny[j] = tracks[temp].Norm[1];	nz[j] = tracks[temp].Norm[2];
	//		px[j] = tracks[temp].X[0];		py[j] = tracks[temp].X[1];		pz[j] = tracks[temp].X[2];	
	//		c[j]  = tracks[temp].avg_view_intensity;
	//		j++;
	//	}
	//	//calculate the light position $l$ and $\gamma$ based on these four points using Matlab solve function
	//	//first store the equation into .m file for fsolve 
	//	std::fstream m_out("c:\\MATLAB\\myfun.m",std::ios::out);
	//	if (!m_out) {
	//		cout << "Can not open c:\\MATLAB\\myfun.m to save data..." << endl;
	//		return false;
	//	}
	//	m_out << "function F = myfun(x)" << endl;
	//	m_out << "F = [" << endl;
	//	
	//	for (int j = 0; j < sample_number; j ++) {
	//		sprintf(eqnbuffer,"%lf - x(1)*( (%lf)*(x(2)-(%lf)) + (%lf)*(x(3)-(%lf)) + (%lf)*(x(4)-(%lf)) )/sqrt( (%lf-x(2))^2 + (%lf-x(3))^2 + (%lf-x(4))^2 )", 
	//			c[j],nx[j],px[j],ny[j],py[j],nz[j],pz[j],px[j],py[j],pz[j]);
	//		m_out << eqnbuffer;
	//		if (j != sample_number - 1) {
	//			m_out << ";" << endl;
	//		}
	//	}
	//	m_out << endl << "];"<<endl;
	//	/*
	//	sprintf(eqnbuffer,"%lf-x(1)*((%lf)*(x(2)-(%lf))+(%lf)*(x(3)-(%lf))+(%lf)*(x(4)-(%lf)))/sqrt((%lf-x(2))^2+(%lf-x(3))^2+(%lf-x(4))^2);", 
	//		c[0],nx[0],px[0],ny[0],py[0],nz[0],pz[0],px[0],py[0],pz[0]);
	//	//sprintf(eqnbuffer,"%lf+x(1)*[(%lf)+(%lf)+(%lf)] - x(1)*[(%lf)*x(2)+(%lf)*x(3)+(%lf)*x(4)];",
	//	//	c[0],nx[0]*px[0],ny[0]*py[0],nz[0]*pz[0],nx[0],ny[0],nz[0]);
	//	m_out << eqnbuffer << endl;
	//	sprintf(eqnbuffer,"%lf-x(1)*((%lf)*(x(2)-(%lf))+(%lf)*(x(3)-(%lf))+(%lf)*(x(4)-(%lf)))/sqrt((%lf-x(2))^2+(%lf-x(3))^2+(%lf-x(4))^2);", 
	//		c[1],nx[1],px[1],ny[1],py[1],nz[1],pz[1],px[1],py[1],pz[1]);
	//	//sprintf(eqnbuffer,"%lf+x(1)*[(%lf)+(%lf)+(%lf)] - x(1)*[(%lf)*x(2)+(%lf)*x(3)+(%lf)*x(4)];",
	//	//	c[1],nx[1]*px[1],ny[1]*py[1],nz[1]*pz[1],nx[1],ny[1],nz[1]);
	//	m_out << eqnbuffer << endl;
	//	sprintf(eqnbuffer,"%lf-x(1)*((%lf)*(x(2)-(%lf))+(%lf)*(x(3)-(%lf))+(%lf)*(x(4)-(%lf)))/sqrt((%lf-x(2))^2+(%lf-x(3))^2+(%lf-x(4))^2);", 
	//		c[2],nx[2],px[2],ny[2],py[2],nz[2],pz[2],px[2],py[2],pz[2]);
	//	//sprintf(eqnbuffer,"%lf+x(1)*[(%lf)+(%lf)+(%lf)] - x(1)*[(%lf)*x(2)+(%lf)*x(3)+(%lf)*x(4)];",
	//	//	c[2],nx[2]*px[2],ny[2]*py[2],nz[2]*pz[2],nx[2],ny[2],nz[2]);
	//	m_out << eqnbuffer << endl;
	//	sprintf(eqnbuffer,"%lf-x(1)*((%lf)*(x(2)-(%lf))+(%lf)*(x(3)-(%lf))+(%lf)*(x(4)-(%lf)))/sqrt((%lf-x(2))^2+(%lf-x(3))^2+(%lf-x(4))^2)", 
	//		c[3],nx[3],px[3],ny[3],py[3],nz[3],pz[3],px[3],py[3],pz[3]);
	//	//sprintf(eqnbuffer,"%lf+x(1)*[(%lf)+(%lf)+(%lf)] - x(1)*[(%lf)*x(2)+(%lf)*x(3)+(%lf)*x(4)]",
	//	//	c[3],nx[3]*px[3],ny[3]*py[3],nz[3]*pz[3],nx[3],ny[3],nz[3]);
	//	m_out << eqnbuffer << endl; m_out << "];"<<endl;
	//	*/
	//	m_out.close();

	//	engEvalString(m_ep, "cd c:\\MATLAB");
	//	sprintf(eqnbuffer,"x0=[%lf, %lf, %lf, %lf];", 200, (double)myrandom()/tracks.size(), (double)myrandom()/tracks.size(), (double)myrandom()/tracks.size());
	//	engEvalString(m_ep, string(eqnbuffer).c_str());
	//	//engEvalString(m_ep, "options=optimset('levenberg-marquardt', 0.05);");
	//	engEvalString(m_ep, "[x_res,fval] = fsolve(@myfun,x0);");
	//	//engEvalString(m_ep, "save('x_res.mat', 'x_res');");
	//	//engEvalString(m_ep, "load('x_res.mat');");
	//	x_res = engGetVariable(m_ep, "x_res");
	//	f_val = engGetVariable(m_ep, "fval");
	//	double *x_data = (double*)mxGetData(x_res);
	//	double *fval_data = (double*)mxGetData(f_val);
	//	Gamma = x_data[1]; LightPosition[0] = x_data[2]; LightPosition[1] = x_data[3]; LightPosition[2] = x_data[4];
	//	cout << x_data[1] << ", " << x_data[2] << ", " << x_data[3] << ", " << x_data[4] << "," << mxGetM(f_val) << ", "<<mxGetN(f_val) << ", " << fval_data[1] << endl;
	//	
	//}
	return true;
}

void InitializeProcess()
{
	//first record the id of track in each listviewpoint
	for (int i = 0; i < nViews; i ++) {
		listViewPoints[i]->visibletrackid.clear();
	}
	double ox[2];
	for (int idx = 0; idx < tracks.size(); idx ++) {
		matrixNorm(tracks[idx].Norm, 3);
		tracks[idx].view_intensity.resize(tracks[idx].views.size());
		double avg_intensity = 0.0;
		for (int i = 0; i < tracks[idx].views.size() ; i++) {
			listViewPoints[tracks[idx].views[i]]->visibletrackid.push_back(idx);
			//record the view intensity of certain track in each view
			matrixProject(tracks[idx].X,listViewPoints[tracks[idx].views[i]]->P,ox,1);
			tracks[idx].view_intensity[i] = CV_IMAGE_ELEM(listViewPoints[tracks[idx].views[i]]->image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5));
			avg_intensity += tracks[idx].view_intensity[i];
		}
		tracks[idx].avg_view_intensity = avg_intensity/tracks[idx].view_intensity.size();
		tracks[idx].SurfaceLambda = 1.0;
	}
}

void FilterModelProcess()
{
	cout << "Start the initial model filtering process...";
	double ox[2];
	//Filter the initial poisson model, delete some vertex and face not in object mask
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		int visible_count = 0;
		//project point back to each view point to see its validation
		for (int i = 0; i < nViews; i ++) {
			//double output_angle = calAngle(ObjTriMesh.normal(v_it).data(), listViewPoints[i]->InverseViewVector);
			//if (output_angle < 100) 
			{
				matrixProject(ObjTriMesh.point(v_it).data(), listViewPoints[i]->P,ox,1);
				if (matrixBound((int)(ox[0]+0.5), 0, listViewPoints[i]->height) && matrixBound((int)(ox[1]+0.5), 0, listViewPoints[i]->width)) {
					if(listViewPoints[i]->sil[(int)(ox[0]+0.5)][(int)(ox[1]+0.5)] == true) {
						visible_count ++;
					}
				}
			}
		}
		if (visible_count < 3) {
			ObjTriMesh.delete_vertex(v_it);
		}
	}
	ObjTriMesh.garbage_collection();
	cout <<"Done..."<<endl;
}

bool MyDataSortPredicate(const pair< double, pair<int, double> >& lhs, const pair< double, pair<int, double> >& rhs) 
{ 
	return lhs.first < rhs.first; 
}

void CalculateVertexIntensity(int step = 0)
{
	cout << "Start to calculate the vertex intensity value...";
	sprintf(buffer, "_%.2d", step);
	string resultfile = (MOptions.DirName+"VertexColor"+string(buffer)+".txt");
	if (FileExisted(resultfile.c_str())) {
		return;
	}
	ObjTriMesh.request_face_normals();	ObjTriMesh.update_face_normals();
	ObjTriMesh.request_vertex_normals();ObjTriMesh.update_vertex_normals();

	double ox[2]; double temp_color, tcolor; int temp_count, tcount; double nx, ny, nz, z1, z2; double cur_view_direction[3];
	//construct the initial normal/light matrix, to learn their relationship
	Vertex_Color_List.clear(); Vertex_Color_List.resize(ObjTriMesh.n_vertices());  std::list< pair<double, pair<int, double> > > tVectexIntensityList;
	mxArray* PointNormalIntensity = mxCreateDoubleMatrix(ObjTriMesh.n_vertices(), 6, mxREAL);
	size_t m_row = mxGetM(PointNormalIntensity);  size_t m_col = mxGetN(PointNormalIntensity);	double* p_data = mxGetPr(PointNormalIntensity);   
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		tVectexIntensityList.clear();
		int vertex_idx = v_it.handle().idx();
		for (int i = 0; i < nViews; i ++) {
			cur_view_direction[0] = listViewPoints[i]->C[0] - ObjTriMesh.point(v_it).data()[0];
			cur_view_direction[1] = listViewPoints[i]->C[1] - ObjTriMesh.point(v_it).data()[1];
			cur_view_direction[2] = listViewPoints[i]->C[2] - ObjTriMesh.point(v_it).data()[2];
			matrixProject(ObjTriMesh.point(v_it).data(), listViewPoints[i]->P,ox,1);
			if (matrixBound((int)(ox[0]+0.5), 0, listViewPoints[i]->height) && matrixBound((int)(ox[1]+0.5), 0, listViewPoints[i]->width)) {
				if(listViewPoints[i]->sil[(int)(ox[0]+0.5)][(int)(ox[1]+0.5)] == true)
				{
					//cout << "View Vector " << i << ": " << cur_view_direction[0] << "," << cur_view_direction[1] << "," << cur_view_direction[2] << ":";  
					//cout << "Inner Product : " << dotProduct(ObjTriMesh.normal(v_it).data(), cur_view_direction, 3) << ";";
					double output_angle = calAngle(ObjTriMesh.normal(v_it).data(), cur_view_direction);
					//cout << "Angle:" << output_angle << endl;
					if (output_angle<150) 
					{
						pair<int, double> tpair = pair<int, double>(i, CV_IMAGE_ELEM(listViewPoints[i]->image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)));
						tVectexIntensityList.push_back(make_pair(output_angle, tpair));
					}
				}
			}
		}
		if (tVectexIntensityList.size() < 1) {
			cout<< "The vertex can not find any visible color." << endl;
		}
		tVectexIntensityList.sort(MyDataSortPredicate);
		//vertex color is the average of the first min(tvil.size, 3) records

		//cout << ObjTriMesh.normal(v_it) << ": " << endl;
		nx = ObjTriMesh.normal(v_it).data()[0]; ny = ObjTriMesh.normal(v_it).data()[1]; nz = ObjTriMesh.normal(v_it).data()[2];
		toSpherical(nx, ny, nz, &z1, &z2);
		p_data[0*m_row + vertex_idx] = z1; p_data[1*m_row + vertex_idx] = z2; //record the two angle of sphere coordinates
		p_data[3*m_row + vertex_idx] = nx; p_data[4*m_row + vertex_idx] = ny; p_data[5*m_row + vertex_idx] = nz; //record the two angle of sphere coordinates

		std::list< pair<double, pair<int, double> > >::iterator tvi_iter;
		int i = 0; tcolor = 0.0; tcount = 0; temp_color = 0.0; temp_count = 0;
		for (tvi_iter = tVectexIntensityList.begin(); tvi_iter!=tVectexIntensityList.end(); tvi_iter ++) {
			pair<double, pair<int, double> > tpair = *tvi_iter;
			if (temp_count < 6) { // calculate average vertex color
				temp_color += tpair.second.second;
				temp_count ++;
			}
			if (tpair.first < 60) { // for data fitting
				tcolor +=tpair.second.second;
				tcount ++;
			}
		}
		Vertex_Color_List[vertex_idx] = (temp_color/temp_count);
		if (tcount > 0) { // records the intensity at this sphere coordinates z1, z2
			p_data[2*m_row + vertex_idx] = temp_color/temp_count;//tcolor/tcount;
		} else {
			p_data[2*m_row + vertex_idx] = temp_color/temp_count;//-1000000; // flag for invalid sample
		}
		//VertexIntensityRecords.push_back(tVectexIntensityList);
	}

	engPutVariable(m_ep, "PNI", PointNormalIntensity);
	engEvalString(m_ep, "save('PointNormalIntensity.mat','PNI');");
	char tempbuffer[255];
	sprintf(tempbuffer,"px = NormalIntensityFitting('PointNormalIntensity.mat',%d);", FittingOrder);
	engEvalString(m_ep, tempbuffer);
	sprintf(tempbuffer,"save('px_%d.mat','px');", FittingOrder);
	engEvalString(m_ep, tempbuffer);

	mxDestroyArray(PointNormalIntensity);

	
	fstream fout(resultfile.c_str(), ios::out);
	if (!fout) {
		cout << "Can not open " << resultfile << " to save data..." << endl;
	}
	fout << Vertex_Color_List.size() << endl;
	for (int i = 0; i < Vertex_Color_List.size(); i ++) {
		fout << Vertex_Color_List[i] << endl;
	}
	fout << endl;
	fout.close();
	cout << "Done ...." << endl;
	return;
}

void RecoverVertexPhotometricNormal(int step = 0)
{
	cout << "Start to recover vertex photometric normal...";
	sprintf(buffer, "_%.2d", step);
	string resultfile = (MOptions.DirName+"VertexPSNormal"+string(buffer)+".txt");
	if (FileExisted(resultfile.c_str())) {
		return;
	}
	Vertex_PSNormal_list.clear();

	//for each vertex, calculate its photometric normal according to its color and initial normal
	double nx, ny, nz, z1, z2; char tempbuffer[255]; 
	mxArray* ps_res = NULL; mxArray* f_val = NULL;
	sprintf(tempbuffer,"load('px_%d.mat','px');", FittingOrder);
	engEvalString(m_ep, tempbuffer); 
	mxArray* QueryPointNormalIntensity = mxCreateDoubleMatrix(ObjTriMesh.n_vertices(), 3, mxREAL);
	size_t m_row = mxGetM(QueryPointNormalIntensity);  size_t m_col = mxGetN(QueryPointNormalIntensity);	double* p_data = mxGetPr(QueryPointNormalIntensity);  
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		int vertex_idx = v_it.handle().idx();
		nx = ObjTriMesh.normal(v_it).data()[0]; ny = ObjTriMesh.normal(v_it).data()[1]; nz = ObjTriMesh.normal(v_it).data()[2];
		toSpherical(nx, ny, nz, &z1, &z2);
		p_data[0*m_row + vertex_idx] = z1; p_data[1*m_row + vertex_idx] = z2; //record the two angle of sphere coordinates
		p_data[2*m_row + vertex_idx] = Vertex_Color_List[vertex_idx];
	}
	engPutVariable(m_ep, "QPNI", QueryPointNormalIntensity);
	engEvalString(m_ep, "save('QueryPointNormalIntensity.mat','QPNI');");
	engEvalString(m_ep, "delete('RecoveredPointNormalIntensity.mat');");
	engEvalString(m_ep, "RecoverAllPhotometricNormalMultiLevel(px, 'QueryPointNormalIntensity.mat', 'RecoveredPointNormalIntensity.mat');");
	//sprintf(tempbuffer,"RecoverPhotometricNormal(px, 'QueryPointNormalIntensity.mat', 'RecoveredPointNormalIntensity.mat');");
	while(!FileExisted("c:\\MATLAB\\RecoveredPointNormalIntensity.mat")) 
	{
		::Sleep(10000); cout << "Wait 1000 ms for the photometric normal recovery process" << endl; 
	}

	Vertex_PSNormal_list.clear(); Vertex_PSNormal_list.resize(ObjTriMesh.n_vertices());
	engEvalString(m_ep, "load('RecoveredPointNormalIntensity.mat');");
	ps_res = engGetVariable(m_ep, "RPNI");
	double *rps_data = (double*)mxGetData(ps_res);
	for (int i = 0; i < m_row; i ++) {
		z1 = rps_data[0*m_row + i] ; z2 = rps_data[1*m_row + i] = z2; //record the two angle of sphere coordinates
		to3DNormal(z1, z2, &nx, &ny, &nz);
		vector<double> trps_vec; trps_vec.clear();
		trps_vec.push_back(nx); trps_vec.push_back(ny);  trps_vec.push_back(nz);
		Vertex_PSNormal_list[i] = (trps_vec);
	}
	mxDestroyArray(QueryPointNormalIntensity);
	mxDestroyArray(ps_res);
	
	fstream fout(resultfile.c_str(), ios::out);
	if (!fout) {
		cout << "Can not open " << resultfile << " to save data..." << endl;
	}
	fout << Vertex_PSNormal_list.size() << endl;
	for (int i = 0; i < Vertex_PSNormal_list.size(); i ++) {
		fout << Vertex_PSNormal_list[i][0] << " " << Vertex_PSNormal_list[i][1] << " "<< Vertex_PSNormal_list[i][2] << endl;
	}
	fout.close();
	cout << "Done..."<<endl;
}

bool SaveSplineToFile(Nag_2dSpline &inSpline, const char* outfilename)
{
	fstream fout(outfilename, ios::out);
	if (!fout) {
		cout << "Can not open " << outfilename << " to save data..." << endl;
		 return false;
	}
	fout << inSpline.nx << "   " << inSpline.ny << endl;
	/* Store the knots spline.lamda[0]...spline.lamda[nx-1]
	 * and spline.mu[0]...spline.mu[ny-1].
	 */
	for (int i = 0; i < inSpline.nx;  i++) {
		fout << inSpline.lamda[i] << "  ";
	}	fout << endl;
	for (int i = 0; i < inSpline.ny;  i++) {
		fout << inSpline.mu[i] << "  ";
	}	fout << endl;
  /* Store spline.c, the bicubic spline coefficients. */
	for (int i = 0; i < (inSpline.nx-4)*(inSpline.ny-4); i++) {
		fout << inSpline.c[i] << "  ";
	}	fout << endl;
	fout.close();
	cout << "Store spline function into file " << outfilename << " Done..." << endl;
	return true;
}

bool LoadSplineFromFile(Nag_2dSpline &outSpline, const char* infilename)
{
	fstream fin(infilename, ios::in);
	if (!fin) {
		cout << "Can not open " << infilename << " to save data..." << endl;
		 return false;
	}
	fin >> outSpline.nx >> outSpline.ny;
	if (!(outSpline.c = NAG_ALLOC((outSpline.nx-4)*(outSpline.ny-4), double)) ||
		!(outSpline.lamda = NAG_ALLOC(outSpline.nx, double)) ||
		!(outSpline.mu = NAG_ALLOC(outSpline.ny, double)))
	{
		cout << "Storage allocation failed." << endl;
		return false;
	}
	/* Read the knots spline.lamda[0]...spline.lamda[nx-1]
	 * and spline.mu[0]...spline.mu[ny-1].
	 */
	for (int i = 0; i < outSpline.nx;  i++) {
		fin >> outSpline.lamda[i];
	}
	for (int i = 0; i < outSpline.ny;  i++) {
		fin >> outSpline.mu[i];
	}
  /* Read spline.c, the bicubic spline coefficients. */
	for (int i = 0; i < (outSpline.nx-4)*(outSpline.ny-4); i++) {
		fin >> outSpline.c[i];
	}
	fin.close();
	cout << "Load spline function from file " << infilename << " Done..." << endl;
	return true;
}

bool CubicSplineFittingData(const char* datafilename, const char* outsplinefilename, int fittingmethod = 1, double range_value = 1.0)
{
	// fitting the cubic spline function and store the final result in file
	// start the nag c cubic spline fitting process
	NagError fail;  INIT_FAIL(fail);	Nag_Start start = Nag_Cold;		Nag_2dSpline splinex;
	double fpcx, s, warmstartinf, xhi, xlo, yhi, ylo;	Integer data_size, nxest, nyest, cxrank;
	double *in_z1=0, *in_z2=0, *in_int=0, *weights=0, *fit_int = 0;
	// load temp data of normal and intensity
	fstream tfin(datafilename, ios::in);
	if (!tfin) {
		cout << "Can not open " << datafilename << " to load data..." << endl;		return false;
	}
	tfin >> data_size;
	if(data_size>=16) {
		if(	!(in_int = NAG_ALLOC(data_size, double)) || !(weights = NAG_ALLOC(data_size, double)) || !(fit_int = NAG_ALLOC(data_size, double)) ||
			!(in_z1 = NAG_ALLOC(data_size, double))	 ||	!(in_z2 = NAG_ALLOC(data_size, double)) ) {
				cout<<"Allocation input data memory failure"<<endl;		return false;
		}
	}
	else {
		cout<<"Invalid size of input data: "<<data_size<<endl;		return false;
	}
	xlo = 1000; xhi = -1000; ylo = xlo; yhi = xhi;
	double iz1, iz2, iinten, wei; 
	for(size_t num = 0; num < data_size; num++) {
		tfin>>iz1>>iz2>>iinten>>wei;
		xlo = min(xlo, iz1-0.0001); xhi = max(xhi, iz1+0.0001); ylo = min(ylo, iz2-0.0001); yhi = max(yhi, iz2+0.0001);
		in_z1[num] = iz1; in_z2[num] = iz2; in_int[num] = iinten; weights[num] = wei; 
	}
	tfin.close();

	////////////////////////////////////////////////////////////////////////////////
	/* Initialise spline */
	splinex.lamda = 0; splinex.mu = 0; splinex.c = 0; 
	/* Scalars */
	double       sigma, sum, eps = 1e-10; // about 180*0.2/pi = 11.5 degree
	Integer      iadres, nc, np, npoint, pz1, pz2, rank;
	double pres = s = 1e35; //smooth factor, the upper bounder of error, initial status
	/* Arrays */
	double       *dl = 0, *lamda = 0, *mu = 0;	Integer      *point = 0;

	fstream tfout;
	cout << "Spline surface fitting using method "<< fittingmethod << " and range value "<<range_value<<endl;
	switch (fittingmethod) {
	case 0: // nag_2d_spline_fit_scat, but can not get the right fitting result, maybe the data size is too large
		nxest = (long int)std::floor(std::sqrt((double)data_size)/5); nyest = (long int)std::floor(std::sqrt((double)data_size)/5); // upper bound of the number of knots
		start = Nag_Cold;
		// start to fit the surface of intensity by cubic spline
		INIT_FAIL(fail);
		nag_2d_spline_fit_scat(start, data_size, in_z1, in_z2, in_int, weights, s, nxest, nyest, 
			&fpcx, &cxrank, &warmstartinf, &splinex, &fail);
		cout<<"The output knot points number are: "<<splinex.nx<<" , "<<splinex.ny<<" for data point size: "<<data_size<<endl;
		cout<<"Rank deficiency cx = "<<(splinex.nx-4)*(splinex.ny-4)-cxrank<<endl;
		cout<<"Average Squared residuals is fpcx = "<<fpcx/data_size<<endl;
		if(splinex.nx == 8 && splinex.ny == 8) {
			cout<<"The spline for cx is the least squares bi-cubic polynomial"<<endl;
		}
		if (fail.code != NE_NOERROR) {
			cout<<"Error from nag_2d_spline_fit_scat in cx. "<<fail.message<<endl;	return false;
		} else {	cout << ".";
			while(fail.code == NE_NOERROR) {
				pres = s; s = fpcx * SplineFittingDecreasePara; 
				nag_2d_spline_fit_scat(start, data_size, in_z1, in_z2, in_int, weights, s, nxest, nyest, 
					&fpcx, &cxrank, &warmstartinf, &splinex, &fail);	cout << ".";
			}	cout << ".";
		}
		cout << endl << "The spline fitting parameter is: " << pres << endl;
		nag_2d_spline_fit_scat(start, data_size, in_z1, in_z2, in_int, weights, pres, nxest, nyest, 
			&fpcx, &cxrank, &warmstartinf, &splinex, &fail);
		cout<<"The output knot points number are: "<<splinex.nx<<" , "<<splinex.ny<<" for data point size: "<<data_size<<endl;
		cout<<"Rank deficiency cx = "<<(splinex.nx-4)*(splinex.ny-4)-cxrank<<endl;
		cout<<"Average Squared residuals is fpcx = "<<fpcx/data_size<<endl;
		if(splinex.nx == 8 && splinex.ny == 8) {
			cout<<"The spline for cx is the least squares bi-cubic polynomial"<<endl;
		}
		break;
	case 1: // nag_2d_spline_fit_panel and nag_2d_panel_sort
		pz1 = 8+std::floor((xhi - xlo)/range_value); 
		pz2 = 8+std::floor((yhi - ylo)/range_value);
		nc = (pz1 - 4) * (pz2 - 4);
		np = (pz1 - 7) * (pz2 - 7);
		npoint = data_size+(pz1-7)*(pz2-7);
		/* Allocate memory */
		if (!(dl = NAG_ALLOC(nc, double))		|| !(point = NAG_ALLOC(npoint, Integer)) || 
			!(lamda = NAG_ALLOC(pz1, double))	|| !(mu = NAG_ALLOC(pz2, double)) ) {
			cout<< "Allocation failure" << endl;	return false;
		}
		for (int i = 5; i <= pz2 - 4; ++ i) {
			mu[i-1] = min(yhi, ylo + range_value * (i-4)); 
		} //mu[0] = mu[1] = mu[2] = mu[3] = -DQPI;mu[pz2-4] = mu[pz2-3] = mu[pz2-2] = mu[pz2-1] = DQPI;
		SaveArray(mu, pz2, "mu.txt");
		for (int i = 5; i <= pz1 - 4; ++ i) {
			lamda[i-1] = min(xhi, xlo + range_value * (i-4));
		} //lamda[0] = lamda[1] = lamda[2] = lamda[3] = 0.00; lamda[pz1-4] = lamda[pz1-3] = lamda[pz1-2] = lamda[pz1-1] = DQPI;
		SaveArray(lamda, pz1, "lambda.txt");
		nag_2d_panel_sort(pz1, pz2, lamda, mu, data_size, in_z1, in_z2, point, &fail);
		if (fail.code != NE_NOERROR) {
			cout << endl << "Error from nag_2d_panel_sort: " << fail.message << endl; return false;
		}	SaveArray(point, npoint, "point.txt");
		splinex.nx = pz1; splinex.ny = pz2;
		if (!(splinex.c = NAG_ALLOC((splinex.nx-4)*(splinex.ny-4), double)) ||
			!(splinex.lamda = NAG_ALLOC(splinex.nx, double)) || !(splinex.mu = NAG_ALLOC(splinex.ny, double))) {
				cout<< "Allocation failure" << endl;	return false;
		}
		for (int i = 0; i < splinex.nx; i++) {
			splinex.lamda[i] = lamda[i];
		}
		for (int i = 0; i < splinex.ny;  i++) {
			splinex.mu[i] = mu[i];
		}
		/* nag_2d_spline_fit_panel (e02dac).
         * Least-squares surface fit, bicubic splines
         */
		nag_2d_spline_fit_panel(data_size, in_z1, in_z2, in_int, weights, point, dl, eps, &sigma, &rank, &splinex, &fail);
		if (fail.code != NE_NOERROR) {
			cout<<endl<<"Error from nag_2d_spline_fit_panel in cx. "<<fail.message<<endl;	
			return false;
		}
		cout << endl << "Rank is: " << rank << "; ";
		if (rank == (splinex.nx-4)*(splinex.ny-4)) {
			cout << "the unique solution is obtained; " << endl;
		}
		cout <<"The fitting error is: " << sigma << endl;
		SaveArray(splinex.mu, splinex.ny, "splinemu.txt"); SaveArray(splinex.lamda, splinex.nx, "splinelambda.txt");
		break;
	default:
		break;
	}
#ifdef _DEBUG
	// Test the fitting result, store in file;
	nag_2d_spline_eval(data_size, in_z1, in_z2, fit_int, &splinex, &fail);
	if (fail.code != NE_NOERROR) {
		cout<<"Error from nag_2d_spline_eval in cx. "<<fail.message<<endl;	
		return false;
	}
	tfout.open("TestCubicSplineFitting.txt", ios::out);
	if (!tfout) {
		cout << "Can not open " << "TestCubicSplineFitting.txt" << " to save data..." << endl;		return false;
	}
	for (int i = 0; i < data_size; ++i)
	{
		tfout << in_z1[i] << "   " << in_z2[i] << "   " << in_int[i] << "   " << fit_int[i] << "   " << std::abs(in_int[i]-fit_int[i]) << endl;
	}
	tfout.close();
#endif
	if (in_z1) 			NAG_FREE(in_z1);
	if (in_z2) 			NAG_FREE(in_z2);
	if (in_int)			NAG_FREE(in_int);
	if (fit_int)		NAG_FREE(fit_int);
	if (weights)		NAG_FREE(weights);
	// finish cubic spline fitting for normal direction z1, z2 and its intensity

	if (SaveSplineToFile(splinex, outsplinefilename)) {
		return true;
	}
	return false;
}

bool PointInsideTriangle(int* P, int* A, int* B, int* C)
{
	double epsilon = 0.0001;
	double AB_A = (B[0]-A[0])/(A[1]-B[1]+epsilon), AB_C = -(A[0]+AB_A*A[1]);
	double BC_A = (C[0]-B[0])/(B[1]-C[1]+epsilon), BC_C = -(B[0]+BC_A*B[1]);
	double AC_A = (C[0]-A[0])/(A[1]-C[1]+epsilon), AC_C = -(A[0]+AC_A*A[1]);
	double PAJudge = (P[0]+BC_A*P[1]+BC_C)*(A[0]+BC_A*A[1]+BC_C); // judge the result of both P and A to line BC
	double PBJudge = (P[0]+AC_A*P[1]+AC_C)*(B[0]+AC_A*B[1]+AC_C); // judge the result of both P and B to line AC
	double PCJudge = (P[0]+AB_A*P[1]+AB_C)*(C[0]+AB_A*C[1]+AB_C); // judge the result of both P and C to line AB
	int num_zero = 0, num_pos = 0, num_neg = 0;
	PAJudge>0?num_pos++:(PAJudge==0?num_zero++:num_neg++);
	PBJudge>0?num_pos++:(PBJudge==0?num_zero++:num_neg++);
	PCJudge>0?num_pos++:(PCJudge==0?num_zero++:num_neg++);
	if (num_neg == 0) {
		return true;
	} else if (num_zero == 2) {
		return true;
	} else {
		return false;
	}
}

bool CalculateVertexVisibleToViewPoint(const char* resultfile)
{
	//calculate the visible result of each vertex to viewpoint
	ObjTriMesh.request_face_normals(); ObjTriMesh.update_face_normals();
	vector< vector<int> > vertex_visible_view; vertex_visible_view.resize(ObjTriMesh.n_vertices());
	for (int i = 0; i < ObjTriMesh.n_vertices(); ++i) {
		vertex_visible_view[i].clear();
	}
	vector< vector<int> > visible_fid; vector< vector<double> > nr_distance; double ox[2];
	for (int idx = 0; idx < nViews; ++ idx) {
		OpenMesh::Vec3f Camera_center(listViewPoints[idx]->C[0], listViewPoints[idx]->C[1], listViewPoints[idx]->C[2]);
		if (visible_fid.size() != listViewPoints[idx]->height) {
			visible_fid.resize(listViewPoints[idx]->height);			nr_distance.resize(listViewPoints[idx]->height);
		}
		for (int j = 0; j < listViewPoints[idx]->height; ++j) {
			if (visible_fid[j].size() != listViewPoints[idx]->width) {
				visible_fid[j].resize(listViewPoints[idx]->width);	nr_distance[j].resize(listViewPoints[idx]->width);
			}
			for (int k = 0; k < listViewPoints[idx]->width; ++k) {
				visible_fid[j][k] = -1; nr_distance[j][k] = 1.0e10;
			}
		} // end of initialization

		//cv::Mat silImage; silImage.create(listViewPoints[idx]->height, listViewPoints[idx]->width, CV_8UC1); silImage.setTo(0);
		//int temp_count = 0;
		for (MyMesh::FaceIter f_it = ObjTriMesh.faces_begin(); f_it != ObjTriMesh.faces_end(); ++ f_it) {
			MyMesh::ConstFaceVertexIter cfv_it = ObjTriMesh.cfv_iter(f_it);  int A[2], B[2], C[2], P[2];
			OpenMesh::Vec3f p1 = ObjTriMesh.point(cfv_it); ++cfv_it;
			matrixProject(p1.data(), listViewPoints[idx]->P,ox,1); int ph1 = (int)(ox[0]+0.5); int pw1 = (int)(ox[1]+0.5);	A[1] = ph1; A[0] = pw1;
			double dis1 = (Camera_center - p1).length();
			OpenMesh::Vec3f p2 = ObjTriMesh.point(cfv_it); ++cfv_it;
			matrixProject(p2.data(), listViewPoints[idx]->P,ox,1); int ph2 = (int)(ox[0]+0.5); int pw2 = (int)(ox[1]+0.5);	B[1] = ph2; B[0] = pw2;
			double dis2 = (Camera_center - p2).length();
			OpenMesh::Vec3f p3 = ObjTriMesh.point(cfv_it);
			matrixProject(p3.data(), listViewPoints[idx]->P,ox,1); int ph3 = (int)(ox[0]+0.5); int pw3 = (int)(ox[1]+0.5);	C[1] = ph3; C[0] = pw3;
			double dis3 = (Camera_center - p3).length();
			double avg_dis = min(min(dis1,dis2),dis3);	int visible_flag = true;  
			int max_ph = max(max(ph1, ph2), ph3); int max_pw = max(max(pw1, pw2), pw3);
			int min_ph = min(min(ph1, ph2), ph3); int min_pw = min(min(pw1, pw2), pw3);

			if (OpenMesh::dot(ObjTriMesh.normal(f_it), Camera_center - p1)<0) {
				continue;
			}

			for (int h = max(1,min_ph-1); h <= min(max_ph+1, listViewPoints[idx]->height-1); ++ h) {
				for (int w = max(1,min_pw-1); w <= min(max_pw+1, listViewPoints[idx]->width-1); ++ w) {
					P[1] = h; P[0] = w;
					if (PointInsideTriangle(P, A, B, C)) {
						if (nr_distance[h][w] > avg_dis) {
							nr_distance[h][w] = avg_dis;
							visible_fid[h][w] = f_it.handle().idx();
							//silImage.at<unsigned char>(h,w) = 255;
						}
					}
				}
			}
		} //end of visible result calculation of vertex to one viewpoint
		//imwrite("test2.png", silImage);

		for (int h = 0; h < visible_fid.size(); ++ h) {
			for (int w = 0; w < visible_fid[h].size(); ++ w) {
				int fid = visible_fid[h][w];
				if (fid != -1) {
					MyMesh::ConstFaceVertexIter cfv_it = ObjTriMesh.cfv_iter(MyMesh::FaceHandle(fid));
					int ida = cfv_it.handle().idx(); ++cfv_it;
					int idb = cfv_it.handle().idx(); ++cfv_it;
					int idc = cfv_it.handle().idx();
					if (vertex_visible_view[ida].size() < 1 || (vertex_visible_view[ida].size() > 0&&vertex_visible_view[ida][vertex_visible_view[ida].size()-1] != idx)) {
						vertex_visible_view[ida].push_back(idx);
					}
					if (vertex_visible_view[idb].size() < 1 || (vertex_visible_view[idb].size() > 0&&vertex_visible_view[idb][vertex_visible_view[idb].size()-1] != idx)) {
						vertex_visible_view[idb].push_back(idx);
					}
					if (vertex_visible_view[idc].size() < 1 || (vertex_visible_view[idc].size() > 0&&vertex_visible_view[idc][vertex_visible_view[idc].size()-1] != idx)) {
						vertex_visible_view[idc].push_back(idx);
					}
				}
			}
		}//record the visible viewpoint result to vertex	
	}

	fstream fout(resultfile, ios::out);
	if (!fout) {
		cout << "Can not open " << resultfile << " to save result." << endl;
		return false;
	}
	for (int i = 0; i < vertex_visible_view.size(); ++ i) {
		fout << vertex_visible_view[i].size()<<" ";
		for (int j = 0; j < vertex_visible_view[i].size(); ++ j) {
			fout << vertex_visible_view[i][j] << " ";
		}
		fout << endl;
	}
	fout.close();
 	return true;
}


bool UpdateMeshVertexIntensity(const char* vertexvisiblefile, const char* vertexintensityfilename)
{
	//load pre-calculate vertex visible viewpoint result
	vector< vector<int> > vertex_visible_view; vertex_visible_view.resize(ObjTriMesh.n_vertices());
	for (int i = 0; i < ObjTriMesh.n_vertices(); ++i) {
		vertex_visible_view[i].clear();
	}
	fstream fin(vertexvisiblefile, ios::in);
	if (!fin) {
		cout << "Can not open " << vertexvisiblefile << " to load result." << endl;
		return false;
	}
	int num, view_id;
	for (int i = 0; i < ObjTriMesh.n_vertices(); ++ i) {
		fin >> num;
		for (int j = 0; j < num; ++ j) {
			fin>>view_id; vertex_visible_view[i].push_back(view_id);
		}
	}
	fin.close();

	ObjTriMesh.update_face_normals(); ObjTriMesh.update_vertex_normals(); double ox[2];
	//construct the initial normal/light matrix, to learn their relationship
	fstream fout; fout.open(vertexintensityfilename, ios::out);
	if (!fout) {
		cout << "Can not open " << vertexintensityfilename << " to save data..." << endl;
	}
	fout << ObjTriMesh.n_vertices() << endl;
	std::list< pair<double, pair<int, double> > > tVectexIntensityList;   vector<OpenMesh::Vec3f> Camera_center_list;
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		tVectexIntensityList.clear();	Camera_center_list.clear();		int vertex_idx = v_it.handle().idx();
		for (int i = 0; i < vertex_visible_view[vertex_idx].size(); ++ i) {
			int view_id = vertex_visible_view[vertex_idx][i];
			OpenMesh::Vec3f Camera_center(listViewPoints[view_id]->C[0], listViewPoints[view_id]->C[1], listViewPoints[view_id]->C[2]);
			OpenMesh::Vec3f View_vector = Camera_center - ObjTriMesh.point(v_it);
			if (OpenMesh::dot(ObjTriMesh.normal(v_it), View_vector)>0) {
				double output_angle = calAngle(ObjTriMesh.normal(v_it).data(), View_vector.data());
				matrixProject(ObjTriMesh.point(v_it).data(), listViewPoints[view_id]->P,ox,1);
				pair<int, double> tpair = pair<int, double>(view_id, CV_IMAGE_ELEM(listViewPoints[view_id]->image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)));
				tVectexIntensityList.push_back(make_pair(output_angle, tpair));
				Camera_center_list.push_back(Camera_center);
			}
		}
		tVectexIntensityList.sort(MyDataSortPredicate); int vi_size = tVectexIntensityList.size();
		std::list< pair<double, pair<int, double> > >::iterator tvi_iter = tVectexIntensityList.begin();  
		pair<double, pair<int, double> > tpair;  int tsize = min(vi_size, 20);
		fout<<tsize<<" ";
		for (int k = 0; k < tsize; tvi_iter ++, k ++) {
			tpair = *tvi_iter;
			fout<<tpair.second.second<<" "<<Camera_center_list[k].data()[0]<<" "<<Camera_center_list[k].data()[1]<<" "<<Camera_center_list[k].data()[2]<<" ";
		}
		fout << endl;
	}
	fout.close();
	return true;
}

bool UpdateMeshVertexIntensity(const char* vertexintensityfilename)
{
	ObjTriMesh.update_face_normals(); ObjTriMesh.update_vertex_normals(); double ox[2];
	//construct the initial normal/light matrix, to learn their relationship
	fstream fout; fout.open(vertexintensityfilename, ios::out);
	if (!fout) {
		cout << "Can not open " << vertexintensityfilename << " to save data..." << endl;
	}
	fout << ObjTriMesh.n_vertices() << endl;
	std::list< pair<double, pair<int, double> > > tVectexIntensityList;		vector<OpenMesh::Vec3f> Camera_center_list;
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
		tVectexIntensityList.clear();	Camera_center_list.clear();		int vertex_idx = v_it.handle().idx();
		for (int i = 0; i < nViews; i ++) {
			OpenMesh::Vec3f Camera_center(listViewPoints[i]->C[0], listViewPoints[i]->C[1], listViewPoints[i]->C[2]);
			OpenMesh::Vec3f View_vector = Camera_center - ObjTriMesh.point(v_it);
			matrixProject(ObjTriMesh.point(v_it).data(), listViewPoints[i]->P,ox,1);
			if (matrixBound((int)(ox[0]+0.5), 0, listViewPoints[i]->height) && matrixBound((int)(ox[1]+0.5), 0, listViewPoints[i]->width)) {
				double output_angle = calAngle(ObjTriMesh.normal(v_it).data(), View_vector.data()); //listViewPoints[i]->InverseViewVector
				//cout << "View Vector " << i << ": " << cur_view_direction[0] << "," << cur_view_direction[1] << "," << cur_view_direction[2] << ":";  
				//cout << "Inner Product : " << dotProduct(ObjTriMesh.normal(v_it).data(), cur_view_direction, 3) << ";";
				//cout << "Angle:" << output_angle << endl;
				if(listViewPoints[i]->sil[(int)(ox[0]+0.5)][(int)(ox[1]+0.5)] == true && OpenMesh::dot(ObjTriMesh.normal(v_it), View_vector)>0)
				{
					pair<int, double> tpair = pair<int, double>(i, CV_IMAGE_ELEM(listViewPoints[i]->image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)));
					tVectexIntensityList.push_back(make_pair(output_angle, tpair));
					Camera_center_list.push_back(Camera_center);
				}
			}
		}
		tVectexIntensityList.sort(MyDataSortPredicate); int vi_size = tVectexIntensityList.size();
		std::list< pair<double, pair<int, double> > >::iterator tvi_iter = tVectexIntensityList.begin();  
		pair<double, pair<int, double> > tpair;  int tsize = min(vi_size, 20);
		fout<<tsize<<" ";
		for (int k = 0; k < tsize; tvi_iter ++, k ++) {
			tpair = *tvi_iter;
			fout<<tpair.second.second<<" "<<Camera_center_list[k].data()[0]<<" "<<Camera_center_list[k].data()[1]<<" "<<Camera_center_list[k].data()[2]<<" ";
		}
		fout << endl;
	}
	fout.close();
	return true;
}

bool UpdateVertexPSNormal(const char* filename)
{
	ObjTriMesh.update_face_normals(); ObjTriMesh.update_vertex_normals();
	// Save the output ps normal for vertex;
	fstream fout(filename, ios::out);
	if (!fout) {
		cout << "Can not open " << filename << " to save data..." << endl; 
		return false;
	}
	fout << ObjTriMesh.n_vertices() << endl; 
	for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin(); v_it != ObjTriMesh.vertices_end(); ++ v_it) {
		fout << ObjTriMesh.normal(v_it)[0] << "   " << ObjTriMesh.normal(v_it)[1] << "   " << ObjTriMesh.normal(v_it)[2] << endl;
	}
	fout << endl;
	fout.close();
	cout << "Done ...." << endl;
	return true;
}

bool CubicSplineFittingPSNormal(const char* vertexintensityfilename, const char* psnormalfilename, int step = 0, int choice = 1, double range_value = 0.5)
{
	ObjTriMesh.request_face_normals();		ObjTriMesh.update_face_normals();
	ObjTriMesh.request_vertex_normals();	ObjTriMesh.update_vertex_normals();

	cout << "Start to calculate the vertex intensity value ";
	sprintf(buffer, "_%.2d", step);
	string tempdatafile = (MOptions.DirName+"TempNormalColor"+string(buffer)+".txt");
	double timer_start = (double)cv::getTickCount();
	//if (!FileExisted(tempdatafile.c_str())) 
	{
		vector<double> vec_z1, vec_z2, vec_weight; int nocolor_point_num = 0;
		double ox[2], cur_view_direction[3]; int temp_count; double temp_color, nx, ny, nz, z1, z2, vertex_weight;
		//construct the initial normal/light matrix, to learn their relationship
		Vertex_Color_List.clear(); Vertex_Color_List.resize(ObjTriMesh.n_vertices());  
		vec_z1.clear(); vec_z1.resize(ObjTriMesh.n_vertices());	vec_z2.clear(); vec_z2.resize(ObjTriMesh.n_vertices()); 
		vec_weight.clear(); vec_weight.resize(ObjTriMesh.n_vertices());
		std::list< pair<double, pair<int, double> > > tVectexIntensityList;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
			tVectexIntensityList.clear();
			int vertex_idx = v_it.handle().idx();
			for (int i = 0; i < nViews; i ++) {
				cur_view_direction[0] = listViewPoints[i]->C[0] - ObjTriMesh.point(v_it).data()[0];
				cur_view_direction[1] = listViewPoints[i]->C[1] - ObjTriMesh.point(v_it).data()[1];
				cur_view_direction[2] = listViewPoints[i]->C[2] - ObjTriMesh.point(v_it).data()[2];
				matrixProject(ObjTriMesh.point(v_it).data(), listViewPoints[i]->P,ox,1);
				if (matrixBound((int)(ox[0]+0.5), 0, listViewPoints[i]->height) && matrixBound((int)(ox[1]+0.5), 0, listViewPoints[i]->width)) {
					double output_angle = calAngle(ObjTriMesh.normal(v_it).data(), cur_view_direction); //listViewPoints[i]->InverseViewVector
					//cout << "View Vector " << i << ": " << cur_view_direction[0] << "," << cur_view_direction[1] << "," << cur_view_direction[2] << ":";  
					//cout << "Inner Product : " << dotProduct(ObjTriMesh.normal(v_it).data(), cur_view_direction, 3) << ";";
					//cout << "Angle:" << output_angle << endl;
					if(listViewPoints[i]->sil[(int)(ox[0]+0.5)][(int)(ox[1]+0.5)] == true && output_angle<150)
					{
						pair<int, double> tpair = pair<int, double>(i, CV_IMAGE_ELEM(listViewPoints[i]->image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)));
						tVectexIntensityList.push_back(make_pair(output_angle, tpair));
					}
				}
			}
			if (tVectexIntensityList.size() < 1) {
				//cout<< "The vertex can not find any visible color." << endl;
				nocolor_point_num ++;
				vec_z1[vertex_idx] = z1;		vec_z2[vertex_idx] = z2;
				Vertex_Color_List[vertex_idx] = -1;
				continue;
			}
			tVectexIntensityList.sort(MyDataSortPredicate);
			//vertex color is the average of the first min(tvil.size, 3) records
			//cout << ObjTriMesh.normal(v_it) << ": " << endl;
			toSpherical(ObjTriMesh.normal(v_it).data()[0], ObjTriMesh.normal(v_it).data()[1], ObjTriMesh.normal(v_it).data()[2], &z1, &z2);
			std::list< pair<double, pair<int, double> > >::iterator tvi_iter = tVectexIntensityList.begin();
			pair<double, pair<int, double> > tpair = *tvi_iter;			double avg_angle = tpair.first;
			int i = 0; temp_color = 0.0; temp_count = 0;   double sum_valid_angle = 0.0;
			for (tvi_iter = tVectexIntensityList.begin(); tvi_iter!=tVectexIntensityList.end(); tvi_iter ++) {
				tpair = *tvi_iter;
				if (abs(tpair.first - avg_angle) < 3) {
					temp_color += tpair.second.second; temp_count ++;
					sum_valid_angle += tpair.first;  
					avg_angle = sum_valid_angle/temp_count;
				} else
					break;
			}
			vec_z1[vertex_idx] = z1;		vec_z2[vertex_idx] = z2;  //vec_weight[vertex_idx] = vertex_weight;
			Vertex_Color_List[vertex_idx] = (temp_color/temp_count);
		}
		cout << nocolor_point_num << " out of " << ObjTriMesh.n_vertices() << " don't find any color." << endl;
		double max_weight = -100;
		for (MyMesh::VertexIter v_it = ObjTriMesh.vertices_begin();v_it != ObjTriMesh.vertices_end(); ++v_it) {
			// calculate the weight for each sample, using the difference between neighbor vertices
			double avg_nei_intensity = 0.0; int nei_count = 0;
			for (MyMesh::ConstVertexVertexIter cvv_it = ObjTriMesh.vv_iter(v_it.handle()); cvv_it; ++ cvv_it) {
				avg_nei_intensity += Vertex_Color_List[cvv_it.handle().idx()]; nei_count ++;
			}
			avg_nei_intensity = avg_nei_intensity/nei_count;
			if (Vertex_Color_List[v_it.handle().idx()] == -1) {
				Vertex_Color_List[v_it.handle().idx()] = avg_nei_intensity;
			}
			vec_weight[v_it.handle().idx()] = std::exp(-std::pow((Vertex_Color_List[v_it.handle().idx()] - avg_nei_intensity),2.0)/50);
			if (vec_weight[v_it.handle().idx()] > max_weight) {
				max_weight = vec_weight[v_it.handle().idx()];
			}
		}
		for (int i = 0; i < vec_weight.size(); ++ i) {
			vec_weight[i] = vec_weight[i]/max_weight;
		}
		// store data into file for nagc cubic spline fitting
		fstream tfout(tempdatafile.c_str(), ios::out);
		if (!tfout) {
			cout << "Can not open " << tempdatafile << " to save data..." << endl;
			return false;
		}
		tfout << Vertex_Color_List.size() << endl;
		for (int i = 0; i < Vertex_Color_List.size(); i ++) {
			tfout << vec_z1[i] << "  " << vec_z2[i] << "  " << Vertex_Color_List[i] << "  " << vec_weight[i] << endl;
		}
		tfout.close();
		printf("\nInitial normal intensity extraction time = %lfs ",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());
	} cout << "Done..." << endl;

	// load the normal and intensity file
	Integer data_size; double *in_z1=0, *in_z2=0, *in_int=0, *weights=0;
	// load temp data of normal and intensity
	fstream tfin(tempdatafile.c_str(), ios::in);
	if (!tfin) {
		cout << "Can not open " << tempdatafile << " to load data..." << endl;		return false;
	}
	tfin >> data_size;
	if(data_size>=16) {
		if(	!(in_int = NAG_ALLOC(data_size, double)) || !(weights = NAG_ALLOC(data_size, double)) ||
			!(in_z1 = NAG_ALLOC(data_size, double))	 ||	!(in_z2 = NAG_ALLOC(data_size, double)) ) {
				cout<<"Allocation input data memory failure"<<endl;		return false;
		}
	}
	else {
		cout<<"Invalid size of input data: "<<data_size<<endl;		return false;
	}
	double xlo = 1000, xhi = -1000, ylo = xlo, yhi = xhi;
	double iz1, iz2, iinten, wei; 
	for(size_t num = 0; num < data_size; num++) {
		tfin>>iz1>>iz2>>iinten>>wei;
		xlo = min(xlo, iz1-0.0001); xhi = max(xhi, iz1+0.0001); ylo = min(ylo, iz2-0.0001); yhi = max(yhi, iz2+0.0001);
		in_z1[num] = iz1; in_z2[num] = iz2; in_int[num] = iinten; weights[num] = wei; 
	}
	tfin.close();

	sprintf(buffer, "_%.2d", step);
	string tempsplineresultfile = (MOptions.DirName+"TempSplineResult"+string(buffer)+".txt");
	timer_start = (double)cv::getTickCount();
	CubicSplineFittingData(tempdatafile.c_str(), tempsplineresultfile.c_str(), choice, range_value);
	printf("\nCubic spline calculation Time = %lfs\n",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());
	
	NagError fail;  INIT_FAIL(fail);
	Nag_2dSpline splinex;	LoadSplineFromFile(splinex, tempsplineresultfile.c_str());
	xlo = max(xlo, splinex.lamda[3]);	xhi = min(xhi, splinex.lamda[splinex.nx-4]);
	ylo = max(ylo, splinex.mu[3]);		yhi = min(yhi, splinex.mu[splinex.ny-4]);
	// start the fitting process of z1, z2 according to the intensity, in multi-level scheme for acceleration
	double *recx = 0, *recy = 0, *fgx = 0; 
	double coarse_unit = 0.0005, fine_unit = 0.0001;
	int coarse_range = 300, fine_range = 30, coarse_diff = 5, fine_diff = 3; double fine_step = 0.3;
	vector<double> Fitted_z1, Fitted_z2, Fitted_intensity; Fitted_z1.clear(); Fitted_z2.clear(); Fitted_intensity.clear(); // the final fitting results
	// start fitting in the coarse level
	int npx = std::ceil((xhi - xlo)/coarse_unit)+1; int npy = std::ceil((yhi - ylo)/coarse_unit)+1;
	//construct the result fitted data for nag_2d_spline_eval_rect 
	if (npx >= 1 && npy >= 1) {
		if (!(recx = NAG_ALLOC(npx,double)) || !(recy = NAG_ALLOC(npy,double)) || !(fgx = NAG_ALLOC(npx*npy,double)) ) {
			cout<<"Allocation output file memory failure"<<endl;	return false;
		}
	}
	else {
		cout<<"Invalid npx or npy"<<endl; return false;
	}
	for (int i = 0; i < npx; i++) {
		recx[i] = min(xhi, xlo + coarse_unit*i);
	} SaveArray(recx, npx, "recx.txt");
	for (int i = 0; i < npy; i++) {
		recy[i] = min(yhi, ylo + coarse_unit*i);
	} SaveArray(recy, npy, "recy.txt");
	//end of construction
	//start of nag_2d_spline_eval_rect calculation, but it will cause error when error point at the boundary of this region.
	nag_2d_spline_eval_rect(npx, npy, recx, recy, fgx, &splinex, &fail);
	if (fail.code != NE_NOERROR) {
		cout<<"Error from nag_2d_spline_eval_rect in cx. "<<fail.message<<endl; return false;
	}

	timer_start = (double)cv::getTickCount();
	int invalid_count = 0; double max_fitting_int_diff = -100;
	for (int idx = 0; idx < data_size; ++ idx) {
		// for each input intensity and normal, find the new one with the smallest intensity difference in certain range +/-12 degree
		double min_diff = 10000, out_z1, out_z2, out_int;
		int i_min = std::floor((in_z2[idx]-ylo)/coarse_unit) - coarse_range; int i_max = std::ceil((in_z2[idx]-ylo)/coarse_unit) + coarse_range;
		int j_min = std::floor((in_z1[idx]-xlo)/coarse_unit) - coarse_range; int j_max = std::ceil((in_z1[idx]-xlo)/coarse_unit) + coarse_range;
		for (int i = max(0, i_min); i < min(npy, i_max); ++ i) {
			for (int j = max(0, j_min); j < min(npx, j_max); ++ j ) {
				if (std::abs(in_int[idx]-fgx[npy*j+i]) < min_diff) {
					min_diff = std::abs(in_int[idx]-fgx[npy*j+i]);
					out_z1 = recx[j]; out_z2 = recy[i]; out_int = fgx[npy*j+i];
				}
			}
		}
		if (std::abs(in_int[idx]-out_int) > max_fitting_int_diff) {
			max_fitting_int_diff = std::abs(in_int[idx]-out_int);
		}
		Fitted_z1.push_back(out_z1); Fitted_z2.push_back(out_z2); Fitted_intensity.push_back(out_int);
	}
	cout << "The max intensity fitting difference is: " << max_fitting_int_diff << endl;
	//for (int idx = 0; idx < data_size; ++ idx) {
	//	vector<double> tz1; tz1.clear(); vector<double> tz2; tz2.clear(); vector<double> tinten; tinten.clear();
	//	int i_min = std::floor((in_z2[idx]-ylo)/coarse_unit) - coarse_range; int i_max = std::ceil((in_z2[idx]-ylo)/coarse_unit) + coarse_range;
	//	for (int i = max(0, i_min); i < min(npy, i_max); ++ i) {
	//		int j_min = std::floor((in_z1[idx]-xlo)/coarse_unit) - coarse_range; int j_max = std::ceil((in_z1[idx]-xlo)/coarse_unit) + coarse_range;
	//		for (int j = max(0, j_min); j < min(npx, j_max); ++ j ) {
	//			if (std::abs(in_int[idx]-fgx[npy*j+i]) < coarse_diff) {
	//				tz1.push_back(recx[j]); tz2.push_back(recy[i]); tinten.push_back(fgx[npy*j+i]);
	//			}
	//		}
	//	}
	//	if (tz1.size() < 1) {
	//		//cout << "No candidates found for vertex " << idx << " in the coarse level." << endl;
	//		Fitted_z1.push_back(in_z1[idx]); Fitted_z2.push_back(in_z2[idx]); Fitted_intensity.push_back(in_int[idx]); invalid_count ++;
	//	} else {
	//		// find the best match from all candidates
	//		double max_dis = -1000; int max_id = -1;
	//		for (int ii = 0; ii < tz1.size(); ++ ii) {
	//			double temp_dis = std::sqrt(std::powl((in_z1[idx]-tz1[ii]), 2.0) + std::powl((in_z2[idx]-tz2[ii]), 2.0));
	//			if (max_dis < temp_dis) {
	//				max_dis = temp_dis; max_id = ii;
	//			}
	//		}
	//		Fitted_z1.push_back(tz1[max_id]); Fitted_z2.push_back(tz2[max_id]); Fitted_intensity.push_back(tinten[max_id]);
	//	}
	//}
	//cout << invalid_count << " vertices are not fitted in the coarse level. " << endl;
	if (recx)	NAG_FREE(recx);
	if (recy)	NAG_FREE(recy);
	if (fgx)	NAG_FREE(fgx);
	// finish fitting in the coarse level, results are stored in coarse_z1, coarse_z2, coarse_inten vector

	//// start fitting in the fine level, results are store in Fitted_z1, Fitted_z2, Fitted_intensity vector
	//double low_z1 = xlo;  double low_z2, upper_z1, upper_z2; invalid_count = 0;
	//while (low_z1 < xhi) {
	//	upper_z1 = min(low_z1+fine_step, xhi);	
	//	low_z2 = ylo; 
	//	while (low_z2 < yhi) {
	//		upper_z2 = min(low_z2+fine_step, yhi);
	//		npx = std::ceil((upper_z1-low_z1)/fine_unit); npy = std::ceil((upper_z2-low_z2)/fine_unit);
	//		//construct the result fitted data for nag_2d_spline_eval_rect 
	//		if (npx >= 1 && npy >= 1) {
	//			if (!(recx = NAG_ALLOC(npx,double)) || !(recy = NAG_ALLOC(npy,double)) || !(fgx = NAG_ALLOC(npx*npy,double)) ) {
	//				cout<<"Allocation output file memory failure"<<endl;	return false;
	//			}
	//		}
	//		else {
	//			cout<<"Invalid npx or npy"<<endl; return false;
	//		}
	//		for (int i = 0; i < npx; i++) {
	//			recx[i] = low_z1 + fine_unit*i;
	//		}
	//		for (int i = 0; i < npy; i++) {
	//			recy[i] = low_z2 + fine_unit*i;
	//		}//end of construction
	//		//start of nag_2d_spline_eval_rect calculation, but it will cause error when error point at the boundary of this region.
	//		nag_2d_spline_eval_rect(npx, npy, recx, recy, fgx, &splinex, &fail);
	//		if (fail.code != NE_NOERROR) {
	//			cout<<"Error from nag_2d_spline_eval_rect in cx. "<<fail.message<<endl;		return false;
	//		}

	//		for (int idx = 0; idx < data_size; ++ idx) {
	//			// need to judge whether the fitted result in coarse level is in this fine level range
	//			if (Fitted_z1[idx] < low_z1 || Fitted_z1[idx] > upper_z1 || Fitted_z2[idx] < low_z2 || Fitted_z2[idx] > upper_z2) {
	//				continue;
	//			}

	//			vector<double> tz1; tz1.clear(); vector<double> tz2; tz2.clear(); vector<double> tinten; tinten.clear();
	//			int i_min = std::floor((Fitted_z2[idx]-low_z2)/fine_unit) - fine_range; int i_max = std::ceil((Fitted_z2[idx]-low_z2)/fine_unit) + fine_range;
	//			for (int i = max(0, i_min); i < min(npy, i_max); ++ i) {
	//				int j_min = std::floor((Fitted_z1[idx]-low_z1)/fine_unit) - fine_range; int j_max = std::ceil((Fitted_z1[idx]-low_z1)/fine_unit) + fine_range;
	//				for (int j = max(0, j_min); j < min(npx, j_max); ++ j ) {
	//					if (std::abs(Fitted_intensity[idx]-fgx[npy*j+i]) < fine_diff) {
	//						tz1.push_back(recx[j]); tz2.push_back(recy[i]); tinten.push_back(fgx[npy*j+i]);
	//					}
	//				}
	//			}
	//			if (tz1.size() < 1) {
	//				//cout << "No candidates found for vertex " << idx << " in the fine level." << endl; // the result in Fitted_z1 remain the previous result
	//				invalid_count ++; 
	//			} else {
	//				// find the best match from all candidates
	//				double max_dis = -1000; int max_id = -1;
	//				for (int ii = 0; ii < tz1.size(); ++ ii) {
	//					double temp_dis = std::sqrt(std::powl((Fitted_z1[idx]-tz1[ii]), 2.0) + std::powl((Fitted_z2[idx]-tz2[ii]), 2.0));
	//					if (max_dis < temp_dis) {
	//						max_dis = temp_dis; max_id = ii;
	//					}
	//				}
	//				Fitted_z1[idx] = tz1[max_id]; Fitted_z2[idx] = tz2[max_id]; Fitted_intensity[idx] = tinten[max_id];
	//			}
	//		}
	//		if (recx)	NAG_FREE(recx);
	//		if (recy)	NAG_FREE(recy);
	//		if (fgx)	NAG_FREE(fgx);
	//		low_z2 += fine_step;
	//	}
	//	low_z1 += fine_step; 
	//}
	//cout << invalid_count << " vertices are not fitted in the fine level. " << endl;
	//if (recx)	NAG_FREE(recx);
	//if (recy)	NAG_FREE(recy);
	//if (fgx)	NAG_FREE(fgx);
	// store the final result color file and ps normal file.
	printf("\nPS normal fitting Time = %lfs ",((double)cv::getTickCount()-timer_start)/cv::getTickFrequency());
	
	// Save the vertex color in image;
	fstream fout; fout.open(vertexintensityfilename, ios::out);
	if (!fout) {
		cout << "Can not open " << vertexintensityfilename << " to save data..." << endl;
	}
	fout << data_size << endl;
	for (int i = 0; i < data_size; i ++) {
		fout << Fitted_intensity[i] << "  " << in_int[i] << endl;
	}
	fout << endl;
	fout.close();

	double x, y, z, xf, yf, zf;
	// Save the fitted vertex color from spline fitting;
	string fittedintensityfile = (MOptions.DirName+"FittingDataComparison"+string(buffer)+".txt");
	fout.open(fittedintensityfile.c_str(), ios::out);
	if (!fout) {
		cout << "Can not open " << fittedintensityfile << " to save data..." << endl;
	}
	for (int i = 0; i < data_size; i ++) {
		to3DNormal(in_z1[i], in_z2[i], &x, &y, &z);
		to3DNormal(Fitted_z1[i], Fitted_z2[i], &xf, &yf, &zf);
		/*fout << in_z1[i] << "   " << Fitted_z1[i] << "   " << abs(in_z1[i] - Fitted_z1[i]) << "  "
			 <<	in_z2[i] << "   " << Fitted_z2[i] << "   " << abs(in_z2[i] - Fitted_z2[i]) << "  "
			 <<	in_int[i]<< "   " << Fitted_intensity[i] << "   " << abs(in_int[i] - Fitted_intensity[i]) << endl;*/
		fout << x << "   " << xf << "   " << abs(x - xf) << "  "
			<<	y << "   " << yf << "   " << abs(y - yf) << "  "
			<<	z << "   " << zf << "   " << abs(z - zf) << "  "
			<<	in_int[i]<< "   " << Fitted_intensity[i] << "   " << abs(in_int[i] - Fitted_intensity[i]) << endl;
	}
	fout << endl;
	fout.close();

	// Save the output ps normal for vertex;
	fout.open(psnormalfilename, ios::out);
	if (!fout) {
		cout << "Can not open " << psnormalfilename << " to save data..." << endl; 
	}
	fout << data_size << endl; 
	for (int i = 0; i < data_size; i ++) {
		to3DNormal(Fitted_z1[i], Fitted_z2[i], &x, &y, &z);
		fout << x << "   " << y << "   " << z << endl;
	}
	fout << endl;
	fout.close();
	cout << "Done ...." << endl;

	if (in_z1) 			NAG_FREE(in_z1);
	if (in_z2) 			NAG_FREE(in_z2);
	if (in_int)			NAG_FREE(in_int);
	if (weights)		NAG_FREE(weights);
	return true;
}

bool CubicSplineFittingPSNormal(int iter_step = 0, int choice = 1, double range_value = 0.5)
{
	char tbuffer[255];
	sprintf(tbuffer, "_%.2d", iter_step);
	string intensityfile = (MOptions.DirName+"VertexColor"+string(tbuffer)+".txt");
	string psnormalfile = (MOptions.DirName+"VertexPSNormal"+string(tbuffer)+".txt");
	return CubicSplineFittingPSNormal(intensityfile.c_str(), psnormalfile.c_str(), iter_step, choice, range_value);
}

void AddAllVisibleColor(int iter_step)
{
	std::cout<<std::endl<<"Start the visible color process of step "<<iter_step << std::endl;
	double ox[2], mx[2];int r, g, b, gray; 
	double view_angle, input_angle, output_angle;
	total_angle_num = (180-0)/theta_step + 1;
#ifdef BRDFFITTING
	total_angle_num = total_angle_num * 3;
#endif
	char temp_char[10];
	sprintf(temp_char, "%.3d", iter_step);
	string TempFolderName = MOptions.DirName + "TempResult\\";
	string est_light_filename = TempFolderName + "LambertainEstimation"+string(temp_char)+".txt";
	std::fstream estlight_out(est_light_filename.c_str(), std::ios::out);
	if (!estlight_out) {
		std::cout << "Can not open " << est_light_filename << " to save data..." << std::endl;
	}
	for (int idx = 0; idx < tracks.size(); idx ++) {
		tracks[idx].vec_visible_color = new double[total_angle_num]; 
		for (int ii = 0; ii < total_angle_num; ii ++) tracks[idx].vec_visible_color[ii] = 0.0;

		for (int i = 0; i < tracks[idx].views.size(); i ++) {
			// first calculate the projected coordinates
			matrixProject(tracks[idx].X,listViewPoints[tracks[idx].views[i]]->P,ox,1);
			if (abs(matrixLength(tracks[idx].Norm, 3) - 1.0) > 0.00001) {
				matrixNorm(tracks[idx].Norm, 3);
			}
			//calculate the angle between light direction and view direction
			view_angle = calAngle(listViewPoints[tracks[idx].views[i]]->LightDirection, listViewPoints[tracks[idx].views[i]]->viewingVector); 
			input_angle = calAngle(listViewPoints[tracks[idx].views[i]]->LightDirection, tracks[idx].Norm);
			//input_angle = (input_angle>90)?(180-input_angle):input_angle;
			if (view_angle > _max_angle) {
				_max_angle = view_angle;
			}
			int ii = (int)(input_angle / theta_step);

#ifdef BRDFFITTING
			output_angle = calAngle(listViewPoints[tracks[idx].views[i]]->viewingVector, tracks[idx].Norm);
			b = CV_IMAGE_ELEM(listViewPoints[tracks[idx].views[i]]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+0);
			g = CV_IMAGE_ELEM(listViewPoints[tracks[idx].views[i]]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+1);
			r = CV_IMAGE_ELEM(listViewPoints[tracks[idx].views[i]]->m_Image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5)*3+2);
			tracks[idx].addViewColor(tracks[idx].views[i],r,g,b, view_angle, input_angle, output_angle);
			tracks[idx].vec_visible_color[ii*3+0] = r;
			tracks[idx].vec_visible_color[ii*3+1] = g;
			tracks[idx].vec_visible_color[ii*3+2] = b;
#else
			gray = CV_IMAGE_ELEM(listViewPoints[tracks[idx].views[i]]->image, uchar,(int)(ox[0]+0.5),(int)(ox[1]+0.5));
			tracks[idx].addViewColor(tracks[idx].views[i], gray, gray, gray, view_angle, input_angle, input_angle);
			tracks[idx].vec_visible_color[ii] = gray;
#endif
			tracks[idx].AddViewLightInfo(tracks[idx].views[i], listViewPoints[tracks[idx].views[i]]->LightDirection, listViewPoints[tracks[idx].views[i]]->LightLuminance);

			// next find the neighbor CameraView of current view that the track point maybe seen, record its color and view information
			for (int j = 0; j < nViews; j ++) {
				if ((theta_lower<thetas[tracks[idx].views[i]*nViews+j] && thetas[tracks[idx].views[i]*nViews+j]<theta_upper) && 
					(distance_lower*medianD < camera_dis[tracks[idx].views[i]*nViews+j] && camera_dis[tracks[idx].views[i]*nViews+j] < distance_upper*medianD)) {
						matrixProject(tracks[idx].X,listViewPoints[j]->P,mx,1);
						//calculate the angle between light direction and view direction
						if (abs(matrixLength(tracks[idx].Norm, 3) - 1.0) > 0.00001) {
							matrixNorm(tracks[idx].Norm, 3);
						}
						view_angle = calAngle(listViewPoints[j]->LightDirection, listViewPoints[j]->viewingVector); 
						input_angle = calAngle(listViewPoints[j]->LightDirection, tracks[idx].Norm);
						//input_angle = (input_angle>90)?(180-input_angle):input_angle;
						if (view_angle > _max_angle) {
							_max_angle = view_angle;
						}
						int ii = (int)(input_angle / theta_step);
#ifdef BRDFFITTING
						output_angle = calAngle(tracks[idx].Norm, listViewPoints[j]->viewingVector);
						b = CV_IMAGE_ELEM(listViewPoints[j]->m_Image, uchar,(int)(mx[0]+0.5),(int)(mx[1]+0.5)*3+0);
						g = CV_IMAGE_ELEM(listViewPoints[j]->m_Image, uchar,(int)(mx[0]+0.5),(int)(mx[1]+0.5)*3+1);
						r = CV_IMAGE_ELEM(listViewPoints[j]->m_Image, uchar,(int)(mx[0]+0.5),(int)(mx[1]+0.5)*3+2);
						tracks[idx].addViewColor(j, r, g, b, view_angle, input_angle, output_angle);
						tracks[idx].vec_visible_color[ii*3+0] = r;
						tracks[idx].vec_visible_color[ii*3+1] = g;
						tracks[idx].vec_visible_color[ii*3+2] = b;
#else
						gray = CV_IMAGE_ELEM(listViewPoints[j]->image, uchar,(int)(mx[0]+0.5),(int)(mx[1]+0.5));
						tracks[idx].addViewColor(j, gray, gray, gray, view_angle, input_angle, input_angle);
						tracks[idx].vec_visible_color[ii] = gray;
#endif
						tracks[idx].AddViewLightInfo(j, listViewPoints[j]->LightDirection, listViewPoints[j]->LightLuminance);
				}
			}
		}
#ifdef BRDFFITTING
		sort(tracks[idx].vec_visible_info.begin(),tracks[idx].vec_visible_info.end(), VIFLess_ViewAngle);
#else
		sort(tracks[idx].vec_visible_info.begin(),tracks[idx].vec_visible_info.end(), VIFLess_InputAngle);
		//average the lambertain parameters for this track point
		tracks[idx].SurfaceLambda = 0.0; int temp_count = 0;
		for (int j = 0; j < tracks[idx].vec_visible_info.size(); j ++) {
			if (tracks[idx].vec_visible_info[j].Lambda > 0) {
				tracks[idx].SurfaceLambda += tracks[idx].vec_visible_info[j].Lambda;
				temp_count ++ ;
			}
		}
		tracks[idx].SurfaceLambda = tracks[idx].SurfaceLambda/temp_count;
		estlight_out << "The average lambertain lambda for track " << idx << " is " << tracks[idx].SurfaceLambda << std::endl;
#endif
	}
	estlight_out.close();
	cout << "The max angle between light direction and view direction is: " << _max_angle << endl;
	string vis_color_filename = TempFolderName + "VisibleColor"+string(temp_char)+".txt";
	std::fstream viscolor_out(vis_color_filename.c_str(), std::ios::out);
	if (!viscolor_out) {
		std::cout << "Can not open " << vis_color_filename << " to save data..." << std::endl;
	}
	for (int idx = 0; idx < tracks.size(); idx ++) {
		for (int i = 0; i < total_angle_num; i ++) {
			viscolor_out << tracks[idx].vec_visible_color[i] << " ";
		}
		viscolor_out << std::endl;
	}
	viscolor_out.close();
	std::cout<<std::endl<<"...done..."<<std::endl;
}

double EstimateMaxColorDistance()
{
	// estimate the max color vector distance between 
	boost::uniform_int<> distribution(0, tracks.size()-1) ;
	boost::mt19937 engine(time(NULL)) ;
	boost::variate_generator<boost::mt19937, boost::uniform_int<> > myrandom (engine, distribution);

	int SampleNumber = 500;_max_dist = -1;
	for (int i = 1; i < SampleNumber; i ++) {
		int t1 = myrandom(); int t2 = myrandom();
		double clr_dis = matrixDist(tracks[t1].vec_visible_color, tracks[t2].vec_visible_color,total_angle_num);
		if (clr_dis > _max_dist) {
			_max_dist = clr_dis;
		}
	}
	cout << "The estimated max color distance is " << _max_dist << endl;
	return _max_dist;
}

bool CompareTwoTrack(const track& t_a, const track& t_b)
{
	// criteria: captured color difference in different view_angle.
	double clr_dis = matrixDist(t_a.vec_visible_color, t_b.vec_visible_color, total_angle_num);
	if (clr_dis < _max_dist) {
		return true;
	}
	else
		return false;// if the two tracks belong to same partition, have same material; else, return false
}

/*
void PartitionTracks(int iter_step)
{
	std::cout<<"Partition tracks of step "<<iter_step<<std::endl;
	// partition tracks into parts with different materials and lambertain lambda
	vector <int> Part_label;
	//EstimateMaxColorDistance();
	//cv::partition(tracks, Part_label, CompareTwoTrack); // using cv::partition

	bool par_flag = true; int repeat_count = 0;
	cv::Mat kmeans_label(tracks.size(), 1, CV_8UC1);
	cv::TermCriteria termcrit(CV_TERMCRIT_ITER, 300, 1.0); termcrit.maxCount = 300; termcrit.epsilon = 1.0; termcrit.type = 1;
	cv::Mat* tracks_matp = NULL;
	std::cout << "Current K number is: ";
	while (par_flag && K_Num > 1) {
		//tracks_matp = new cv::Mat(tracks.size(),total_angle_num,CV_32F); // using all visible color information
		tracks_matp = new cv::Mat(tracks.size(),1,CV_64F); // using only lambertain lambda
		for (int i = 0; i < tracks.size(); i ++) {
			//for (int j = 0; j < total_angle_num; j ++) {
			//	tracks_matp->at<double>(i,j) = tracks[i].vec_visible_color[j];// using all visible color information
			//}
			tracks_matp->at<double>(i,0) = tracks[i].SurfaceLambda; // using only lambertain lambda
		}
		tracks_matp->convertTo(*tracks_matp, CV_32F);
		cv::kmeans(*tracks_matp, K_Num, kmeans_label, termcrit, 50, cv::KMEANS_PP_CENTERS);
		//double **Averge_Center = kutility::allocate<double>(K_Num,total_angle_num); // using all visible color information
		double *Average_Lambda = kutility::allocate<double>(K_Num); // using only lambertain lambda
		double *Stat_Count = kutility::allocate<double>(K_Num);
		for (int i = 0; i < K_Num; i ++) {
			Stat_Count[i] = 0;
			//for (int j = 0; j < total_angle_num; j ++) {
			//	Averge_Center[i][j] = 0.0;
			//}
			Average_Lambda[i] = 0.0;
		}
		for (int i = 0; i < tracks.size(); i ++) {
			//std::cout << kmeans_label.at<int>(i,1) << std::endl;
			//matrixAdd(tracks[i].vec_visible_color, Averge_Center[kmeans_label.at<int>(i,0)], total_angle_num);
			Average_Lambda[kmeans_label.at<int>(i,0)] += tracks[i].SurfaceLambda;
			Stat_Count[kmeans_label.at<int>(i,0)] = Stat_Count[kmeans_label.at<int>(i,0)] + 1;
		}
		for (int i = 0; i < K_Num; i ++) {
			//matrixDiv(Averge_Center[i], total_angle_num, Stat_Count[i]);
			Average_Lambda[i] = Average_Lambda[i]/Stat_Count[i];
		}
		int cur_K = K_Num;
		for (int i = 1; i < K_Num; i ++) {
			//if (matrixDist(Averge_Center[i],Averge_Center[0], total_angle_num) < _max_dist/4) {
			//	cur_K --;
			//}
			if (abs(Average_Lambda[i]-Average_Lambda[0]) < Average_Lambda[i] * 0.1) {
				cur_K --;
			}
		}
		if (cur_K == K_Num) {
			repeat_count ++;
			if (repeat_count > 2) {
				par_flag = false;
			}
		} else {
			repeat_count = 0;
			par_flag = true;
		}
		//kutility::deallocate(Averge_Center, K_Num);
		kutility::deallocate(Average_Lambda);
		kutility::deallocate(Stat_Count);
		K_Num = cur_K;
		std::cout << K_Num << " ";
	}
	//tracks_matp = new cv::Mat(tracks.size(),total_angle_num,CV_32F); // using all visible color information
	tracks_matp = new cv::Mat(tracks.size(),1,CV_64F); // using only lambertain lambda
	for (int i = 0; i < tracks.size(); i ++) {
		//for (int j = 0; j < total_angle_num; j ++) {
		//	tracks_matp->at<double>(i,j) = tracks[i].vec_visible_color[j];// using all visible color information
		//}
		tracks_matp->at<double>(i,0) = tracks[i].SurfaceLambda; // using only lambertain lambda
	}
	tracks_matp->convertTo(*tracks_matp, CV_32F);
	cv::kmeans(*tracks_matp, K_Num, kmeans_label, termcrit, 50, cv::KMEANS_PP_CENTERS);
 	std::cout << "Finally the object will be clustered into " << K_Num << " parts according to reflectance properties" << std::endl;
	//average the lambertain parameters for each kinds of material
	double *Average_Lambda = kutility::allocate<double>(K_Num);
	double *Stat_Count = kutility::allocate<double>(K_Num);
	for (int i = 0; i < K_Num; i ++) {
		Stat_Count[i] = 0;	Average_Lambda[i] = 0.0;
	}
	for (int idx = 0 ; idx < tracks.size(); idx ++) {
		//std::cout << kmeans_label.at<int>(idx,0) << " " << kmeans_label.at<int>(idx,1) << " ";
		Average_Lambda[(int)(kmeans_label.at<int>(idx,0))] += tracks[idx].SurfaceLambda;
		Stat_Count[(int)(kmeans_label.at<int>(idx,0))] +=1;
	}
	for (int i = 0; i < K_Num; i++) {
		Average_Lambda[i] = Average_Lambda[i]/Stat_Count[i];
		std::cout << "The " << i << "-th fitted Lambertain parameters is " << Average_Lambda[i] << " , for " << Stat_Count[i] <<" track points."<< std::endl;
	}
	char temp_char[10];
	sprintf(temp_char, "%.3d", iter_step);
	string TempFolderName = MOptions.DirName + "TempResult\\";
	string est_light_filename = TempFolderName + "LambertainUpdate"+string(temp_char)+".txt";
	std::fstream estlight_out(est_light_filename.c_str(), std::ios::out);
	if (!estlight_out) {
		std::cout << "Can not open " << est_light_filename << " to save data..." << std::endl;
	}
	for (int kk = 0; kk < nViews; kk ++) {
		listViewPoints[kk]->visiblematerialcount.clear();
		listViewPoints[kk]->visiblematerialcount.resize(K_Num);
		for (int kkk = 0; kkk < K_Num; kkk ++) {
			listViewPoints[kk]->visiblematerialcount[kkk] = 0;
		}
	}
	for (int idx = 0 ; idx < tracks.size(); idx ++) {
		tracks[idx].material_label = (int)(kmeans_label.at<int>(idx,0));
		for (int i = 0; i < tracks[idx].views.size(); i ++) {
			listViewPoints[tracks[idx].views[i]]->visiblematerialcount[tracks[idx].material_label] ++;
		}
		tracks[idx].Fitted_SurfaceLambda = Average_Lambda[tracks[idx].material_label];
		for (int i = 0; i < tracks[idx].vec_visible_info.size(); i ++) {
			tracks[idx].vec_visible_info[i].Lambda = tracks[idx].Fitted_SurfaceLambda;
		}
		//std::cout << "Track " << idx << ": Label is " << tracks[idx].material_label << ", Lambda is: " << tracks[idx].SurfaceLambda << ", new Lambda is: " << tracks[idx].Fitted_SurfaceLambda << std::endl;
		estlight_out<< "Track " << idx << ": Label is " << tracks[idx].material_label << ", Lambda is: " << tracks[idx].SurfaceLambda << ", new Lambda is: " << tracks[idx].Fitted_SurfaceLambda << std::endl;
	}
	estlight_out.close();
	kutility::deallocate(Stat_Count);
	kutility::deallocate(Average_Lambda);
}
*/

void CalculatePhotometricNormal(int iter_step)
{
	std::cout<<"Start the photometric norm update process of step "<<iter_step<<std::endl;
	//Update the track norm according to fitted lambertain lambda, visible color and light direction and luminance
	cv::Mat PNorm = cv::Mat(3,1,CV_64F); cv::Mat SumLLT; cv::Mat SumLC;
	cv::Mat LightD = cv::Mat(3,3,CV_64F);
	cv::Mat Vcolor = cv::Mat(3,1,CV_64F);
	//first update the cos(theta) according to the fitted lambda
	// save the estimated light information of all the CameraView
	char temp_char[10];
	sprintf(temp_char, "%.3d", iter_step);
	string TempFolderName = MOptions.DirName + "TempResult\\";
	string est_light_filename = TempFolderName + "PhotometricNorm"+string(temp_char)+".txt";
	std::fstream estlight_out(est_light_filename.c_str(), std::ios::out);
	if (!estlight_out) {
		std::cout << "Can not open " << est_light_filename << " to save data..." << std::endl;
	}
	for (int idx = 0; idx < tracks.size(); idx ++) {
		SumLLT.zeros(3,3,CV_64F); SumLC.zeros(3,1,CV_64F);
		if (tracks[idx].vec_visible_info.size() < 3) {
			continue;
		}
		for (int i = 0; i < SumLLT.rows; i ++) {
			for (int j = 0; j < SumLLT.cols; j ++) {
				SumLLT.at<double>(i,j) = 0.0;
			}
			SumLC.at<double>(i,0) = 0.0;
		}
		for (int i = 0; i < tracks[idx].vec_visible_info.size(); i ++) {
			cv::Mat LightD = cv::Mat(3,1,CV_64F, tracks[idx].vec_visible_info[i].LightDirection); 
			LightD = LightD * tracks[idx].vec_visible_info[i].LightLuminance[0]; //Mat_Print("light direction", LightD);
			//Mat_Print("Prev sumllt", SumLLT); Mat_Print("temp llt",LightD * LightD.t() * tracks[idx].Fitted_SurfaceLambda);
			SumLLT = SumLLT + LightD * LightD.t() * tracks[idx].Fitted_SurfaceLambda; //Mat_Print("Current sumllt", SumLLT);
			//Mat_Print("Prev sumlc", SumLC); Mat_Print("Temp lc",LightD*tracks[idx].vec_visible_info[i].view_color[0]);
			SumLC = SumLC + LightD*tracks[idx].vec_visible_info[i].view_color[0]; //Mat_Print("Current sumlc", SumLC);
		}
		PNorm = SumLLT.inv() * SumLC; //Mat_Print("Estimated photometric norm", PNorm); 
		double temp_pn[3]; temp_pn[0] = PNorm.at<double>(0,0); temp_pn[1] = PNorm.at<double>(1,0); temp_pn[2] = PNorm.at<double>(2,0);
		double new_lambda = matrixLength(temp_pn, 3);		matrixNorm(temp_pn, 3);
		//std::cout<<"The lambda of track "<<idx<<" changes from "<<tracks[idx].SurfaceLambda<<" to "<<new_lambda<<std::endl;
		estlight_out<<"The lambda of track "<<idx<<" changes from "<<tracks[idx].SurfaceLambda<<" to "<<new_lambda<<std::endl;
		//std::cout << "Norm of track " << idx << " changes from (" << tracks[idx].Norm[0] << ", " << tracks[idx].Norm[1] << ", " << tracks[idx].Norm[2] ;
		estlight_out << "Norm of track " << idx << " changes from (" << tracks[idx].Norm[0] << ", " << tracks[idx].Norm[1] << ", " << tracks[idx].Norm[2] ;
		matrixCopy(temp_pn, tracks[idx].PhotometricNorm, 3);
		//std::cout << ") to (" << tracks[idx].PhotometricNorm[0]  << ", " << tracks[idx].PhotometricNorm[1]  << ", " << tracks[idx].PhotometricNorm[2]  << ")." << std::endl;
		estlight_out << ") to (" << tracks[idx].PhotometricNorm[0]  << ", " << tracks[idx].PhotometricNorm[1]  << ", " << tracks[idx].PhotometricNorm[2]  << ")." << std::endl;
	}
	estlight_out.close();
	std::cout<<"...done..."<<std::endl;
	return;
}

void ReleaseResources()
{
	if (paired != NULL) {
		for (int i=0;i<nViews;i++)
		{
			delete[] paired[i];
		}
		delete[] paired;
	}
	if (camera_dis != NULL) {
		delete[] camera_dis;
	}
	if (thetas != NULL) {
		delete[] thetas;
	}
	for (int i=0;i<nViews;i++) 
	{
		delete listViewPoints[i];
		//listViewPoints[i] = NULL;	
	}
	delete[] listViewPoints;
	tracks.clear();
}