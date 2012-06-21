extern int cntF;
extern double *** Rij;
extern double *** Tij;
extern vector<track> tracks;
/*
double crossCorrelationInHomography(int k)
{
	double *x = tracks[k].X;
	double xi[2], xj[2];
	double _n[3],n[3];

	int i,j,ii,jj;
	int viewI, viewJ;
	double Xicam[3];

	double *Kj;
	double d;
	double H[9],H1[9];

	double fDaisy = 0;

	float* thor0 = new float[daisyNormalsLength];
	float* thor1 = new float[daisyNormalsLength];

	for (i=0;i<tracks[k].size;i++)
	{
	
		viewI = tracks[k].views[i];
		
		to3DNormal(x[3],x[4],&_n[0],&_n[1],&_n[2]);
		matrixMul(listViewPoints[viewI]->R,3,3,_n,3,1,n);
		
		matrixMul(listViewPoints[viewI]->R,3,3,x,3,1,Xicam);
		matrixAdd(listViewPoints[viewI]->T,Xicam,3);

		// image: X-row, Y-col  ->  2D coordinate y first; then x
		//xi[0] = tracks[k].py[i];
		//xi[1] = tracks[k].px[i];
		matrixProject(x,listViewPoints[viewI]->P,xi,0);
		//xi[0] = int(xi[0]);
		//xi[1] = int(xi[1]);
	
		//thor0 = daisyMEM[viewI][tracks[k].px[i]][tracks[k].py[i]];

		for (int dx=-1;dx<=1;dx++)
			for (int dy=-1;dy<=1; dy++)
			{
				xi[0]+=dx;
				xi[1]+=dy;
				if (matrixBound(xi[0],0,listViewPoints[viewI]->width) && 
					matrixBound(xi[1],0,listViewPoints[viewI]->height))
				{
					for (int z=0;z<daisyNormalsLength;z++) thor0[z] = 0;
					//listViewPoints[viewI]->desc_original->get_descriptor(tracks[k].px[i],tracks[k].py[i],0,thor0);
					listViewPoints[viewI]->desc_original->get_descriptor(xi[1],xi[0],0,thor0);

					d = dotProduct(n,Xicam,3);

					for (j=0;j<tracks[k].size;j++)
					{
						viewJ = tracks[k].views[j];
						if (viewI != viewJ)
						{
							Kj = listViewPoints[viewJ]->K;
							for (ii=0;ii<3;ii++)
								for(jj=0;jj<3;jj++) H[ii*3+jj] = Rij[viewI][viewJ][ii*3+jj] + Tij[viewI][viewJ][ii] * n[jj] / d;
							matrixMul(Kj,3,3,H,3,3,H1);
							matrixMul(H1,3,3,listViewPoints[viewI]->Kinv,3,3,H);
							matrixHomographyProject(xi,H,xj);

							if (matrixBound(xj[1],0,listViewPoints[viewJ]->height) &&
								matrixBound(xj[0],0,listViewPoints[viewJ]->width) )
							{

								for (int z=0;z<daisyNormalsLength;z++) thor1[z] = 0;
								listViewPoints[viewJ]->desc_original->get_descriptor(xj[1],xj[0],0,thor1);					
								fDaisy += cross_correlation(thor0, thor1, daisyNormalsLength);

							} else fDaisy ++;
						}
					}
				} else fDaisy++;
				xi[0]-=dx;
				xi[1]-=dy;
			}
	}
	deallocate(thor0);
	deallocate(thor1);
	
	fDaisy/=tracks[k].size*(tracks[k].size-1)*9;
	return fDaisy;
}
*/


double getColorValue(uchar *im, int h, int w, double x, double y)
{
	return im[((int)(x+0.5))*w+(int)(y+0.5)];

	// bilinear interpolation
	if (x>=h-2 || y>=w-2) return im[((int)(x+0.5))*w+(int)(y+0.5)];
	int x0 = (int)x;
	int x1 = x0+1;
	int y0 = (int) y;
	int y1 = y0+1;
	double f00 = im[x0*w+y0];
	double f01 = im[x0*w+y1];
	double f10 = im[x1*w+y0];
	double f11 = im[x1*w+y1];
	double dx = x-x0;
	double dy = y-y0;
	double b1 = f00;
	double b2 = f10-f00;
	double b3 = f01-f00;
	double b4 = f00 - f10 - f01 + f11;
	double res = b1+b2*dx+b3*dy+b4*dx*dy;
	return res;
}

double initNormalsF(int k, double *x)
{
	//printf("%d %lf %lf %lf %lf %lf\n", k, x[0], x[1], x[2], x[3], x[4]);
	cntF++;

	double f = 0.0;

	double xi[2], xj[2], XicamOi[3];
	double _n[3],n[3];
	double dx, dy, inbetweenAngle;

	int i,j,ii,jj;
	
	bool valid = true;
	for (i=0;i<tracks[k].size; i++)
	{	
		matrixProject(x,listViewPoints[tracks[k].views[i]]->P,xi,1);
		dx = fabs(xi[0] - tracks[k].px[i]);
		dy = fabs(xi[1] - tracks[k].py[i]);

		if (dx>gamma) 
		{
			f+= dx*1000000000;
			valid = false;
		}
		if (dy>gamma) 
		{
			f+= dy*1000000000;		
			valid = false;
		}

		to3DNormal(x[3],x[4],&n[0],&n[1],&n[2]);
		
		matrixSubtract(listViewPoints[tracks[k].views[i]]->C,x,XicamOi,3);
		matrixNorm(XicamOi,3);

		inbetweenAngle = dotProduct(n,XicamOi,3);
		if (inbetweenAngle < angleThreshold) 
		{
			f += (1.0-inbetweenAngle) * 1000000000;
			valid = false;
		}
	}
	if (!valid) return f;
	
	int viewI, viewJ;
	double Xicam[3];

	double *Kj;
	double d;
	double H[9],H1[9];



	double* thor0 = new double[patchSize*patchSize];
	double* thor1 = new double[patchSize*patchSize];
	
	for (i=0;i<tracks[k].size;i++)
	{
		viewI = tracks[k].views[i];

		to3DNormal(x[3],x[4],&_n[0],&_n[1],&_n[2]);
		matrixMul(listViewPoints[viewI]->R,3,3,_n,3,1,n);

		matrixMul(listViewPoints[viewI]->R,3,3,x,3,1,Xicam);
		matrixAdd(listViewPoints[viewI]->T,Xicam,3);

		// image: X-row, Y-col  ->  2D coordinate y first; then x
		xi[0] = tracks[k].py[i];
		xi[1] = tracks[k].px[i];
		//matrixProject(x,listViewPoints[viewI]->P,xi,0);
		//xi[0] = int(xi[0]+0.5);
		//xi[1] = int(xi[1]+0.5);

		//printf("%d %lf %lf %d %d\n", viewI, xi[1], xi[0], tracks[k].px[i], tracks[k].py[i]);		

		d = dotProduct(n,Xicam,3);
		for (j=0;j<tracks[k].size;j++)
		{
			viewJ = tracks[k].views[j];
			if (viewI != viewJ)
			{
				Kj = listViewPoints[viewJ]->K;
				for (ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++) H[ii*3+jj] = Rij[viewI][viewJ][ii*3+jj] + Tij[viewI][viewJ][ii] * n[jj] / d;
				matrixMul(Kj,3,3,H,3,3,H1);
				matrixMul(H1,3,3,listViewPoints[viewI]->Kinv,3,3,H);

				int index = 0;
				for (int k=0;k<patchSize*patchSize;k++)
				{
					thor0[k] = thor1[k] = 0;
				}
				for (int dx=-patchSize2;dx<=patchSize2 ;dx++)
					for (int dy=-patchSize2;dy<=patchSize2 ; dy++)
					{						
						xi[0]+=patchGridRes*dx;
						xi[1]+=patchGridRes*dy;
						if (matrixBound(xi[0],0,listViewPoints[viewI]->width) && matrixBound(xi[1],0,listViewPoints[viewI]->height))
							//if (listViewPoints[viewI]->sil[(int)xi[1]][(int)xi[0]])
						{
							thor0[index] = getColorValue(listViewPoints[viewI]->ucharXImage, 
								listViewPoints[viewI]->xImage->height, listViewPoints[viewI]->xImage->width,
								xi[1]*scaleXd, xi[0]*scaleXd);
							matrixHomographyProject(xi,H,xj);

							//if (dx==0 && dy==0)
							//printf("%d %lf %lf %d %d\n", viewJ, xj[1], xj[0], tracks[k].px[j], tracks[k].py[j]);
							if (matrixBound(xj[1],0,listViewPoints[viewJ]->height) &&
								matrixBound(xj[0],0,listViewPoints[viewJ]->width) )
								//if (listViewPoints[viewJ]->sil[(int)xj[1]][(int)xj[0]])
							{


								//f+= abs((double)listViewPoints[viewI]->ucharImage[((int)xi[1])*listViewPoints[viewI]->width+(int)xi[0]] - (double)listViewPoints[viewJ]->ucharImage[((int)xj[1])*listViewPoints[viewJ]->width+(int)xj[0]]);
								//thor1[index] = listViewPoints[viewJ]->ucharXImage[((int)(xj[1]*scaleXd))*listViewPoints[viewI]->xImage->width+(int)(xj[0]*scaleXd)];
								thor1[index] = getColorValue(listViewPoints[viewJ]->ucharXImage, 
								listViewPoints[viewJ]->xImage->height, listViewPoints[viewJ]->xImage->width,
								xj[1]*scaleXd, xj[0]*scaleXd);
								//f+=abs(thor0[index] - thor1[index]);
							} else 
							{
								//f+= daisyNormalsLength;				
							}				
						}
						xi[0]-=patchGridRes*dx;
						xi[1]-=patchGridRes*dy;
						index ++;
					}
				//for (int k=0;k<pSizeX*pSizeY;k++) f+=abs(thor0[k]-thor1[k]);
				f-= cross_correlation(thor0,thor1,patchSize*patchSize);


			}
		}

		//printf("\n");
	}

	f = f/(tracks[k].size*(tracks[k].size-1));
	deallocate(thor0);
	deallocate(thor1);
	

	//printf("time = %9d, k = %6d, f = %19.9lf, x0 = %9.9lf, x1 = %9.9lf, x2 = %9.9lf, x3 = %9.9lf, x4 = %9.9lf\n", cntF, k, f, 
	//	x[0], x[1], x[2], x[3], x[4]);

	return f;
}


lbfgsfloatval_t LBFGS_F_function(int k, lbfgsfloatval_t *x)
{
	//printf("%d %lf %lf %lf %lf %lf\n", k, x[0], x[1], x[2], x[3], x[4]);
	cntF++;

	
	//int i;
    //lbfgsfloatval_t f = 0.0;

    //for (i = 0;i < 4;i += 2) {
    //    lbfgsfloatval_t t1 = 1.0 - x[i];
    //    lbfgsfloatval_t t2 = 10.0 * (x[i+1] - x[i] * x[i]);
    //    f += t1 * t1 + t2 * t2;
    //}
	//f+= (1-x[4])*(1-x[4]);
	

	
	lbfgsfloatval_t f = 0.0;

	double xi[2], xj[2], XicamOi[3];
	double _n[3],n[3];
	double dx, dy, inbetweenAngle;

	int i,j,ii,jj;
	
	for (i=0;i<tracks[k].size; i++)
	{	
		matrixProject(x,listViewPoints[tracks[k].views[i]]->P,xi,1);
		dx = fabs(xi[0] - tracks[k].px[i]);
		dy = fabs(xi[1] - tracks[k].py[i]);

		if (dx>gamma) f+= dx*1e9;
		if (dy>gamma) f+= dy*1e9;		

		to3DNormal(x[3],x[4],&n[0],&n[1],&n[2]);
		
		matrixSubtract(listViewPoints[tracks[k].views[i]]->C,x,XicamOi,3);
		matrixNorm(XicamOi,3);

		inbetweenAngle = dotProduct(n,XicamOi,3);
		if (inbetweenAngle < angleThreshold) f += (1.0-inbetweenAngle) * 1e9;
	}
	
	int viewI, viewJ;
	double Xicam[3];

	double *Kj;
	double d;
	double H[9],H1[9];

	//double fDaisy = 0;

	float* thor0 = new float[daisyNormalsLength];
	float* thor1 = new float[daisyNormalsLength];

	for (i=0;i<tracks[k].size;i++)
	{
	
		viewI = tracks[k].views[i];
		
		to3DNormal(x[3],x[4],&_n[0],&_n[1],&_n[2]);
		matrixMul(listViewPoints[viewI]->R,3,3,_n,3,1,n);
		
		matrixMul(listViewPoints[viewI]->R,3,3,x,3,1,Xicam);
		matrixAdd(listViewPoints[viewI]->T,Xicam,3);

		// image: X-row, Y-col  ->  2D coordinate y first; then x
		xi[0] = tracks[k].py[i];
		xi[1] = tracks[k].px[i];
		//matrixProject(x,listViewPoints[viewI]->P,xi,0);
		//xi[0] = int(xi[0]+.5);
		//xi[1] = int(xi[1]+.5);
	
		//thor0 = daisyMEM[viewI][tracks[k].px[i]][tracks[k].py[i]];
		
		
		listViewPoints[viewI]->desc_original->get_unnormalized_descriptor(xi[1],xi[0],0,thor0);
		
		//printf("%d %lf %lf %d %d\n", viewI, xi[1], xi[0], tracks[k].px[i], tracks[k].py[i]);
		
		//for (int t=0; t<daisyNormalsLength; t++)
		//{
		//	printf("%10.9lf ", thor0[t]);
		//	if ((t+1)%daisyParasNormals[3]==0) printf("\n");
		//}
		

		d = dotProduct(n,Xicam,3);

		for (j=0;j<tracks[k].size;j++)
		{
			viewJ = tracks[k].views[j];
			if (viewI != viewJ)
			{
				Kj = listViewPoints[viewJ]->K;
				for (ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++) H[ii*3+jj] = Rij[viewI][viewJ][ii*3+jj] + Tij[viewI][viewJ][ii] * n[jj] / d;
				matrixMul(Kj,3,3,H,3,3,H1);
				matrixMul(H1,3,3,listViewPoints[viewI]->Kinv,3,3,H);
				matrixHomographyProject(xi,H,xj);

				//printf("%d %lf %lf %d %d\n", viewJ, xj[1], xj[0], tracks[k].px[j], tracks[k].py[j]);
				if (matrixBound(xj[1],0,listViewPoints[viewJ]->height) &&
					matrixBound(xj[0],0,listViewPoints[viewJ]->width) 
					//&& fabs(xj[1]-tracks[k].px[j])<=gamma &&
					//fabs(xj[0]-tracks[k].py[j])<=gamma
					)
				{
					
					listViewPoints[viewJ]->desc_original->get_unnormalized_descriptor(xj[1],xj[0],0,thor1);										
					
					//for (int t=0; t<daisyNormalsLength; t++)
					//{
					//	printf("%10.9lf ", thor1[t]);
					//	if ((t+1)%daisyParasNormals[3]==0) printf("\n");
					//}
					
					//double tmp = matrixL1(thor0, thor1, daisyNormalsLength);
					//f+= tmp;
					f-= cross_correlation(thor0, thor1, daisyNormalsLength);

				} else 
				{
					//printf("out of bound\n");
					f+= daisyNormalsLength;				
					//fDaisy --;
				}				
			}
		}
		//printf("\n");
	}

	f = f/(tracks[k].size*(tracks[k].size-1));

	deallocate(thor0);
	deallocate(thor1);
	

	//printf("time = %9d, k = %6d, f = %19.9lf, x0 = %9.9lf, x1 = %9.9lf, x2 = %9.9lf, x3 = %9.9lf, x4 = %9.9lf\n", cntF, k, f, 
	//	x[0], x[1], x[2], x[3], x[4]);

	//fDaisy/=tracks[k].size*(tracks[k].size-1);
	//if (fDaisy<0.7) f+= tracks[k].size*(tracks[k].size-1)*daisyLength;

	return f;
}




lbfgsfloatval_t evaluate(
    void *instance,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{

	int k = (int)instance;	
	//lbfgsfloatval_t f = LBFGS_F_function(k,x);
	lbfgsfloatval_t f = initNormalsF(k,x);
	for (int i=0;i<n;i++)
	{
		x[i] += deltaX;
		//lbfgsfloatval_t fnew = LBFGS_F_function(k,x);
		lbfgsfloatval_t fnew = initNormalsF(k,x);
		g[i] = (fnew-f)/deltaX;
		x[i] -= deltaX;

	}
	//printf("\n");
	return f;

}

int lbfgsFunction(int k)
{	
	if (k%10000==0)
	printf("Optimizing point %d / %d\n",k,tracks.size());
	//_getch();
	//printf("%d / %d \n",k , tracks.size());
	lbfgsfloatval_t fx;
	lbfgsfloatval_t *x = tracks[k].X;
	lbfgs_parameter_t param;

	fx = 1000000000;	
	double n[3];
	
	int viewId = tracks[k].views[0];		
	matrixSubtract(listViewPoints[viewId]->C,x,n,3);
	matrixNorm(n,3);
	toSpherical(n[0],n[1],n[2],&x[3],&x[4]);	
	
	//printf("\n");
	

	int dx_ = 0;
	int dy_ = 0;	
	for (int dx = -angleSearchSpace; dx <= angleSearchSpace; dx+=angleSearchStep)
	{
		x[3] += dx * RADIAN;
		for (int dy = -angleSearchSpace; dy<= angleSearchSpace; dy+=angleSearchStep)
		{			
			x[4] += dy * RADIAN;
			double ff = initNormalsF(k,x);
			//double ff = LBFGS_F_function(k,x);
			
			if (ff<fx)
			{
				fx = ff;
				dx_ = dx;
				dy_ = dy;
			}			
			x[4] -= dy * RADIAN;
		}
		x[3] -= dx * RADIAN;
	}
	x[3] += dx_ * RADIAN;
	x[4] += dy_ * RADIAN;
	

	//printf("%d %lf %lf %lf %lf %lf fx = %lf\n", k, x[0], x[1], x[2], x[3], x[4], fx);

	//tracks[k].fx = fx;
	//return 0;

	//printf("\n");
	
	//if (fx>=tracks[k].size*(tracks[k].size-1)*daisyLength)  tracks[k].valid = false;
	

	/* Initialize the parameters for the L-BFGS optimization. */
	lbfgs_parameter_init(&param);
	//param.orthantwise_c = 1;
	//param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;

	/*
	Start the L-BFGS optimization; this will invoke the callback functions
	evaluate() and progress() when necessary.
	*/
	int ret = 0;
	if (fx<=0)
	ret = lbfgs(5, x, &fx, evaluate, NULL, (void *)k, &param);
	tracks[k].fx = fx;

	/* Report the result. */
	//printf("L-BFGS optimization terminated with status code = %d\n", ret);
	//printf("%d:     fx = %f, code = %d\n", k, fx, ret);

	//lbfgs_free(x);

	//for (int i=0;i<5;i++) tracks[k].X[i] = x[i];
	//printf("k = %6d, x0 = %9.9lf, x1 = %9.9lf, x2 = %9.9lf, x3 = %9.9lf, x4 = %9.9lf\n", k, 
	//	tracks[k].X[0], tracks[k].X[1], tracks[k].X[2], tracks[k].X[3], tracks[k].X[4]);

	//tracks[k].valid = crossCorrelationInHomography(k)>cross_correlation_threshold_homography;

	return ret;
}
