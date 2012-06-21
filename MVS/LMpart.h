
void triangulateTrack(track *trk)
{	
	// direct linear transform method
	cv::Mat A = cv::Mat(trk->size*2,4,CV_64F);
	for (int i=0;i<trk->size;i++)
	{	
		double *P = listViewPoints[trk->views[i]]->P;
		for (int j=0;j<4;j++)
		{
			A.at<double>(i*2  ,j) = trk->py[i]*P[8+j] - P[j]; 
			A.at<double>(i*2+1,j) = trk->px[i]*P[8+j] - P[4+j];
		}		
		A.row(i*2)/= norm(A.row(i*2));
		A.row(i*2+1)/= norm(A.row(i*2+1));
	}
	cv::SVD svd = cv::SVD(A);
	for (int i=0;i<3;i++) trk->X[i] = svd.vt.at<double>(3,i)/svd.vt.at<double>(3,3);
}


int LMA_function(void *p, int m, int n, const double *X, double *fvec, int iflag)
{

	track* trk = (track *)p;
	trk->reliableViews = 0;

	double x[2];

	double absx, absy;
	for (int i = 0; i < m; i++)
	{		
		matrixProject(X,listViewPoints[trk->views[i]]->P,x,1);
		absx = fabs(x[0] - trk->px[i]);
		absy = fabs(x[1] - trk->py[i]);
		fvec[i] = sqrt(absx*absx + absy*absy);

		if (matrixBound(x[0],0,listViewPoints[trk->views[i]]->height) &&
			matrixBound(x[1],0,listViewPoints[trk->views[i]]->width) &&
			absx<=gamma && absy <= gamma && fvec[i]<gamma)
		{			
			trk->reliableViews ++;
			trk->reliable[i] = true;
		} else
		{
			fvec[i] *= 10;
			trk->reliable[i] = false;
		}
	}

	return 0;
}

void optimizeByLMA(track *trk)
{
	int m = trk->size ; //no. of observation variables		
	int n = 3; //no. of paremeters 	

	double* fvec = new double[m]; //no need to populate 
	double ftol = 1e-08; //tolerance
	double xtol = 1e-08; //tolerance
	double gtol = 1e-08; //tolerance
	int maxfev = 400; //maximum function evaluations
	double epsfcn=1e-08; //tolerance
	double* diag=new double[n]; //some internal thing
	int mode=1; //some internal thing
	double factor=1; // a default recommended value
	int nprint=0; //don't know what it does
	int info=0; //output variable
	int nfev=0; //output variable will store no. of function evals
	double* fjac=new double[m*n]; //output array of jacobian
	int ldfjac=m; //recommended setting
	int* ipvt=new int[n]; //for internal use
	double* qtf=new double[n]; //for internal use
	double* wa1=new double[n]; //for internal use
	double* wa2=new double[n]; //for internal use
	double* wa3=new double[n]; //for internal use
	double* wa4=new double[m]; //for internal use

	info = lmdif(LMA_function, (void *)trk , m, n, trk->X, fvec, ftol, xtol, gtol, maxfev, epsfcn, 
		diag, mode, factor, nprint, &nfev, fjac, ldfjac, 
		ipvt, qtf, wa1, wa2, wa3, wa4);
	
	deallocate(fvec);
	deallocate(diag);
	deallocate(fjac);
	deallocate(ipvt);
	deallocate(qtf);
	deallocate(wa1);
	deallocate(wa2);
	deallocate(wa3);
	deallocate(wa4);
}
