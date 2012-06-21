// FROM: Li, Jianguo
#ifndef _RECITIFY_HPP
#define _RECITIFY_HPP

enum {HOR_STEREO = 0, VER_STEREO};

// P = K * [R T] is the base formulation for all the operators 

//f: 2*1 or 1*2 vector of focal length
//c: 2*1 or 1*2 vector of principal point
//alpha: skew
//R,T: extrinsic camera parameters
void project_points2(CvMat* X, CvMat*& Xp, CvMat* R, 
					 double* T, double* f, double* c, double alpha)
{
	int n = X->cols;
	double inv_Z = 0;

	CvMat* Y = cvCreateMat(3, n, CV_64F);
	cvMatMul(R, X, Y);
	
	Xp = cvCreateMat(2, n, CV_64F);	
#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
	for(int i=0; i<n; i++)
	{
		cvmSet(Y, 0, i, cvmGet(Y, 0, i) + T[0]);
		cvmSet(Y, 1, i, cvmGet(Y, 1, i) + T[1]);
		cvmSet(Y, 2, i, cvmGet(Y, 2, i) + T[2]);
		
		inv_Z = 1/cvmGet(Y, 2, i);
		cvmSet(Y, 0, i, cvmGet(Y, 0, i)*inv_Z);
		cvmSet(Y, 1, i, cvmGet(Y, 1, i)*inv_Z);

		double tmp = cvmGet(Y, 0, i) + alpha * cvmGet(Y, 1, i);
		cvmSet(Y, 0, i, tmp);

		cvmSet(Xp, 0, i, cvmGet(Y, 0, i)*f[0] + c[0]);
		cvmSet(Xp, 1, i, cvmGet(Y, 1, i)*f[1] + c[1]);
	}
	cvReleaseMat(&Y);
}

void columnWiseMean(CvMat* X, CvMat*& Y)
{
	int nc = X->cols;
	int nr = X->rows;
	Y = cvCreateMat(nr, 1, CV_64F);
	for(int i=0; i<nr; ++i)
	{
		CvMat row;
		cvGetRow(X, &row, i);
		CV_MAT_ELEM(*Y, double, i, 0) = cvMean(&row);
	}
}

//MyRectify: Implemented the rectify algorithm in Bouguet's Matlab toolbox
//K,R,T: camera paramters
//K_l_rect, K_r_rect: rectified intrinsic parameters
//R_l_rect, R_r_rect: rotation matrixes from original images to rectified images
//T_new: new translation between rectified images
//type_stere: horizontal or vertical
void icvRectify(CvMat* K_l, CvMat* K_r, CvMat* R_l, CvMat* R_r, CvMat* T_l, CvMat* T_r, 
				double* D1, double * D2, CvSize image_size, 
				CvMat* K_l_rect, CvMat* K_r_rect, CvMat* R_l_rect, CvMat* R_r_rect,
				CvMat* T_new, int &type_stereo)
{
	double _R[9], Tmp9[9];
	CvMat R_ref = cvMat(3, 3, CV_64F, _R);
	CvMat R_l_t = cvMat(3, 3, CV_64F, Tmp9);
	cvTranspose(R_l, &R_l_t);
	cvMatMul(R_r, &R_l_t, &R_ref);

	double Tmp3[3], _T[3];
	CvMat Tprj = cvMat(3,1, CV_64F, Tmp3);
	CvMat T_ref = cvMat(3,1, CV_64F, _T);
	cvMatMul(&R_ref, T_l, &Tprj);
	cvSub(T_r, &Tprj, &T_ref);

	//Bring the 2 cameras in the same orientation by rotating them "minimally": 
	double _om[3] = {0};
	CvMat om = cvMat(3, 1, CV_64F, _om);
	cvRodrigues2( &R_ref, &om);
	cvScale(&om, &om, -0.5);

	double _r_r[9];
	CvMat r_r = cvMat(3, 3, CV_64F, _r_r);
	cvRodrigues2( &om, &r_r );

	double _r_l[9];
	CvMat r_l = cvMat(3, 3, CV_64F, _r_l);
	cvTranspose(&r_r, &r_l);

	double _t[3];
	CvMat t = cvMat(3, 1, CV_64F, _t);
	cvMatMul(&r_r, &T_ref, &t);

	//Rotate both cameras so as to bring the translation vector in alignment with the (1;0;0) axis:
	double _uu[3] = {0};
	CvMat uu = cvMat(3, 1, CV_64F, _uu);
	if( fabs(_t[0]) > fabs(_t[1]) )
	{
		// Horizontal epipolar lines
		type_stereo = HOR_STEREO;
		_uu[0] = 1;
	}
	else
	{
		// Vertical epipolar lines
		type_stereo = VER_STEREO;
		_uu[1] = 1;
	}

	if( cvDotProduct(&uu, &t) < 0) // switch side of the vector 
		cvScale(&uu, &uu, -1.0);		
	
	double _ww[3];
	CvMat ww = cvMat(3, 1, CV_64F, _ww);
	cvCrossProduct( &t, &uu, &ww );
	// support in-place
	cvNormalize(&ww, &ww);

	double scale = acos(fabs(cvDotProduct(&t, &uu))/(cvNorm(&t) * cvNorm(&uu) ));
	cvScale(&ww, &ww, scale);

	double _R2[9];
	CvMat R2 = cvMat(3, 3, CV_64F, _R2);
	cvRodrigues2(&ww, &R2);

	//Global rotations to be applied to both views:
	double _R_R[9];
	double _R_L[9];
	CvMat R_R = cvMat(3, 3, CV_64F, _R_R);
	CvMat R_L = cvMat(3, 3, CV_64F, _R_L);
	cvMatMul(&R2, &r_r, &R_R);
	cvMatMul(&R2, &r_l, &R_L);

	//Note: R_L is the motion of the points in space
	//So: X2 = R_L*X where X: coord in the old reference frame, X2: coord in the
	//new ref frame.
	cvCopy(&R_L, R_l_rect);
	cvCopy(&R_R, R_r_rect);

	//The resulting rigid motion between the two cameras after image rotations (substitutes of om, R and T):
	// R_new = eye(3);
	// om_new = zeros(3,1);
	cvMatMul(&R_R, &T_ref, T_new);

	// Computation of the *new* intrinsic parameters for both left and right cameras
	double fc_left[2], fc_right[2], cc_left[2], cc_right[2], alpha_c_left, alpha_c_right;
	fc_left[0] = cvmGet(K_l, 0,0);
	fc_left[1] = cvmGet(K_l, 1,1);
	fc_right[0]= cvmGet(K_r, 0,0);
	fc_right[1]= cvmGet(K_r, 1,1);

	cc_left[0] = cvmGet(K_l, 0,2);
	cc_left[1] = cvmGet(K_l, 1,2);
	cc_right[0] = cvmGet(K_r, 0,2);
	cc_right[1] = cvmGet(K_r, 1,2);

	alpha_c_left = cvmGet(K_l, 0,1)/cvmGet(K_l, 0,0);
	alpha_c_right= cvmGet(K_r, 0,1)/cvmGet(K_r, 0,0);

	int ny = image_size.height;
	int nx = image_size.width;
	
	// Vertical focal length *MUST* be the same for both images 
	// (here, we are trying to find a focal length that retains as much information contained in the original distorted images):
	double fc_y_left_new, fc_y_right_new;
	double kc_left = 0; 
	double kc_right = 0;
	if( D1 != NULL )
		kc_left = D1[0];
	if( D2 != NULL )
		kc_right = D2[0];
	if (kc_left < 0) // kc is the distortion
		fc_y_left_new = fc_left[1] * (1 + kc_left*(nx*nx + ny*ny)/(4*fc_left[1]*fc_left[1]));
	else
		fc_y_left_new = fc_left[1];
	
	if (kc_right < 0)
		fc_y_right_new = fc_right[1] * (1 + kc_right*(nx*nx + ny*ny)/(4*fc_right[1]*fc_right[1]));
	else
		fc_y_right_new = fc_right[1];
	double fc_y_new = std::min(fc_y_left_new, fc_y_right_new);

	// For simplicity, let's pick the same value for the horizontal focal length as the vertical focal length (resulting into square pixels):
	double fc_left_new[2], fc_right_new[2];
	fc_left_new[0] = fc_right_new[0] = int(fc_y_new+0.5);
	fc_left_new[1] = fc_right_new[1] = int(fc_y_new+0.5);
	
	// Select the new principal points to maximize the visible area in the rectified images
	double vX[4], vY[4];
	vX[0] = 0; vX[1] = nx-1; vX[2] = nx-1; vX[3] = 0;
	vY[0] = 0; vY[1] = 0; vY[2] =ny-1; vY[3] = ny-1;

	// normalize pixels
	double _vNmlPxlLeft[3*4] = {0};
	double _vNmlPxlRigh[3*4] = {0};
	CvMat vNmlPxlLeft = cvMat(3, 4, CV_64F, _vNmlPxlLeft);
	CvMat vNmlPxlRigh = cvMat(3, 4, CV_64F, _vNmlPxlRigh);
	for(int i=0; i<4; i++)
	{
		_vNmlPxlLeft[0*4+i] = (vX[i] - cc_left[0])/fc_left[0];
		_vNmlPxlLeft[1*4+i] = (vY[i] - cc_left[1])/fc_left[1];
		_vNmlPxlLeft[0*4+i] = _vNmlPxlLeft[0*4+i] - alpha_c_right * _vNmlPxlLeft[1*4+i];
		_vNmlPxlLeft[2*4+i] = 1;

		_vNmlPxlRigh[0*4+i] = (vX[i] - cc_right[0])/fc_right[0];
		_vNmlPxlRigh[1*4+i] = (vY[i] - cc_right[1])/fc_right[1];
		_vNmlPxlRigh[0*4+i] = _vNmlPxlRigh[0*4+i] - alpha_c_right * _vNmlPxlRigh[1*4+i];
		_vNmlPxlRigh[2*4+i] = 1;
	}

	//project_points2
	double zero3[3] = {0};
	double zero2[3] = {0};
	CvMat *Xp_Left, *Xp_Right;
	project_points2(&vNmlPxlLeft, Xp_Left, &R_L, zero3, fc_left_new, zero2, 0);
	project_points2(&vNmlPxlRigh, Xp_Right, &R_R, zero3, fc_right_new, zero2, 0);

	CvMat *mean_pt_left, *mean_pt_right;
	columnWiseMean(Xp_Left, mean_pt_left);
	columnWiseMean(Xp_Right, mean_pt_right);
	
	double cc_left_new[2], cc_right_new[2];
	cc_left_new[0] = (nx-1)/2.0;
	cc_left_new[1] = (ny-1)/2.0;
	cc_right_new[0] = (nx-1)/2.0;
	cc_right_new[1] = (ny-1)/2.0;

	cc_left_new[0] = cc_left_new[0] - cvmGet(mean_pt_left, 0, 0);
	cc_left_new[1] = cc_left_new[1] - cvmGet(mean_pt_left, 1, 0);
	cc_right_new[0] = cc_right_new[0] - cvmGet(mean_pt_right, 0, 0);
	cc_right_new[1] = cc_right_new[1] - cvmGet(mean_pt_right, 1, 0);

	cvReleaseMat(&mean_pt_left);
	cvReleaseMat(&mean_pt_right);
	cvReleaseMat(&Xp_Left);
	cvReleaseMat(&Xp_Right);

	//For simplicity, set the principal points for both cameras to be the average of the two principal points.
	if( type_stereo == HOR_STEREO )
	{
		//-- Horizontal stereo
		double cc_y_new = (cc_left_new[1] + cc_right_new[1])/2;
		cc_left_new[1] = cc_y_new;
		cc_right_new[1] = cc_y_new;
	}
	else
	{
		//-- Vertical stereo
		double cc_x_new = (cc_left_new[0] + cc_right_new[0])/2;
		cc_left_new[0] = cc_x_new;
		cc_right_new[0] = cc_x_new;
	}
	
	// Of course, we do not want any skew or distortion after rectification:
	double alpha_c_left_new = 0;
	double alpha_c_right_new = 0;
	double kc_left_new[5] = {0};
	double kc_right_new[5] = {0};

	//The resulting left and right camera matrices
	cvmSet(K_l_rect, 0, 0, fc_left_new[0]);
	cvmSet(K_l_rect, 0, 1, fc_left_new[0]*alpha_c_left_new);
	cvmSet(K_l_rect, 0, 2, cc_left_new[0]);
	cvmSet(K_l_rect, 1, 0, 0);
	cvmSet(K_l_rect, 1, 1, fc_left_new[1]);
	cvmSet(K_l_rect, 1, 2, cc_left_new[1]);
	cvmSet(K_l_rect, 2, 0, 0);
	cvmSet(K_l_rect, 2, 1, 0);
	cvmSet(K_l_rect, 2, 2, 1);

	cvmSet(K_r_rect, 0, 0, fc_right_new[0]);
	cvmSet(K_r_rect, 0, 1, fc_right_new[0]*alpha_c_right_new);
	cvmSet(K_r_rect, 0, 2, cc_right_new[0]);
	cvmSet(K_r_rect, 1, 0, 0);
	cvmSet(K_r_rect, 1, 1, fc_right_new[1]);
	cvmSet(K_r_rect, 1, 2, cc_right_new[1]);
	cvmSet(K_r_rect, 2, 0, 0);
	cvmSet(K_r_rect, 2, 1, 0);
	cvmSet(K_r_rect, 2, 2, 1);
}

void
icvInitUndistortRectifyMap( const CvMat* A, const CvMat* distCoeffs,
						  const CvMat *R, const CvMat* Ar, CvArr* mapxarr, CvArr* mapyarr )
{
	CV_FUNCNAME( "cvInitUndistortMap" );

	__CV_BEGIN__;

	double a[9], ar[9], r[9], ir[9], k[5]={0,0,0,0,0};
	int coi1 = 0, coi2 = 0;
	CvMat mapxstub, *_mapx = (CvMat*)mapxarr;
	CvMat mapystub, *_mapy = (CvMat*)mapyarr;
	CvMat _a = cvMat( 3, 3, CV_64F, a );
	CvMat _k = cvMat( 4, 1, CV_64F, k );
	CvMat _ar = cvMat( 3, 3, CV_64F, ar );
	CvMat _r = cvMat( 3, 3, CV_64F, r );
	CvMat _ir = cvMat( 3, 3, CV_64F, ir );
	int i, j;
	double fx, fy, u0, v0, k1, k2, k3, p1, p2;
	CvSize size;

	CV_CALL( _mapx = cvGetMat( _mapx, &mapxstub, &coi1 ));
	CV_CALL( _mapy = cvGetMat( _mapy, &mapystub, &coi2 ));

	if( coi1 != 0 || coi2 != 0 )
		CV_ERROR( CV_BadCOI, "The function does not support COI" );

	if( CV_MAT_TYPE(_mapx->type) != CV_32FC1 )
		CV_ERROR( CV_StsUnsupportedFormat, "Both maps must have 32fC1 type" );

	if( !CV_ARE_TYPES_EQ( _mapx, _mapy ))
		CV_ERROR( CV_StsUnmatchedFormats, "" );

	if( !CV_ARE_SIZES_EQ( _mapx, _mapy ))
		CV_ERROR( CV_StsUnmatchedSizes, "" );

	if( A )
	{
		if( !CV_IS_MAT(A) || A->rows != 3 || A->cols != 3  ||
			(CV_MAT_TYPE(A->type) != CV_32FC1 && CV_MAT_TYPE(A->type) != CV_64FC1) )
			CV_ERROR( CV_StsBadArg, "Intrinsic matrix must be a valid 3x3 floating-point matrix" );
		cvConvert( A, &_a );
	}
	else
		cvSetIdentity( &_a );

	if( Ar )
	{
		CvMat Ar33;
		if( !CV_IS_MAT(Ar) || Ar->rows != 3 || (Ar->cols != 3 && Ar->cols != 4) ||
			(CV_MAT_TYPE(Ar->type) != CV_32FC1 && CV_MAT_TYPE(Ar->type) != CV_64FC1) )
			CV_ERROR( CV_StsBadArg, "The new intrinsic matrix must be a valid 3x3 floating-point matrix" );
		cvGetCols( Ar, &Ar33, 0, 3 );
		cvConvert( &Ar33, &_ar );
	}
	else
		cvSetIdentity( &_ar );

	if( !CV_IS_MAT(R) || R->rows != 3 || R->cols != 3  ||
		(CV_MAT_TYPE(R->type) != CV_32FC1 && CV_MAT_TYPE(R->type) != CV_64FC1) )
		CV_ERROR( CV_StsBadArg, "Rotaion/homography matrix must be a valid 3x3 floating-point matrix" );

	if( distCoeffs )
	{
		CV_ASSERT( CV_IS_MAT(distCoeffs) &&
			(distCoeffs->rows == 1 || distCoeffs->cols == 1) &&
			(distCoeffs->rows*distCoeffs->cols*CV_MAT_CN(distCoeffs->type) == 4 ||
			distCoeffs->rows*distCoeffs->cols*CV_MAT_CN(distCoeffs->type) == 5) &&
			(CV_MAT_DEPTH(distCoeffs->type) == CV_64F ||
			CV_MAT_DEPTH(distCoeffs->type) == CV_32F) );
		_k = cvMat( distCoeffs->rows, distCoeffs->cols,
			CV_MAKETYPE(CV_64F, CV_MAT_CN(distCoeffs->type)), k );
		cvConvert( distCoeffs, &_k );
	}
	else
		cvZero( &_k );

	cvConvert( R, &_r );    // rectification matrix
	cvMatMul( &_ar, &_r, &_r ); // Ar*R
	cvInvert( &_r, &_ir );  // inverse: R^-1*Ar^-1

	u0 = a[2]; v0 = a[5];
	fx = a[0]; fy = a[4];
	k1 = k[0]; k2 = k[1]; k3 = k[4];
	p1 = k[2]; p2 = k[3];

	// size = cvGetMatSize(_mapx);
	size = cvSize(_mapx->width, _mapx->height);

#ifdef WITH_OPENMP
#pragma omp parallel for private(i,j)
#endif
	for( i = 0; i < size.height; i++ )
	{
		float* mapx = (float*)(_mapx->data.ptr + _mapx->step*i);
		float* mapy = (float*)(_mapy->data.ptr + _mapy->step*i);
		double _x = i*ir[1] + ir[2], _y = i*ir[4] + ir[5], _w = i*ir[7] + ir[8];

		for( j = 0; j < size.width; j++, _x += ir[0], _y += ir[3], _w += ir[6] )
		{
			double w = 1./_w, x = _x*w, y = _y*w;
			double x2 = x*x, y2 = y*y;
			double r2 = x2 + y2, _2xy = 2*x*y;
			double kr = 1 + ((k3*r2 + k2)*r2 + k1)*r2;
			double u = fx*(x*kr + p1*_2xy + p2*(r2 + 2*x2)) + u0;
			double v = fy*(y*kr + p1*(r2 + 2*y2) + p2*_2xy) + v0; 
			mapx[j] = (float)u;
			mapy[j] = (float)v;
		}
	}
	__CV_END__;
};

#endif

