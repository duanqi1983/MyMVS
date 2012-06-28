#ifndef _MATRIX_HPP
#define _MATRIX_HPP

#include "kutility/kutility.def"

#define equals(a,b) fabs(a-b)<0.000000001

inline void matrixMul(double *a, int a_row, int a_col, double *b, int b_row, int b_col, double *c)
{
	for (int i=0;i<a_row;i++)
		for (int j=0;j<b_col;j++)
		{
			c[i*b_col+j] = 0;
			for (int k=0;k<a_col;k++) c[i*b_col+j] += a[i*a_col+k] * b[k*b_col+j];
		}
}

template <class T>
inline double matrixLength(T *a, int size)
{
	double res = 0;
	for (int i=0;i<size;i++) res+=a[i]*a[i];
	return sqrt(res);
}

inline void matrixDiv(double *a, int size, double x)
{
	for (int i=0;i<size;i++) a[i]/=x;
}

inline void matrixNorm(double *a, int size)
{
	matrixDiv(a,size,matrixLength(a,size));
}

inline void matrixCopy(double *src, double *dst, int size)
{
	for (int i=0;i<size;i++) dst[i] = src[i];
}

inline void matrixCopy(int *src, int *dst, int size)
{
	for (int i=0;i<size;i++) dst[i] = src[i];
}

inline void matrixProject(const double *X, double *P, double *y, int swap)
{
	double x[3];
	x[0] = P[0]*X[0] + P[1]*X[1] + P[2]*X[2] + P[3];
	x[1] = P[4]*X[0] + P[5]*X[1] + P[6]*X[2] + P[7];
	x[2] = P[8]*X[0] + P[9]*X[1] + P[10]*X[2] + P[11];
	if (swap==1)
	{
		y[0] = x[1]/x[2];
		y[1] = x[0]/x[2];
	} else
	{
		y[0] = x[0]/x[2];
		y[1] = x[1]/x[2];
	}
}


inline void matrixProject(double *X, double *P, double *y, int swap)
{
	double x[3];
	x[0] = P[0]*X[0] + P[1]*X[1] + P[2]*X[2] + P[3];
	x[1] = P[4]*X[0] + P[5]*X[1] + P[6]*X[2] + P[7];
	x[2] = P[8]*X[0] + P[9]*X[1] + P[10]*X[2] + P[11];
	if (swap==1)
	{
		y[0] = x[1]/x[2];
		y[1] = x[0]/x[2];
	} else
	{
		y[0] = x[0]/x[2];
		y[1] = x[1]/x[2];
	}
}

inline void matrixProject(float *X, double *P, double *y, int swap)
{
	double x[3];
	x[0] = P[0]*X[0] + P[1]*X[1] + P[2]*X[2] + P[3];
	x[1] = P[4]*X[0] + P[5]*X[1] + P[6]*X[2] + P[7];
	x[2] = P[8]*X[0] + P[9]*X[1] + P[10]*X[2] + P[11];
	if (swap==1)
	{
		y[0] = x[1]/x[2];
		y[1] = x[0]/x[2];
	} else
	{
		y[0] = x[0]/x[2];
		y[1] = x[1]/x[2];
	}
}

inline void matrixProject(double *X, double *P, int *y, int swap)
{
	double x[3];
	x[0] = P[0]*X[0] + P[1]*X[1] + P[2]*X[2] + P[3];
	x[1] = P[4]*X[0] + P[5]*X[1] + P[6]*X[2] + P[7];
	x[2] = P[8]*X[0] + P[9]*X[1] + P[10]*X[2] + P[11];
	if (swap==1)
	{
		y[0] = x[1]/x[2];
		y[1] = x[0]/x[2];
	} else
	{
		y[0] = x[0]/x[2];
		y[1] = x[1]/x[2];
	}
}

inline bool matrixBound(double x, int a, int b)
{
	if (x<a) return false;
	if (x>=b) return false;
	return true;
}

inline bool matrixBound(int x, int a, int b)
{
	if (x<a) return false;
	if (x>=b) return false;
	return true;
}

inline double dotProduct(double *a, double *b, int size)
{
	double res = 0;
	for (int i=0; i<size; i++) res+= a[i]*b[i];
	return res;
}

inline double dotProduct(const float *a, double *b, int size)
{
	double res = 0;
	for (int i=0; i<size; i++) res+= a[i]*b[i];
	return res;
}

inline double dotProduct(const float *a, const float *b, int size)
{
	double res = 0;
	for (int i=0; i<size; i++) res+= a[i]*b[i];
	return res;
}

inline double calAngle(double *a, double *b)
{	
	return acos(dotProduct(a,b,3)) * DEGREE;
}

inline double calAngle(const float *a, double *b)
{	
	return acos(dotProduct(a,b,3)) * DEGREE;
}

inline double calAngle(const float *a, const float *b)
{	
	return acos(dotProduct(a,b,3)) * DEGREE;
}

inline void matrixHomographyProject(double *src, double *H, double *dst)
{
	double x[3];
	x[0] = H[0] * src[0] + H[1] * src[1] + H[2];
	x[1] = H[3] * src[0] + H[4] * src[1] + H[5];
	x[2] = H[6] * src[0] + H[7] * src[1] + H[8];
	dst[0] = x[0]/x[2];
	dst[1] = x[1]/x[2];
}

inline double matrixDist(double *x, double *y, int size)
{	
	double l = 0;
	for (int i=0;i<size;i++) l += (x[i]-y[i])*(x[i]-y[i]);
	return sqrt(l);
	
}

inline void matrixAdd(double *src, double *dst, int size)
{
	for (int i=0;i<size;i++) dst[i] += src[i];
}

inline void matrixSubtract(double *src1, double *src2, double *dst, int size)
{
	for (int i=0;i<size;i++) dst[i] = src1[i] - src2[i];
}


// L1 distance
inline double matrixL1(float *thor0, float *thor1, int size)
{
	double l = 0;	
	for (int i=0;i<size;i++) l += fabs(thor0[i]-thor1[i]);
	return l; 
}

//cross_correlation
template <class T>
inline double cross_correlation(T *v1, T *v2, int size)
{
	double v1Norm = matrixLength(v1,size);
	double v2Norm = matrixLength(v2,size);

	double cc = 0;
	for (int i=0;i<size;i++) cc+=v1[i]*v2[i];		
	if (equals(v1Norm,0) || equals(v2Norm,0))
		if (equals(cc,0)) return 1; 
		else return 0;
	cc/=v1Norm;
	cc/=v2Norm;
	return cc;

}

inline void toSpherical(double x, double y, double z, double *z1, double *z2)
{
	*z1 = acos(z);
	*z2 = atan2(y,x);
}

inline void to3DNormal(double z1, double z2, double *x, double *y, double *z)
{
	*x = sin(z1) * cos(z2);
	*y = sin(z1) * sin(z2);
	*z = cos(z1);
}

#endif