/*

	GENOPS.C
		Generalized arithmetic operators overloading built-in functions.

	written by Douglas M. Schwarz
	schwarz@servtech.com
	26 December 1998
	Last modified: 2 April 1999

	Copyright 1998-1999 by Douglas M. Schwarz.  All rights reserved.

*/


/*

Build MEX file by entering the appropriate command at the MATLAB prompt
(-D<name> option is equivalent to #define <name> in source file):

mex genops.c -DPLUS_MEX -output plus
mex genops.c -DMINUS_MEX -output minus
mex genops.c -DTIMES_MEX -output times
mex genops.c -DRDIVIDE_MEX -output rdivide
mex genops.c -DLDIVIDE_MEX -output ldivide
mex genops.c -DPOWER_MEX -output power
mex genops.c -DEQ_MEX -output eq
mex genops.c -DNE_MEX -output ne
mex genops.c -DLT_MEX -output lt
mex genops.c -DGT_MEX -output gt
mex genops.c -DLE_MEX -output le
mex genops.c -DGE_MEX -output ge

*/

/*	This file has been formatted for a tab equal to 4 spaces. */

#if defined(EQ_MEX) || defined(NE_MEX) || defined(LT_MEX) || defined(GT_MEX) \
		|| defined(LE_MEX) || defined(GE_MEX)
#define	RELOP_MEX
#endif

#include "mex.h"
#include "matrix.h"
#ifdef POWER_MEX
#include <math.h>
#define PI 3.141592653589793
#endif

bool allequal(int, const int *, const int *);
void removeZeroImag(double *, double *, int, const int *, int, mxArray **);

#define	xMat  prhs[0]
#define	yMat  prhs[1]
#define	zMat  plhs[0]

#define	min(A,B)  ((A) < (B) ? (A) : (B))
#define	max(A,B)  ((A) > (B) ? (A) : (B))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double		*xrp, *xip, *yrp, *yip;
#ifndef RELOP_MEX
	double		*zr, *zi, *zip;
#endif
	double		*zrp, *zrend;
	int			xnd, ynd, numElements = 1;
	const int	*xdim, *ydim;
	bool		xcmplx, ycmplx;
	mxClassID	yclass;
	int			*s, ndim, *sx, *sy, i, *cpsx, *cpsy;
	int			*subs, *s1, *cpsx2, *cpsy2;
	int			ix = 0, iy = 0;
	mxArray		*args[3], *result[1];
#if defined(RDIVIDE_MEX) || defined(LDIVIDE_MEX)
	double		denom;
#endif
#ifdef POWER_MEX
	double		mag, theta, phi, magx;
	int			flops = 0;
#endif
	
	
	if (nrhs != 2)
		mexErrMsgTxt("Incorrect number of inputs.");
	
	if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
	
	xnd = mxGetNumberOfDimensions(xMat);
	ynd = mxGetNumberOfDimensions(yMat);
	xdim = mxGetDimensions(xMat);
	ydim = mxGetDimensions(yMat);
	
	yclass = mxGetClassID(yMat);
	
/*	If the built-in function in MATLAB can handle the arguments
	then use that. */
	if (yclass != mxDOUBLE_CLASS || 
		(xnd == 2  &&  xdim[0] == 1  &&  xdim[1] == 1) || 
		(ynd == 2  &&  ydim[0] == 1  &&  ydim[1] == 1) || 
		(xnd == ynd  &&  allequal(xnd,xdim,ydim)))
	{
#ifdef PLUS_MEX
		args[0] = mxCreateString("plus");
#elif defined(MINUS_MEX)
		args[0] = mxCreateString("minus");
#elif defined(TIMES_MEX)
		args[0] = mxCreateString("times");
#elif defined(RDIVIDE_MEX)
		args[0] = mxCreateString("rdivide");
#elif defined(LDIVIDE_MEX)
		args[0] = mxCreateString("ldivide");
#elif defined(POWER_MEX)
		args[0] = mxCreateString("power");
#elif defined(EQ_MEX)
		args[0] = mxCreateString("eq");
#elif defined(NE_MEX)
		args[0] = mxCreateString("ne");
#elif defined(LT_MEX)
		args[0] = mxCreateString("lt");
#elif defined(GT_MEX)
		args[0] = mxCreateString("gt");
#elif defined(LE_MEX)
		args[0] = mxCreateString("le");
#elif defined(GE_MEX)
		args[0] = mxCreateString("ge");
#endif
		args[1] = (mxArray *)xMat;
		args[2] = (mxArray *)yMat;
		mexCallMATLAB(1, result, 3, args, "builtin");
		mxDestroyArray(args[0]);
		zMat = result[0];
	}
	else  /* X and Y are both N-D and different dimensionality. */
	{
		ndim = max(xnd,ynd);
		sx = (int *)mxMalloc(sizeof(int)*ndim);
		sy = (int *)mxMalloc(sizeof(int)*ndim);
		s =  (int *)mxMalloc(sizeof(int)*ndim);
		s1 = (int *)mxMalloc(sizeof(int)*ndim);
		*(cpsx = (int *)mxMalloc(sizeof(int)*ndim)) = 1;
		*(cpsy = (int *)mxMalloc(sizeof(int)*ndim)) = 1;
		subs =   (int *)mxMalloc(sizeof(int)*ndim);
		cpsx2 =  (int *)mxMalloc(sizeof(int)*ndim);
		cpsy2 =  (int *)mxMalloc(sizeof(int)*ndim);
		for (i = 0; i < ndim; i++)
		{
			subs[i] = 0;
			sx[i] = (i < xnd) ? xdim[i] : 1;
			sy[i] = (i < ynd) ? ydim[i] : 1;
			if (sx[i] == sy[i])
				s[i] = sx[i];
			else if (sx[i] == 1)
				s[i] = sy[i];
			else if (sy[i] == 1)
				s[i] = sx[i];
			else
			{
				mxFree(sx);
				mxFree(sy);
				mxFree(s);
				mxFree(s1);
				mxFree(cpsx);
				mxFree(cpsy);
				mxFree(subs);
				mxFree(cpsx2);
				mxFree(cpsy2);
				mexErrMsgTxt("Array dimensions are not appropriate.");
			}
			s1[i] = s[i] - 1;
			numElements *= s[i];
		}
				
		for (i = 0; i < ndim-1; i++)
		{
			cpsx[i+1] = cpsx[i]*sx[i]--;
			cpsy[i+1] = cpsy[i]*sy[i]--;
			cpsx2[i] = cpsx[i]*sx[i];
			cpsy2[i] = cpsy[i]*sy[i];
		}
		cpsx2[ndim-1] = cpsx[ndim-1]*(--sx[ndim-1]);
		cpsy2[ndim-1] = cpsy[ndim-1]*(--sy[ndim-1]);
		
		xcmplx = mxIsComplex(xMat);
		ycmplx = mxIsComplex(yMat);
		
		if (!xcmplx && !ycmplx)  /* X and Y both N-D, both real. */
		{
#ifdef POWER_MEX
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxCOMPLEX);
			zrp = zr = mxGetPr(zMat);
			zip = zi = mxGetPi(zMat);
#elif defined(RELOP_MEX)
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
			mxSetLogical(zMat);
			zrp = mxGetPr(zMat);
#else
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
			zrp = mxGetPr(zMat);
#endif
			xrp = mxGetPr(xMat);
			yrp = mxGetPr(yMat);
			zrend = zrp + numElements;
			while (zrp < zrend)
			{
#ifdef PLUS_MEX
				*zrp++ = *xrp + *yrp;
#elif defined(MINUS_MEX)
				*zrp++ = *xrp - *yrp;
#elif defined(TIMES_MEX)
				*zrp++ = *xrp * *yrp;
#elif defined(RDIVIDE_MEX)
				*zrp++ = *xrp / *yrp;
#elif defined(LDIVIDE_MEX)
				*zrp++ = *yrp / *xrp;
#elif defined(POWER_MEX)
				if (*xrp < 0.0 && *yrp != floor(*yrp))
				{
					mag = pow(-*xrp,*yrp);
					theta = PI * *yrp;
					*zrp++ = mag*cos(theta);
					*zip++ = mag*sin(theta);
					flops += 18;
				}
				else
				{
					*zrp++ = pow(*xrp,*yrp);
					*zip++ = 0.0;
					flops++;
				}
#elif defined(EQ_MEX)
				*zrp++ = (*xrp == *yrp);
#elif defined(NE_MEX)
				*zrp++ = (*xrp != *yrp);
#elif defined(LT_MEX)
				*zrp++ = (*xrp < *yrp);
#elif defined(GT_MEX)
				*zrp++ = (*xrp > *yrp);
#elif defined(LE_MEX)
				*zrp++ = (*xrp <= *yrp);
#elif defined(GE_MEX)
				*zrp++ = (*xrp >= *yrp);
#endif
				for (i = 0; i < ndim; i++)
				{
					if (subs[i] == s1[i])
					{
						subs[i] = 0;
						if (sx[i])
							xrp -= cpsx2[i];
						if (sy[i])
							yrp -= cpsy2[i];
					}
					else
					{
						subs[i]++;
						if (sx[i])
							xrp += cpsx[i];
						if (sy[i])
							yrp += cpsy[i];
						break;
					}
				}
			}
#ifdef POWER_MEX
			mexAddFlops(flops);
			removeZeroImag(zr, zi, ndim, (const int *)s, numElements, &zMat);
#elif !defined(RELOP_MEX)
			mexAddFlops(numElements);
#endif
		}
		else if (!xcmplx && ycmplx)  /* X and Y both N-D, X real, Y complex. */
		{
#ifdef POWER_MEX
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxCOMPLEX);
			zrp = zr = mxGetPr(zMat);
			zip = zi = mxGetPi(zMat);
#elif defined(RELOP_MEX)
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
			mxSetLogical(zMat);
			zrp = mxGetPr(zMat);
#else
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxCOMPLEX);
			zrp = mxGetPr(zMat);
			zip = mxGetPi(zMat);
#endif
			xrp = mxGetPr(xMat);
			yrp = mxGetPr(yMat);
			yip = mxGetPi(yMat);
			zrend = zrp + numElements;
			while (zrp < zrend)
			{
#ifdef PLUS_MEX
				*zrp++ = *xrp + *yrp;
				*zip++ = *yip;
#elif defined(MINUS_MEX)
				*zrp++ = *xrp - *yrp;
				*zip++ = -*yip;
#elif defined(TIMES_MEX)
				*zrp++ = *xrp * *yrp;
				*zip++ = *xrp * *yip;
#elif defined(RDIVIDE_MEX)
				denom = *yrp * *yrp + *yip * *yip;
				*zrp++ = (*xrp * *yrp)/denom;
				*zip++ = (-*xrp * *yip)/denom;
#elif defined(LDIVIDE_MEX)
				*zrp++ = *yrp / *xrp;
				*zip++ = *yip / *xrp;
#elif defined(POWER_MEX)
				if (*yip == 0.0)
				{
					if (*xrp < 0.0 && *yrp != floor(*yrp))
					{
						mag = pow(-*xrp,*yrp);
						theta = PI * *yrp;
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 18;
					}
					else
					{
						*zrp++ = pow(*xrp,*yrp);
						*zip++ = 0.0;
						flops++;
					}
				}
				else
				{
					if (*xrp < 0.0)
					{
						mag = pow(-*xrp,*yrp)*exp(-PI * *yip);
						theta = *yip * log(-*xrp) + PI * *yrp;
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 18;
					}
					else
					{
						mag = pow(*xrp,*yrp);
						theta = *yip * log(*xrp);
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 13;
					}
				}
#elif defined(EQ_MEX)
				*zrp++ = (*xrp == *yrp) && (*yip == 0.0);
#elif defined(NE_MEX)
				*zrp++ = (*xrp != *yrp) || (*yip != 0.0);
#elif defined(LT_MEX)
				*zrp++ = (*xrp < *yrp);
#elif defined(GT_MEX)
				*zrp++ = (*xrp > *yrp);
#elif defined(LE_MEX)
				*zrp++ = (*xrp <= *yrp);
#elif defined(GE_MEX)
				*zrp++ = (*xrp >= *yrp);
#endif
				for (i = 0; i < ndim; i++)
				{
					if (subs[i] == s1[i])
					{
						subs[i] = 0;
						if (sx[i])
							xrp -= cpsx2[i];
						if (sy[i])
						{
							yrp -= cpsy2[i];
							yip -= cpsy2[i];
						}
					}
					else
					{
						subs[i]++;
						if (sx[i])
							xrp += cpsx[i];
						if (sy[i])
						{
							yrp += cpsy[i];
							yip += cpsy[i];
						}
						break;
					}
				}
			}
#if defined(PLUS_MEX) || defined(MINUS_MEX)
			mexAddFlops(2*numElements);
#elif defined(TIMES_MEX) || defined(RDIVIDE_MEX) || defined(LDIVIDE_MEX)
			mexAddFlops(6*numElements);
#elif defined(POWER_MEX)
			mexAddFlops(flops);
			removeZeroImag(zr, zi, ndim, (const int *)s, numElements, &zMat);
#endif
		}
		else if (xcmplx && !ycmplx)  /* X and Y both N-D, X complex, Y real. */
		{
#ifdef POWER_MEX
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxCOMPLEX);
			zrp = zr = mxGetPr(zMat);
			zip = zi = mxGetPi(zMat);
#elif defined(RELOP_MEX)
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
			mxSetLogical(zMat);
			zrp = mxGetPr(zMat);
#else
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxCOMPLEX);
			zrp = mxGetPr(zMat);
			zip = mxGetPi(zMat);
#endif
			xrp = mxGetPr(xMat);
			xip = mxGetPi(xMat);
			yrp = mxGetPr(yMat);
			zrend = zrp + numElements;
			while (zrp < zrend)
			{
#ifdef PLUS_MEX
				*zrp++ = *xrp + *yrp;
				*zip++ = *xip;
#elif defined(MINUS_MEX)
				*zrp++ = *xrp - *yrp;
				*zip++ = *xip;
#elif defined(TIMES_MEX)
				*zrp++ = *xrp * *yrp;
				*zip++ = *xip * *yrp;
#elif defined(RDIVIDE_MEX)
				*zrp++ = *xrp / *yrp;
				*zip++ = *xip / *yrp;
#elif defined(LDIVIDE_MEX)
				denom = *xrp * *xrp + *xip * *xip;
				*zrp++ = (*xrp * *yrp)/denom;
				*zip++ = (-*xip * *yrp)/denom;
#elif defined(POWER_MEX)
				if (*xip == 0.0)
				{
					if (*xrp < 0.0 && *yrp != floor(*yrp))
					{
						mag = pow(-*xrp,*yrp);
						theta = PI * *yrp;
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 18;
					}
					else
					{
						*zrp++ = pow(*xrp,*yrp);
						*zip++ = 0.0;
						flops++;
					}
				}
				else
				{
					mag = pow(*xrp * *xrp + *xip * *xip,0.5 * *yrp);
					theta = *yrp*atan2(*xip,*xrp);
					*zrp++ = mag*cos(theta);
					*zip++ = mag*sin(theta);
					flops += 18;
				}
#elif defined(EQ_MEX)
				*zrp++ = (*xrp == *yrp) && (*xip == 0.0);
#elif defined(NE_MEX)
				*zrp++ = (*xrp != *yrp) || (*xip != 0.0);
#elif defined(LT_MEX)
				*zrp++ = (*xrp < *yrp);
#elif defined(GT_MEX)
				*zrp++ = (*xrp > *yrp);
#elif defined(LE_MEX)
				*zrp++ = (*xrp <= *yrp);
#elif defined(GE_MEX)
				*zrp++ = (*xrp >= *yrp);
#endif
				for (i = 0; i < ndim; i++)
				{
					if (subs[i] == s1[i])
					{
						subs[i] = 0;
						if (sx[i])
						{
							xrp -= cpsx2[i];
							xip -= cpsx2[i];
						}
						if (sy[i])
							yrp -= cpsy2[i];
					}
					else
					{
						subs[i]++;
						if (sx[i])
						{
							xrp += cpsx[i];
							xip += cpsx[i];
						}
						if (sy[i])
							yrp += cpsy[i];
						break;
					}
				}
			}
#if defined(PLUS_MEX) || defined(MINUS_MEX)
			mexAddFlops(2*numElements);
#elif defined(TIMES_MEX) || defined(RDIVIDE_MEX) || defined(LDIVIDE_MEX)
			mexAddFlops(6*numElements);
#elif defined(POWER_MEX)
			mexAddFlops(flops);
			removeZeroImag(zr, zi, ndim, (const int *)s, numElements, &zMat);
#endif
		}
		else if (xcmplx && ycmplx)  /* X and Y both N-D, both complex. */
		{
#if defined(RELOP_MEX)
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
			mxSetLogical(zMat);
			zrp = mxGetPr(zMat);
#else
			zMat = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxCOMPLEX);
			zrp = zr = mxGetPr(zMat);
			zip = zi = mxGetPi(zMat);
#endif
			xrp = mxGetPr(xMat);
			xip = mxGetPi(xMat);
			yrp = mxGetPr(yMat);
			yip = mxGetPi(yMat);
			zrend = zrp + numElements;
			while (zrp < zrend)
			{
#ifdef PLUS_MEX
				*zrp++ = *xrp + *yrp;
				*zip++ = *xip + *yip;
#elif defined(MINUS_MEX)
				*zrp++ = *xrp - *yrp;
				*zip++ = *xip - *yip;
#elif defined(TIMES_MEX)
				*zrp++ = *xrp * *yrp - *xip * *yip;
				*zip++ = *xip * *yrp + *xrp * *yip;
#elif defined(RDIVIDE_MEX)
				denom = *yrp * *yrp + *yip * *yip;
				*zrp++ = (*xrp * *yrp + *xip * *yip)/denom;
				*zip++ = (*xip * *yrp - *xrp * *yip)/denom;
#elif defined(LDIVIDE_MEX)
				denom = *xrp * *xrp + *xip * *xip;
				*zrp++ = (*xrp * *yrp + *xip * *yip)/denom;
				*zip++ = (*xrp * *yip - *xip * *yrp)/denom;
#elif defined(POWER_MEX)
				if (*xip == 0.0 && *yip == 0.0)
				{
					if (*xrp < 0.0 && *yrp != floor(*yrp))
					{
						mag = pow(-*xrp,*yrp);
						theta = PI * *yrp;
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 18;
					}
					else
					{
						*zrp++ = pow(*xrp,*yrp);
						*zip++ = 0.0;
						flops++;
					}
				}
				else if (*xip == 0.0)
				{
					if (*xrp < 0.0)
					{
						mag = pow(-*xrp,*yrp)*exp(-PI * *yip);
						theta = *yip * log(-*xrp) + PI * *yrp;
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 18;
					}
					else
					{
						mag = pow(*xrp,*yrp);
						theta = *yip * log(*xrp);
						*zrp++ = mag*cos(theta);
						*zip++ = mag*sin(theta);
						flops += 13;
					}
				}
				else if (*yip == 0.0)
				{
					mag = pow(*xrp * *xrp + *xip * *xip,0.5 * *yrp);
					theta = *yrp * atan2(*xip,*xrp);
					*zrp++ = mag*cos(theta);
					*zip++ = mag*sin(theta);
					flops += 18;
				}
				else
				{
					magx = sqrt(*xrp * *xrp + *xip * *xip);
					phi = atan2(*xip,*xrp);
					mag = pow(magx,*yrp)*exp(-*yip * phi);
					theta = *yip * log(magx) + *yrp * phi;
					*zrp++ = mag*cos(theta);
					*zip++ = mag*sin(theta);
					flops += 18;
				}
#elif defined(EQ_MEX)
				*zrp++ = (*xrp == *yrp) && (*xip == *yip);
#elif defined(NE_MEX)
				*zrp++ = (*xrp != *yrp) || (*xip != *yip);
#elif defined(LT_MEX)
				*zrp++ = (*xrp < *yrp);
#elif defined(GT_MEX)
				*zrp++ = (*xrp > *yrp);
#elif defined(LE_MEX)
				*zrp++ = (*xrp <= *yrp);
#elif defined(GE_MEX)
				*zrp++ = (*xrp >= *yrp);
#endif
				for (i = 0; i < ndim; i++)
				{
					if (subs[i] == s1[i])
					{
						subs[i] = 0;
						if (sx[i])
						{
							xrp -= cpsx2[i];
							xip -= cpsx2[i];
						}
						if (sy[i])
						{
							yrp -= cpsy2[i];
							yip -= cpsy2[i];
						}
					}
					else
					{
						subs[i]++;
						if (sx[i])
						{
							xrp += cpsx[i];
							xip += cpsx[i];
						}
						if (sy[i])
						{
							yrp += cpsy[i];
							yip += cpsy[i];
						}
						break;
					}
				}
			}
#if defined(PLUS_MEX) || defined(MINUS_MEX)
			mexAddFlops(2*numElements);
#elif defined(TIMES_MEX) || defined(RDIVIDE_MEX) || defined(LDIVIDE_MEX)
			mexAddFlops(6*numElements);
#elif defined(POWER_MEX)
			mexAddFlops(flops);
#endif
#ifndef RELOP_MEX
			removeZeroImag(zr, zi, ndim, (const int *)s, numElements, &zMat);
#endif
		}
	}
}


/***********************************************************
*                                                          *
*   Tests to see if the vectors xdim and ydim are equal.   *
*                                                          *
***********************************************************/
bool allequal(int ndim, const int *xdim, const int *ydim)
{
	int		i;
	bool	result = true;
	
	for (i = 0; i < ndim; i++)
		result = result && (xdim[i] == ydim[i]);
	
	return(result);
}


/******************************************************************************
*                                                                             *
*   Tests to see if every imaginary element is identically zero and, if so,   *
*   creates a new array which is real and copies the real elements to it.     *
*                                                                             *
******************************************************************************/
void removeZeroImag(double *zr, double *zi, int ndim, const int *s,
					int numElements, mxArray *plhs[])
{
	double			*zrend, *ziend, *zip, *z1p, *z2p;
	bool			allImZero = true;
	mxArray			*temp;
	
	zip = zi;
	ziend = zi + numElements;
	while (zip < ziend)
	{
		allImZero = allImZero && (*zip++ == 0.0);
		if (!allImZero)
			return;
	}
	
	temp = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
	z1p = zr;
	z2p = mxGetPr(temp);
	zrend = z1p + numElements;
	while (z1p < zrend)
		*z2p++ = *z1p++;
	mxDestroyArray(plhs[0]);
	plhs[0] = temp;
	return;
}
