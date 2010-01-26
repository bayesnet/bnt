/* rep_mult.c  repmat first two operands to the size provided by */
/* the third operand, then perform point multiply                */
/* 3 input, 1 output                                             */
/* C = rep_mult(A, B, sizes)                                     */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double		*xp, *yp, *zp, *pSizes;
	int			xnd, ynd, numElements = 1;
	const int	*xdim, *ydim;
	int         i, j, ndim;
	int			*s, *sx, *sy, *cpsx, *cpsy;
	int			*subs, *s1, *cpsx2, *cpsy2;

	if (nrhs != 3)
		mexErrMsgTxt("Incorrect number of inputs.");
	
	if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
	
	xnd = mxGetNumberOfDimensions(prhs[0]);
	ynd = mxGetNumberOfDimensions(prhs[1]);
	xdim = mxGetDimensions(prhs[0]);
	ydim = mxGetDimensions(prhs[1]);
	ndim = mxGetNumberOfElements(prhs[2]);

	pSizes = mxGetPr(prhs[2]);

	sx = (int *)malloc(sizeof(int)*ndim);
	sy = (int *)malloc(sizeof(int)*ndim);
	s =  (int *)malloc(sizeof(int)*ndim);
	s1 = (int *)malloc(sizeof(int)*ndim);
	*(cpsx = (int *)malloc(sizeof(int)*ndim)) = 1;
	*(cpsy = (int *)malloc(sizeof(int)*ndim)) = 1;
	subs =   (int *)malloc(sizeof(int)*ndim);
	cpsx2 =  (int *)malloc(sizeof(int)*ndim);
	cpsy2 =  (int *)malloc(sizeof(int)*ndim);
	for(i=0; i<ndim; i++){
		subs[i] = 0;
		sx[i] = (i < xnd) ? xdim[i] : 1;
		sy[i] = (i < ynd) ? ydim[i] : 1;
		s[i] = (int)pSizes[i];
		s1[i] = s[i] - 1;
		numElements *= s[i];
	}
				
	for(i=0; i<ndim-1; i++){
		cpsx[i+1] = cpsx[i]*sx[i]--;
		cpsy[i+1] = cpsy[i]*sy[i]--;
		cpsx2[i] = cpsx[i]*sx[i];
		cpsy2[i] = cpsy[i]*sy[i];
	}
	cpsx2[ndim-1] = cpsx[ndim-1]*(--sx[ndim-1]);
	cpsy2[ndim-1] = cpsy[ndim-1]*(--sy[ndim-1]);
	
	plhs[0] = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
	zp = mxGetPr(plhs[0]);
	xp = mxGetPr(prhs[0]);
	yp = mxGetPr(prhs[1]);

	for(j=0; j<numElements; j++){
		*zp++ = *xp * *yp;
		for(i=0; i<ndim; i++){
			if(subs[i] == s1[i]){
				subs[i] = 0;
				if(sx[i])
					xp -= cpsx2[i];
				if(sy[i])
					yp -= cpsy2[i];
			}
			else{
				subs[i]++;
				if(sx[i])
					xp += cpsx[i];
				if(sy[i])
					yp += cpsy[i];
				break;
			}
		}
	}
	free(sx);
	free(sy);
	free(s);
	free(s1);
	free(cpsx);
	free(cpsy);
	free(subs);
	free(cpsx2);
	free(cpsy2);
}
