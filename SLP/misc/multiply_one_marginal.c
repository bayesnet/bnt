/* multiply_one_marginals.c */
/* the first operand can be a joint marginals of nodes set,*/
/* but the second operand nust be a single node's marginal.*/
/* and the result joint marginal has domain like [prhs[0].domain, prhs[2].domain]*/
/* i.e. cat the second domain at the end of the first domain*/
/* the third operands will be the eff_node_sizes */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray     *ptemp, *ptemp1;
	double      *xdom, *ydom, *xp, *yp, *zp, *pr, *pSizes;
	int         N, xElements, numElements;
	int         i, j, nxdom, ndim, ydim;
	int			*xdim, *s, *sx, *sy, *cpsx, *cpsy;
	int			*subs, *s1, *cpsx2, *cpsy2;
	int         rdims[2];
	const char  *field_names[] = {"domain", "T", "mu", "Sigma"};

	if (nrhs != 3)
		mexErrMsgTxt("Incorrect number of inputs.");
	
	if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
	
	if(mxIsEmpty(prhs[0])){
		plhs[0] = mxDuplicateArray(prhs[1]);
		return;
	}

	N = mxGetNumberOfElements(prhs[2]);
	pSizes = mxGetPr(prhs[2]);

	ptemp = mxGetField(prhs[0], 0, "domain");
	nxdom = mxGetNumberOfElements(ptemp);
	xdom = mxGetPr(ptemp);
	ptemp = mxGetField(prhs[1], 0, "domain");
	ydom = mxGetPr(ptemp);
	ndim = nxdom + 1;
	
	rdims[0] = 1;
	rdims[1] = 1;
	plhs[0] = mxCreateStructArray(2, rdims, 4, field_names);
	ptemp = mxCreateDoubleMatrix(1, ndim, mxREAL);
	mxSetField(plhs[0], 0, "domain", ptemp);
	pr = mxGetPr(ptemp);
	for(i=0; i<nxdom; i++){
		pr[i] = xdom[i];
	}
	pr[ndim-1] = *ydom;

	xdim = (int *)malloc(sizeof(int)*nxdom);
	for(i=0; i<nxdom; i++){
		xdim[i] = (int)pSizes[(int)xdom[i]-1];
	}
	ydim = (int)pSizes[(int)*ydom - 1];

	ptemp = mxGetField(prhs[1], 0, "T");
	yp = mxGetPr(ptemp);
	ptemp = mxGetField(prhs[0], 0, "T");
	xp = mxGetPr(ptemp);
	xElements = mxGetNumberOfElements(ptemp);
	if(ydim == 1){
		ptemp1 = mxDuplicateArray(ptemp);
		mxSetField(plhs[0], 0, "T", ptemp1);
		free(xdim);
		return;
	}
	numElements = xElements * ydim;

	sx = (int *)malloc(sizeof(int)*ndim);
	sy = (int *)malloc(sizeof(int)*ndim);
	s =  (int *)malloc(sizeof(int)*ndim);
	s1 = (int *)malloc(sizeof(int)*ndim);
	*(cpsx = (int *)malloc(sizeof(int)*ndim)) = 1;
	*(cpsy = (int *)malloc(sizeof(int)*ndim)) = 1;
	subs =   (int *)malloc(sizeof(int)*ndim);
	cpsx2 =  (int *)malloc(sizeof(int)*ndim);
	cpsy2 =  (int *)malloc(sizeof(int)*ndim);
	for(i=0; i<nxdom; i++){
		subs[i] = 0;
		sx[i] = xdim[i];
		sy[i] = 1;
		s[i] = sx[i];
		s1[i] = s[i] - 1;
	}
	subs[ndim-1] = 0;
	sx[ndim-1] = 1;
	sy[ndim-1] = ydim;
	s[ndim-1] = ydim;
	s1[ndim-1] = s[ndim-1] - 1;
				
	for(i=0; i<ndim-1; i++){
		cpsx[i+1] = cpsx[i]*sx[i]--;
		cpsy[i+1] = cpsy[i]*sy[i]--;
		cpsx2[i] = cpsx[i]*sx[i];
		cpsy2[i] = cpsy[i]*sy[i];
	}
	cpsx2[ndim-1] = cpsx[ndim-1]*(--sx[ndim-1]);
	cpsy2[ndim-1] = cpsy[ndim-1]*(--sy[ndim-1]);
	
	ptemp = mxCreateNumericArray(ndim, s, mxDOUBLE_CLASS, mxREAL);
	mxSetField(plhs[0], 0, "T", ptemp);
	zp = mxGetPr(ptemp);

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
	free(xdim);
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
