/****************************************************
A = mult_by_array(big, small)
implicitely copies small |big|/|small| times 
and then does element-wise multiplication.

i.e.,
C = repmat(small(:), 1, length(big(:))/length(small(:)))
A = reshape(big(:) .* C(:), size(big))

However, this C version avoids the expense of the repmat.

Written by wei.hu@intel.com, 28 Jan 2002.
/****************************************************/


#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double		*sp, *zp;
  int			i, j, NB, NS, xnd, ynd, ndim;
  const int	*xdim, *ydim;
  int			*s, *sx, *sy, *cpsy, *subs, *cpsy2;
  
  if (nrhs != 2)
    mexErrMsgTxt("Incorrect number of inputs.");
  
  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  
  plhs[0] = mxDuplicateArray(prhs[0]);
  zp = mxGetPr(plhs[0]);
  sp = mxGetPr(prhs[1]);
  
  xnd = mxGetNumberOfDimensions(prhs[0]);
  ynd = mxGetNumberOfDimensions(prhs[1]);
  xdim = mxGetDimensions(prhs[0]);
  ydim = mxGetDimensions(prhs[1]);
  ndim = xnd;
  
  NB = mxGetNumberOfElements(prhs[0]);
  NS = mxGetNumberOfElements(prhs[1]);
  
  if(NS == 1){
    for(i=0; i<NB; i++){
      *zp++ *= *sp;
    }
    return;
  }
  
  if(NS == NB){
    for(i=0; i<NB; i++){
      *zp++ *= *sp++;
    }
    return;
  }
  
  sx = (int *)malloc(sizeof(int)*ndim);
  sy = (int *)malloc(sizeof(int)*ndim);
  s =  (int *)malloc(sizeof(int)*ndim);
  *(cpsy = (int *)malloc(sizeof(int)*ndim)) = 1;
  subs =   (int *)malloc(sizeof(int)*ndim);
  cpsy2 =  (int *)malloc(sizeof(int)*ndim);
  for(i=0; i<ndim; i++){
    subs[i] = 0;
    sx[i] = xdim[i];
    sy[i] = (i < ynd) ? ydim[i] : 1;
    s[i] = sx[i] - 1;
  }
  
  for (i = 0; i < ndim-1; i++){
    cpsy[i+1] = cpsy[i]*sy[i]--;
    cpsy2[i] = cpsy[i]*sy[i];
  }
  cpsy2[ndim-1] = cpsy[ndim-1]*(--sy[ndim-1]);
  
  for(j=0; j<NB; j++){
    *zp++ *= *sp;
    for(i=0; i<ndim; i++){
      if(subs[i] == s[i]){
	subs[i] = 0;
	if(sy[i]) sp -= cpsy2[i];
      }
      else{
	subs[i]++;
	if(sy[i]) sp += cpsy[i];
	break;
      }
    }
  }
  free(sx);
  free(sy);
  free(s);
  free(cpsy);
  free(subs);
  free(cpsy2);
}
