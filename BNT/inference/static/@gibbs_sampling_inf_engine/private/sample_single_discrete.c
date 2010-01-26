#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[]) 
{
  double *y, *pr, *dist;
  int k, distSize;
  double r, cumSum;
  
  plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
  y = mxGetPr (plhs[0]);

  pr = mxGetPr (prhs[0]);
  r = pr[0];

  dist = mxGetPr (prhs[1]);
  distSize = mxGetNumberOfElements (prhs[1]);

  for (k = 0, cumSum = 0; (k < distSize) && (r >= cumSum); cumSum += dist[k], k++);

  y[0] = k;
}
