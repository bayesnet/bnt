#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[]) 
{
  double *pn, *pi, *pj, *pm, *y, *ecElts, *pcpt, *famElts, *strideElts,
    *ev, *nsElts;
  int i, k, j, m, n;
  mxArray *ec, *cpt, *fam, *ns;
  int c1, famSize, nsj;
  int strideStride, startInd, stride, pos, numNodes;

  const int BNET = 0;
  const int STATE = 1;
  const int STRIDES = 6;
  const int FAMILIES = 7;
  const int CPT = 8;

  pn = mxGetPr(prhs[3]);
  n = (int) pn[0];
  pi = mxGetPr(prhs[2]);
  i = (int) pi[0];
  pj = mxGetPr(prhs[4]);
  j = (int) pj[0];
  pm = mxGetPr(prhs[5]);
  m = (int) pm[0];
  ev = mxGetPr(prhs[STATE]);
  ns = mxGetField (prhs[BNET], 0, "node_sizes");
  nsElts = mxGetPr (ns);
  numNodes = mxGetM(ns);

  strideStride = mxGetM(prhs[STRIDES]);
  strideElts = mxGetPr(prhs[STRIDES]);


  
  /* Treat the case n = 1 separately */
  if (pn[0] == 1) {

    /* Get the appropriate CPT */
    ec = mxGetField (prhs[BNET], 0, "eclass1");
    ecElts = mxGetPr(ec);
    k = (int) ecElts[i-1];
    cpt = mxGetCell (prhs[8], k-1);
    pcpt = mxGetPr(cpt);

    nsj = (int) nsElts[j-1];

    /* Get the correct family vector */
    /* (Note : MEX is painful) */
    fam = mxGetCell (prhs[FAMILIES], i - 1);
    famSize = mxGetNumberOfElements(fam);
    famElts = mxGetPr(fam);


    /* Figure out starting position and stride */
    startInd = 0;
    for (c1 = 0, pos = k-1; c1 < famSize; c1++, pos+=strideStride) {
      if (famElts[c1] != j) {
	startInd += strideElts[pos]*(ev[(int)famElts[c1]-1]-1);
      }
      else {
	stride = strideElts[pos];
      }
    }
    
    plhs[0] = mxCreateDoubleMatrix (1, nsj, mxREAL);
    y = mxGetPr(plhs[0]);
    for (c1 = 0, pos = startInd; c1 < nsj; c1++, pos+=stride) {
      y[c1] = pcpt[pos];
    }
  }

  /* Handle the case n > 1 */
  else {

    /* Get the appropriate CPT */
    ec = mxGetField (prhs[BNET], 0, "eclass2");
    ecElts = mxGetPr(ec);
    k = (int) ecElts[i-1];
    cpt = mxGetCell (prhs[8], k-1);
    pcpt = mxGetPr(cpt);

    /* Figure out size of slice */
    if (m == 1) {
      nsj = (int) nsElts[j-1];
    }
    else {
      nsj = (int) nsElts[j-1+numNodes];
    }

    /* Figure out family */
    fam = mxGetCell (prhs[FAMILIES], i - 1 + numNodes);
    famSize = mxGetNumberOfElements(fam);
    famElts = mxGetPr(fam);
    
    startInd = 0;
    for (c1 = 0, pos = k-1; c1 < famSize; c1++, pos+=strideStride) {
      int f = (int) famElts[c1];

      if (((f == j+numNodes) && (m == n)) || ((f == j) && (m ==
							    n-1))) {
	stride = strideElts[pos];
      }
      else {
	startInd += strideElts[pos] * (ev[f-1+((n-2)*numNodes)]-1);
      }
    }

    plhs[0] = mxCreateDoubleMatrix(1,nsj, mxREAL);
    y = mxGetPr(plhs[0]);
    for (c1 = 0, pos = startInd; c1 < nsj; c1++, pos+=stride) {
      y[c1] = pcpt[pos];
    }
  }
}
