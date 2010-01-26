#include "mex.h"

/* Helper function that extracts a one-dimensional slice from a cpt */
/*
void multiplySlice(mxArray *bnet, mxArray *state, int i, int nsi, int j,
		   mxArray *strides, mxArray *fam, mxArray *cpts,
		   double *y)
*/
void multiplySlice(const mxArray *bnet, const mxArray *state, int i, int nsi, int j,
		   const mxArray *strides, const mxArray *fam, const mxArray *cpts,
		   double *y)
{
  mxArray *ec, *cpt, *family;
  double *ecElts, *cptElts, *famElts, *strideElts, *ev;
  int c1, k, famSize, startInd, strideStride, pos, stride;
  
  strideStride = mxGetM(strides);
  strideElts = mxGetPr(strides);

  ev = mxGetPr(state);

  /* Get the CPT */
  ec = mxGetField (bnet, 0, "equiv_class");
  ecElts = mxGetPr(ec);
  k = (int) ecElts[j-1];
  cpt = mxGetCell (cpts, k-1);
  cptElts = mxGetPr (cpt);

  /* Get the family vector for this cpt */
  family = mxGetCell (fam, j-1);
  famSize = mxGetNumberOfElements (family);
  famElts = mxGetPr (family);

  /* Figure out starting position and stride */
  startInd = 0;
  for (c1 = 0, pos = k-1; c1 < famSize; c1++, pos +=strideStride) {
    if (famElts[c1] != i) {
      startInd += strideElts[pos]*(ev[(int)famElts[c1]-1]-1);
    }
    else {
      stride = strideElts[pos];
    }
  }

  for (c1 = 0, pos = startInd; c1 < nsi; c1++, pos+=stride) {
    y[c1] *= cptElts[pos];
  }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
		 *prhs[])
{
  double *pi, *nsElts, *y, *childrenElts;
  mxArray *ns, *children;
  double sum;
  int i, nsi, c1, numChildren;

  pi = mxGetPr(prhs[2]);
  i = (int) pi[0];

  ns = mxGetField(prhs[0], 0, "node_sizes");
  nsElts = mxGetPr(ns);
  nsi = (int) nsElts[i-1];

  /* Initialize the posterior */
  plhs[0] = mxCreateDoubleMatrix (1, nsi, mxREAL);
  y = mxGetPr(plhs[0]);
  for (c1 = 0; c1 < nsi; c1++) {
    y[c1] = 1;
  }

  /* Multiply in the cpt of the node i */
  multiplySlice(prhs[0], prhs[1], i, nsi, i, prhs[3], prhs[4],
		prhs[6], y);


  /* Multiply in cpts of children of i */
  children = mxGetCell (prhs[5], i-1);
  numChildren = mxGetNumberOfElements (children);
  childrenElts = mxGetPr (children);
  
  for (c1 = 0; c1 < numChildren; c1++) {
    int j;
    j = (int) childrenElts[c1];
    multiplySlice (prhs[0], prhs[1], i, nsi, j, prhs[3], prhs[4],
		   prhs[6], y);
  }

  sum = 0;
  /* normalize! */
  for (c1 = 0; c1 < nsi; c1++) {
    sum += y[c1];
  }

  for (c1 = 0; c1 < nsi; c1++) {
    y[c1] /= sum;
  }
}








