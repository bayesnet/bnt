/* To compile, type "mex C_quickscore.c" */

#include <stdio.h>
#include "nrutil.h"
#include "nrutil.c"
#include <math.h>
#include "mex.h"

#define MAX(X,Y) (X)>(Y)?(X):(Y)

int two_to_the(int n)
{
  return 1 << n;
}

void int2bin(int num, int nbits, int bits[])
{
  int i, mask;
  mask = 1 << (nbits-1); /* mask = 0010...0 , where the 1 is in col nbits (rightmost = col 1) */
  for (i = 0; i < nbits; i++) {
    bits[i] = ((num & mask) == 0) ? 0 : 1;
    num <<= 1;
  }
}


void quickscore(int ndiseases, int nfindings, const double *fpos, int npos, const double *fneg, int nneg,
                  const double *inhibit, const double *prior, const double *leak, double *prob)
{
  double *Pon, *Poff, **Uon, **Uoff, **post, *pterm, *ptermOff, *ptermOn, temp, p, myp;
  int *bits, nsubsets, *fmask;
  int f, d, i, j, si, size_subset, sign;

  Pon = dvector(0, ndiseases);
  Poff = dvector(0, ndiseases);
  Pon[0] = 1;
  Poff[0] = 0;
  for (i=1; i <= ndiseases; i++) {
    Pon[i] = prior[i-1];
    Poff[i] = 1-Pon[i];
  }

  Uon = dmatrix(0, nfindings-1, 0, ndiseases);
  Uoff = dmatrix(0, nfindings-1, 0, ndiseases);
  d = 0;
  for (f=0; f < nfindings; f++) {
    Uon[f][d] = leak[f];
    Uoff[f][d] = leak[f];
  }
  for (f=0; f < nfindings; f++) {
    for (d=1; d <= ndiseases; d++) {
      Uon[f][d] = inhibit[f + nfindings*(d-1)];
      Uoff[f][d] = 1;
    }
  }
  
  post = dmatrix(0, ndiseases, 0, 1);
  for (d = 0; d <= ndiseases; d++) {
    post[d][0] = 0;
    post[d][1] = 0;
  }
  
  bits = ivector(0, npos-1);
  fmask = ivector(0, nfindings-1);
  pterm = dvector(0, ndiseases);
  ptermOff = dvector(0, ndiseases);
  ptermOn = dvector(0, ndiseases);

  nsubsets = two_to_the(npos);

  for (si = 0; si < nsubsets; si++) {
    int2bin(si, npos, bits);
    for (i=0; i < nfindings; i++) fmask[i] = 0;
    for (i=0; i < nneg; i++) fmask[(int)fneg[i]-1] = 1;
    size_subset = 0;
    for (i=0; i < npos; i++) {
      if (bits[i]) {
	size_subset++;
	fmask[(int)fpos[i]-1] = 1;
      }
    }
    p = 1;
    for (d=0; d <= ndiseases; d++) {
      temp = 1;
      for (j = 0; j < nfindings; j++) {
	if (fmask[j]) temp *= Uoff[j][d];
      }
      ptermOff[d] = temp;

      temp = 1;
      for (j = 0; j < nfindings; j++) {
	if (fmask[j]) temp *= Uon[j][d];
      }
      ptermOn[d] = temp;

      pterm[d] = Poff[d]*ptermOff[d] + Pon[d]*ptermOn[d];
      p *= pterm[d];
    }
    sign = (int) pow(-1, size_subset);
    for (d=0; d <= ndiseases; d++) {
      myp = p / pterm[d];
      post[d][0] += sign*(myp * ptermOff[d]);
      post[d][1] += sign*(myp * ptermOn[d]);
    }
  } /* next si */

  
  for (d=0; d <= ndiseases; d++) {
    post[d][0] *= Poff[d];
    post[d][1] *= Pon[d];
  }
  for (d=0; d <= ndiseases; d++) {
    temp = post[d][0] + post[d][1];
    post[d][0] /= temp;
    post[d][1] /= temp;
    if (d>0) { prob[d-1] = post[d][1]; }
  }

  
  free_dvector(Pon, 0, ndiseases);
  free_dvector(Poff, 0, ndiseases);
  free_dmatrix(Uon, 0, nfindings-1, 0, ndiseases);
  free_dmatrix(Uoff, 0, nfindings-1, 0, ndiseases);
  free_dmatrix(post, 0, ndiseases, 0, 1);
  free_ivector(bits, 0, npos-1);
  free_ivector(fmask, 0, nfindings-1);
  free_dvector(pterm, 0, ndiseases);
  free_dvector(ptermOff, 0, ndiseases);
  free_dvector(ptermOn, 0, ndiseases);
}


void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
                 )
{
  double *fpos, *fneg, *inhibit, *prior, *leak, *prob;
  int npos, nneg, ndiseases, nfindings;
  double *p;

  /* read the input args */
  fpos = mxGetPr(prhs[0]);
  npos = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));

  fneg = mxGetPr(prhs[1]);
  nneg = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));

  inhibit = mxGetPr(prhs[2]); /* inhibit(finding, disease) */
  nfindings = mxGetM(prhs[2]);
  ndiseases = mxGetN(prhs[2]);

  prior = mxGetPr(prhs[3]);

  leak = mxGetPr(prhs[4]);


 /* set the output pointers */
  plhs[0] = mxCreateDoubleMatrix(1, ndiseases, mxREAL);
  prob = mxGetPr(plhs[0]);

  quickscore(ndiseases, nfindings, fpos, npos, fneg, nneg, inhibit, prior, leak, prob);
}
  
