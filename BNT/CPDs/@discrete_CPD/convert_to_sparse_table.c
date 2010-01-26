/* convert_to_sparse_table.c  convert a sparse discrete CPD with evidence into sparse table */
/* convert_to_pot.m located in ../CPDs/discrete_CPD call it */
/* 3 input */
/* CPD      prhs[0] with 1D sparse CPT */
/* domain   prhs[1]                    */
/* evidence prhs[2]                    */
/* 1 output */
/* T        plhs[0] sparse table       */

#include <math.h>
#include "mex.h"

void ind_subv(int index, const int *cumprod, const int n, int *bsubv){
	int i;

	for (i = n-1; i >= 0; i--) {
		bsubv[i] = ((int)floor(index / cumprod[i]));
		index = index % cumprod[i];
	}
}

int subv_ind(const int n, const int *cumprod, const int *subv){
	int i, index=0;

	for(i=0; i<n; i++){
		index += subv[i] * cumprod[i];
	}
	return index;
}

void reset_nzmax(mxArray *spArray, const int old_nzmax, const int new_nzmax){
	double *ptr;
	void   *newptr;
	int    *ir, *jc;
	int    nbytes;

	if(new_nzmax == old_nzmax) return;
	nbytes = new_nzmax * sizeof(*ptr);
	ptr = mxGetPr(spArray);
	newptr = mxRealloc(ptr, nbytes);
	mxSetPr(spArray, newptr);
	nbytes = new_nzmax * sizeof(*ir);
	ir = mxGetIr(spArray);
	newptr = mxRealloc(ir, nbytes);
	mxSetIr(spArray, newptr);
	jc = mxGetJc(spArray);
	jc[0] = 0;
	jc[1] = new_nzmax;
	mxSetNzmax(spArray, new_nzmax);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int     i, j, NS, NZB, count, bdim, match, domain, bindex, sindex, nzCounts=0;
	int     *observed, *bsubv, *ssubv, *bir, *sir, *bjc, *sjc, *mask, *ssize, *bcumprod, *scumprod;
	double  *pDomain, *pSize, *bpr, *spr;
	mxArray *pTemp;

	pTemp = mxGetField(prhs[0], 0, "CPT");
	bpr = mxGetPr(pTemp);
	bir = mxGetIr(pTemp);
	bjc = mxGetJc(pTemp);
	NZB = bjc[1];
	pTemp = mxGetField(prhs[0], 0, "sizes");
	pSize = mxGetPr(pTemp);

	pDomain = mxGetPr(prhs[1]);
	bdim = mxGetNumberOfElements(prhs[1]);

	mask = malloc(bdim * sizeof(int));
	ssize = malloc(bdim * sizeof(int));
	observed = malloc(bdim * sizeof(int));

	for(i=0; i<bdim; i++){
		ssize[i] = (int)pSize[i];
	}

	count = 0;
	for(i=0; i<bdim; i++){
		domain = (int)pDomain[i] - 1;
		pTemp = mxGetCell(prhs[2], domain);
		if(pTemp){
			mask[count] = i;
			ssize[i] = 1;
			observed[count] = (int)mxGetScalar(pTemp) - 1;
			count++;
		}
	}

	if(count == 0){
		pTemp = mxGetField(prhs[0], 0, "CPT");
		plhs[0] = mxDuplicateArray(pTemp);
		free(mask);
		free(ssize);
		free(observed);
		return;
	}

	bsubv = malloc(bdim * sizeof(int));
	ssubv = malloc(count * sizeof(int));
	bcumprod = malloc(bdim * sizeof(int));
	scumprod = malloc(bdim * sizeof(int));

	NS = 1;
	for(i=0; i<bdim; i++){
		NS *= ssize[i];
	}

	plhs[0] = mxCreateSparse(NS, 1, NS, mxREAL);
	spr = mxGetPr(plhs[0]);
	sir = mxGetIr(plhs[0]);
	sjc = mxGetJc(plhs[0]);
	sjc[0] = 0;
	sjc[1] = NS;

	bcumprod[0] = 1;
	scumprod[0] = 1;
	for(i=0; i<bdim-1; i++){
		bcumprod[i+1] = bcumprod[i] * (int)pSize[i];
		scumprod[i+1] = scumprod[i] * ssize[i];
	}

	nzCounts = 0;
	for(i=0; i<NZB; i++){
		bindex = bir[i];
		ind_subv(bindex, bcumprod, bdim, bsubv);
		for(j=0; j<count; j++){
			ssubv[j] = bsubv[mask[j]];
		}
		match = 1;
		for(j=0; j<count; j++){
			if((ssubv[j]) != observed[j]){
				match = 0;
				break;
			}
		}
		if(match){
			spr[nzCounts] = bpr[i];
			sindex = subv_ind(bdim, scumprod, bsubv);
			sir[nzCounts] = sindex;
			nzCounts++;
		}
	}

	reset_nzmax(plhs[0], NS, nzCounts);
	free(mask);
	free(ssize);
	free(observed);
	free(bsubv);
	free(ssubv);
	free(bcumprod);
	free(scumprod);
}

