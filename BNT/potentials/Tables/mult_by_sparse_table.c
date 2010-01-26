/* mult_by_sparse_table.c ../potential/tables*/


/******************************************/
/* 6 input & 1 output                     */
/* Big table    [0]                       */
/* Big domain   [1]                       */
/* big sizes    [2]                       */
/* Small table  [3]                       */
/* small domain [4]                       */
/* small sizes  [5]                       */
/*                                        */
/* New big table[0]                       */
/******************************************/

#include <math.h>
#include <stdlib.h>
#include "mex.h"

int compare(const void* src1, const void* src2){
	int i1 = *(int*)src1 ;
	int i2 = *(int*)src2 ;
	return i1-i2 ;
}

void ind_subv(int index, const int *cumprod, int n, int *bsubv){
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
	int     i, j, count, bdim, sdim, NB, NZB, NZS, position, bindex, sindex, nzCounts=0;
	int     *mask, *result, *bir, *sir, *rir, *bjc, *sjc, *rjc, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double  *pbDomain, *psDomain, *pbSize, *psSize, *bpr, *spr, *rpr;

	pbDomain = mxGetPr(prhs[1]);
	bdim = mxGetNumberOfElements(prhs[1]);
	psDomain = mxGetPr(prhs[4]);
	sdim = mxGetNumberOfElements(prhs[4]);

	pbSize = mxGetPr(prhs[2]);
	psSize = mxGetPr(prhs[5]);

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
	}

	bpr = mxGetPr(prhs[0]);
	bir = mxGetIr(prhs[0]);
	bjc = mxGetJc(prhs[0]);
	NZB = bjc[1];

	spr = mxGetPr(prhs[3]);
	sir = mxGetIr(prhs[3]);
	sjc = mxGetJc(prhs[3]);
	NZS = sjc[1];

	plhs[0] = mxDuplicateArray(prhs[0]);
	rpr = mxGetPr(plhs[0]);
	rir = mxGetIr(plhs[0]);
	rjc = mxGetJc(plhs[0]);
	rjc[0] = 0;
	rjc[1] = NZB;

	if(sdim == 0){
		for(i=0; i<NZB; i++){
			rpr[i] *= *spr;
		}	
		return;
	}

	mask = malloc(sdim * sizeof(int));
	bCumprod = malloc(bdim * sizeof(int));
	sCumprod = malloc(sdim * sizeof(int));
	bsubv = malloc(bdim * sizeof(int));
	ssubv = malloc(sdim * sizeof(int));

	count = 0;
	for(i=0; i<sdim; i++){
		for(j=0; j<bdim; j++){
			if(psDomain[i] == pbDomain[j]){
				mask[count] = j;
				count++;
				break;
			}
		}
	}
	
	bCumprod[0] = 1;
	for(i=0; i<bdim-1; i++){
		bCumprod[i+1] = bCumprod[i] * (int)pbSize[i];
	}
	sCumprod[0] = 1;
	for(i=0; i<sdim-1; i++){
		sCumprod[i+1] = sCumprod[i] * (int)psSize[i];
	}

	for(i=0; i<NZB; i++){
		bindex = bir[i];
		ind_subv(bindex, bCumprod, bdim, bsubv);
		for(j=0; j<sdim; j++){
			ssubv[j] = bsubv[mask[j]];
		}
		sindex = subv_ind(sdim, sCumprod, ssubv);
		result = (int *) bsearch(&sindex, sir, NZS, sizeof(int), compare);
		if(result){
			position = result - sir;
			rpr[nzCounts] = bpr[i] * spr[position];
			rir[nzCounts] = bindex;
			nzCounts++;
		}
	}

	reset_nzmax(plhs[0], NZB, nzCounts);
	free(mask);
	free(bCumprod);
	free(sCumprod);
	free(bsubv);
	free(ssubv);
}
