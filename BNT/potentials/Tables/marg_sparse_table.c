/* marg_sparse_table.c ../potential/tables*/

/******************************************/
/* 5 input & 1 output                     */
/* Big sparse table                       */
/* Big domain                             */
/* Big sizes                              */
/* onto                                   */
/* maximize, if missed, maximize=0        */
/*                                        */
/* small sparse table                     */
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

mxArray* convert_table_to_sparse(const double *Table, const int *sequence, const int nzCounts, const int N){
	mxArray *spTable;
	int     i, temp, *irs, *jcs, count=0;
	double  *sr;

	spTable = mxCreateSparse(N, 1, nzCounts, mxREAL);
    sr  = mxGetPr(spTable);
    irs = mxGetIr(spTable);
    jcs = mxGetJc(spTable);

	jcs[0] = 0;
	jcs[1] = nzCounts;

	for(i=0; i<nzCounts; i++){
		irs[i] = sequence[count];
		count++;
		temp = sequence[count];
		sr[i] = Table[temp];
		count++;
	}
	return spTable;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int        i, j, count, bdim, sdim, NS, NZB, position, bindex, sindex, maximize, nzCounts=0;
	int        *mask, *sequence, *result, *bir, *bjc, *ssize, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double     *sTable, *pbDomain, *psDomain, *pbSize, *bpr, *spr;
	const char *field_names[] = {"domain", "T", "sizes"};

	if(nrhs < 5) maximize = 0;
	else maximize = (int)mxGetScalar(prhs[4]);

	bdim = mxGetNumberOfElements(prhs[1]);
	sdim = mxGetNumberOfElements(prhs[3]);
	pbSize = mxGetPr(prhs[2]);
	pbDomain = mxGetPr(prhs[1]);
	psDomain = mxGetPr(prhs[3]);
	bpr = mxGetPr(prhs[0]);
	bir = mxGetIr(prhs[0]);
	bjc = mxGetJc(prhs[0]);
	NZB = bjc[1];

	if(sdim == 0){
		plhs[0] = mxCreateSparse(1, 1, 1, mxREAL);
		spr = mxGetPr(plhs[0]);
		bir = mxGetIr(plhs[0]);
		bjc = mxGetJc(plhs[0]);
		*spr = 0;
		*bir = 0;
		bjc[0] = 0;
		bjc[1] = 1;
		if(maximize){
			for(i=0; i<NZB; i++){
				*spr = (*spr < bpr[i])? bpr[i] : *spr;
			}
		}
		else{
			for(i=0; i<NZB; i++){
				*spr += bpr[i];
			}
		}	
		return;
	}

	mask = malloc(sdim * sizeof(int));
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
		
	sTable = malloc(NZB * sizeof(double));
	sequence = malloc(NZB * 2 * sizeof(double));
	bCumprod = malloc(bdim * sizeof(int));
	sCumprod = malloc(sdim * sizeof(int));
	bsubv = malloc(bdim * sizeof(int));
	ssubv = malloc(sdim * sizeof(int));
	ssize = malloc(sdim * sizeof(int));

	NS = 1;
	for(i=0; i<count; i++){
		ssize[i] = (int)pbSize[mask[i]];
		NS *= ssize[i];
	}

	for(i=0; i<NZB; i++)sTable[i] = 0;
	
	bCumprod[0] = 1;
	for(i=0; i<bdim-1; i++){
		bCumprod[i+1] = bCumprod[i] * (int)pbSize[i];
	}
	sCumprod[0] = 1;
	for(i=0; i<sdim-1; i++){
		sCumprod[i+1] = sCumprod[i] * ssize[i];
	}

	count = 0;
	for(i=0; i<NZB; i++){
		bindex = bir[i];
		ind_subv(bindex, bCumprod, bdim, bsubv);
		for(j=0; j<sdim; j++){
			ssubv[j] = bsubv[mask[j]];
		}
		sindex = subv_ind(sdim, sCumprod, ssubv);
		result = (int *) bsearch(&sindex, sequence, nzCounts, sizeof(int)*2, compare);
		if(result){
			position = (result - sequence) / 2;
			if(maximize) 
				sTable[position] = (sTable[position] < bpr[i]) ? bpr[i] : sTable[position];
			else sTable[position] += bpr[i];
		}
		else {
			if(maximize) 
				sTable[nzCounts] = (sTable[nzCounts] < bpr[i]) ? bpr[i] : sTable[nzCounts];
			else sTable[nzCounts] += bpr[i];
			sequence[count] = sindex;
			count++;
			sequence[count] = nzCounts;
			nzCounts++;
			count++;
		}
	}
	
	qsort(sequence, nzCounts, sizeof(int) * 2, compare);
	plhs[0] = convert_table_to_sparse(sTable, sequence, nzCounts, NS);

	free(sTable);
	free(sequence);
	free(mask);
	free(bCumprod);
	free(sCumprod);
	free(bsubv);
	free(ssubv);
	free(ssize);
}
