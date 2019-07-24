/* C mex init_pot for in @jtree_sparse_inf_engine directory               */
/* The file enter_evidence.m in directory @jtree_sparse_inf_engine call it*/

/**************************************/
/* init_pot.c has 6 input & 2 output  */
/* engine                             */
/* clqs                               */
/* pots                               */
/* pot_type                           */
/* onodes                             */
/* ndx                                */
/*                                    */
/* clpot                              */
/* seppot                             */
/**************************************/
#include <math.h>
#include <search.h>
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

void compute_fixed_weight(int *weight, const double *pbSize, const int *dmask, const int *bCumprod, const int ND, const int diffdim){
	int i, j;
	int *eff_cumprod, *subv, *diffsize, *diff_cumprod;

	subv = malloc(diffdim * sizeof(int));
	eff_cumprod = malloc(diffdim * sizeof(int));
	diffsize = malloc(diffdim * sizeof(int));
	diff_cumprod = malloc(diffdim * sizeof(int));
	for(i=0; i<diffdim; i++){
		eff_cumprod[i] = bCumprod[dmask[i]];
		diffsize[i] = (int)pbSize[dmask[i]];
	}
	diff_cumprod[0] = 1;
	for(i=0; i<diffdim-1; i++){
		diff_cumprod[i+1] = diff_cumprod[i] * diffsize[i];
	}
	for(i=0; i<ND; i++){
		ind_subv(i, diff_cumprod, diffdim, subv);
		weight[i] = 0;
		for(j=0; j<diffdim; j++){
			weight[i] += eff_cumprod[j] * subv[j];
		}
	}
	free(eff_cumprod);
	free(subv);
	free(diffsize);
	free(diff_cumprod);
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

mxArray* convert_table_to_sparse(const double *bT, const int *index, const int nzCounts, const int NB){
	mxArray  *spTable;
    int      i, *irs, *jcs;
    double   *sr;
    
	spTable = mxCreateSparse(NB, 1, nzCounts, mxREAL);
    sr  = mxGetPr(spTable);
    irs = mxGetIr(spTable);
    jcs = mxGetJc(spTable);

	jcs[0] = 0;
	jcs[1] = nzCounts;

	for(i=0; i<nzCounts; i++){
			sr[i] = bT[i];
			irs[i] = index[i];
    }
	return spTable;	
}

mxArray* convert_ill_table_to_sparse(const double *bigTable, const int *sequence, const int nzCounts, const int NB){
	mxArray *spTable;
	int     i, temp, *irs, *jcs, count=0;
	double  *sr;

	spTable = mxCreateSparse(NB, 1, nzCounts, mxREAL);
    sr  = mxGetPr(spTable);
    irs = mxGetIr(spTable);
    jcs = mxGetJc(spTable);

	jcs[0] = 0;
	jcs[1] = nzCounts;

	for(i=0; i<nzCounts; i++){
		irs[i] = sequence[count];
		count++;
		temp = sequence[count];
		sr[i] = bigTable[temp];
		count++;
	}
	return spTable;
}

void multiply_null_by_fuPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, NB, NS, siz_b, siz_s, ndim, nzCounts=0;
	int     *mask, *sx, *sy, *cpsy, *subs, *s, *cpsy2, *bir, *bjc;
	double  *pbDomain, *psDomain, *pbSize, *psSize, *spr, *bpr, value;
	mxArray *pTemp, *pTemp1;

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	siz_b = mxGetNumberOfElements(pTemp);
	pTemp = mxGetField(smallPot, 0, "domain");
	psDomain = mxGetPr(pTemp);
	siz_s = mxGetNumberOfElements(pTemp);

	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);
	pTemp = mxGetField(smallPot, 0, "sizes");
	psSize = mxGetPr(pTemp);

	NB = 1;
	for(i=0; i<siz_b; i++){
		NB *= (int)pbSize[i];
	}
	NS = 1;
	for(i=0; i<siz_s; i++){
		NS *= (int)psSize[i];
	}

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);

	pTemp1 = mxCreateSparse(NB, 1, NB, mxREAL);
	bpr = mxGetPr(pTemp1);
	bir = mxGetIr(pTemp1);
	bjc = mxGetJc(pTemp1);
	bjc[0] = 0;
	bjc[1] = NB;

	if(NS == 1){
		value = *spr;
		for(i=0; i<NB; i++){
			bpr[i] = value;
			bir[i] = i;
		}
		nzCounts = NB;
		pTemp = mxGetField(bigPot, 0, "T");
		if(pTemp)mxDestroyArray(pTemp);
		reset_nzmax(pTemp1, NB, nzCounts);
		mxSetField(bigPot, 0, "T", pTemp1);
		return;
	}

	if(NS == NB){
		for(i=0; i<NB; i++){
			if(spr[i] != 0){
				bpr[nzCounts] = spr[i];
				bir[nzCounts] = i;
				nzCounts++;
			}
		}
		pTemp = mxGetField(bigPot, 0, "T");
		if(pTemp)mxDestroyArray(pTemp);
		reset_nzmax(pTemp1, NB, nzCounts);
		mxSetField(bigPot, 0, "T", pTemp1);
		return;
	}

	mask = malloc(siz_s * sizeof(int));
	count = 0;
	for(i=0; i<siz_s; i++){
		for(j=0; j<siz_b; j++){
			if(psDomain[i] == pbDomain[j]){
				mask[count] = j;
				count++;
				break;
			}
		}
	}
	
	ndim = siz_b;
	sx = (int *)malloc(sizeof(int)*ndim);
	sy = (int *)malloc(sizeof(int)*ndim);
	for(i=0; i<ndim; i++){
		sx[i] = (int)pbSize[i];
		sy[i] = 1;
	}
	for(i=0; i<count; i++){
		sy[mask[i]] = sx[mask[i]];
	}

	s = (int *)malloc(sizeof(int)*ndim);
	*(cpsy = (int *)malloc(sizeof(int)*ndim)) = 1;
	subs =   (int *)malloc(sizeof(int)*ndim);
	cpsy2 =  (int *)malloc(sizeof(int)*ndim);
	for(i = 0; i < ndim; i++){
		subs[i] = 0;
		s[i] = sx[i] - 1;
	}
			
	for(i = 0; i < ndim-1; i++){
		cpsy[i+1] = cpsy[i]*sy[i]--;
		cpsy2[i] = cpsy[i]*sy[i];
	}
	cpsy2[ndim-1] = cpsy[ndim-1]*(--sy[ndim-1]);

	for(j=0; j<NB; j++){
		if(*spr != 0){
			bpr[nzCounts] = *spr;
			bir[nzCounts] = j;
			nzCounts++;
		}
		for(i = 0; i < ndim; i++){
			if(subs[i] == s[i]){
				subs[i] = 0;
				if(sy[i])
					spr -= cpsy2[i];
			}
			else{
				subs[i]++;
				if(sy[i])
					spr += cpsy[i];
				break;
			}
		}
	}

	pTemp = mxGetField(bigPot, 0, "T");
	if(pTemp)mxDestroyArray(pTemp);
	reset_nzmax(pTemp1, NB, nzCounts);
	mxSetField(bigPot, 0, "T", pTemp1);

	free(sx);
	free(sy);
	free(s);
	free(cpsy);
	free(subs);
	free(cpsy2);
    free(mask);
}

void multiply_null_by_spPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, count1, match, temp, bdim, sdim, diffdim, NB, NS, ND, NZB, NZS, bindex, sindex, nzCounts=0;
	int     *samemask, *diffmask, *sir, *sjc, *bCumprod, *sCumprod, *ssubv, *sequence, *weight;
	double  *bigTable, *pbDomain, *psDomain, *pbSize, *psSize, *spr;
	mxArray *pTemp, *pTemp1;

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	bdim = mxGetNumberOfElements(pTemp);
	pTemp = mxGetField(smallPot, 0, "domain");
	psDomain = mxGetPr(pTemp);
	sdim = mxGetNumberOfElements(pTemp);

	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);
	pTemp = mxGetField(smallPot, 0, "sizes");
	psSize = mxGetPr(pTemp);

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
	}
	NS = 1;
	for(i=0; i<sdim; i++){
		NS *= (int)psSize[i];
	}
	ND = NB / NS;

	if(ND == 1){
		pTemp = mxGetField(bigPot, 0, "T");
		if(pTemp)mxDestroyArray(pTemp);
		pTemp1 = mxGetField(smallPot, 0, "T");
		pTemp = mxDuplicateArray(pTemp1);
		mxSetField(bigPot, 0, "T", pTemp);
		return;
	}

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);
	sir = mxGetIr(pTemp);
	sjc = mxGetJc(pTemp);
	NZS = sjc[1];

	NZB = ND * NZS;

	diffdim = bdim - sdim;
	sequence = malloc(NZB * 2 * sizeof(int));
	bigTable = malloc(NZB * sizeof(double));
	samemask = malloc(sdim * sizeof(int));
	diffmask = malloc(diffdim * sizeof(int));
	bCumprod = malloc(bdim * sizeof(int));
	sCumprod = malloc(sdim * sizeof(int));
	weight = malloc(ND * sizeof(int));
	ssubv = malloc(sdim * sizeof(int));

	count = 0;
	count1 = 0;
	for(i=0; i<bdim; i++){
		match = 0;
		for(j=0; j<sdim; j++){
			if(pbDomain[i] == psDomain[j]){
				samemask[count] = i;
				match = 1;
				count++;
				break;
			}
		}
		if(match == 0){
			diffmask[count1] = i; 
			count1++;
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

	count = 0;
	compute_fixed_weight(weight, pbSize, diffmask, bCumprod, ND, diffdim);
	for(i=0; i<NZS; i++){
		sindex = sir[i];
		ind_subv(sindex, sCumprod, sdim, ssubv);
		temp = 0;
		for(j=0; j<sdim; j++){
			temp += ssubv[j] * bCumprod[samemask[j]];
		}
		for(j=0; j<ND; j++){
			bindex = weight[j] + temp;
			bigTable[nzCounts] = spr[i];
			sequence[count] = bindex;
			count++;
			sequence[count] = nzCounts;
			nzCounts++;
			count++;
		}
	}

	pTemp = mxGetField(bigPot, 0, "T");
	if(pTemp)mxDestroyArray(pTemp);
	qsort(sequence, nzCounts, sizeof(int) * 2, compare);
	pTemp = convert_ill_table_to_sparse(bigTable, sequence, nzCounts, NB);
	mxSetField(bigPot, 0, "T", pTemp);

	free(sequence); 
	free(bigTable);
	free(samemask);
	free(diffmask);
	free(bCumprod);
	free(sCumprod);
	free(weight);
	free(ssubv);
}

void multiply_spPot_by_fuPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, bdim, sdim, NB, NZB, bindex, sindex, nzCounts=0;
	int     *mask, *index, *bir, *bjc, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double  *bigTable, *pbDomain, *psDomain, *pbSize, *psSize, *bpr, *spr, value;
	mxArray *pTemp;

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	bdim = mxGetNumberOfElements(pTemp);
	pTemp = mxGetField(smallPot, 0, "domain");
	psDomain = mxGetPr(pTemp);
	sdim = mxGetNumberOfElements(pTemp);

	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);
	pTemp = mxGetField(smallPot, 0, "sizes");
	psSize = mxGetPr(pTemp);

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
	}

	pTemp = mxGetField(bigPot, 0, "T");
	bpr = mxGetPr(pTemp);
	bir = mxGetIr(pTemp);
	bjc = mxGetJc(pTemp);
	NZB = bjc[1];

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);

	bigTable = malloc(NZB * sizeof(double));
	index = malloc(NZB * sizeof(double));
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
		value = spr[sindex];
		if(value != 0){
			bigTable[nzCounts] = bpr[i] * value;
			index[nzCounts] = bindex;
			nzCounts++;
		}
	}

	pTemp = mxGetField(bigPot, 0, "T");
	if(pTemp)mxDestroyArray(pTemp);
	pTemp = convert_table_to_sparse(bigTable, index, nzCounts, NB);
	mxSetField(bigPot, 0, "T", pTemp);

	free(bigTable);
	free(index);
	free(mask);
	free(bCumprod);
	free(sCumprod);
	free(bsubv);
	free(ssubv);
}

void multiply_spPot_by_spPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, bdim, sdim, NB, NZB, NZS, position, bindex, sindex, nzCounts=0;
	int     *mask, *index, *result, *bir, *sir, *bjc, *sjc, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double  *bigTable, *pbDomain, *psDomain, *pbSize, *psSize, *bpr, *spr, value;
	mxArray *pTemp;

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	bdim = mxGetNumberOfElements(pTemp);
	pTemp = mxGetField(smallPot, 0, "domain");
	psDomain = mxGetPr(pTemp);
	sdim = mxGetNumberOfElements(pTemp);

	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);
	pTemp = mxGetField(smallPot, 0, "sizes");
	psSize = mxGetPr(pTemp);

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
	}

	pTemp = mxGetField(bigPot, 0, "T");
	bpr = mxGetPr(pTemp);
	bir = mxGetIr(pTemp);
	bjc = mxGetJc(pTemp);
	NZB = bjc[1];

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);
	sir = mxGetIr(pTemp);
	sjc = mxGetJc(pTemp);
	NZS = sjc[1];

	bigTable = malloc(NZB * sizeof(double));
	index = malloc(NZB * sizeof(double));
	mask = malloc(sdim * sizeof(int));
	bCumprod = malloc(bdim * sizeof(int));
	sCumprod = malloc(sdim * sizeof(int));
	bsubv = malloc(bdim * sizeof(int));
	ssubv = malloc(sdim * sizeof(int));

	for(i=0; i<NZB; i++){
		bigTable[i] = 0;
	}
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
		value = bpr[i];
		bindex = bir[i];
		ind_subv(bindex, bCumprod, bdim, bsubv);
		for(j=0; j<sdim; j++){
			ssubv[j] = bsubv[mask[j]];
		}
		sindex = subv_ind(sdim, sCumprod, ssubv);
		result = (int *) bsearch(&sindex, sir, NZS, sizeof(int), compare);
		if(result){
			position = result - sir;
			value *= spr[position];
			bigTable[nzCounts] = value;
			index[nzCounts] = bindex;
			nzCounts++;
		}
	}

	pTemp = mxGetField(bigPot, 0, "T");
	if(pTemp)mxDestroyArray(pTemp);
	pTemp = convert_table_to_sparse(bigTable, index, nzCounts, NB);
	mxSetField(bigPot, 0, "T", pTemp);

	free(bigTable);
	free(index);
	free(mask);
	free(bCumprod);
	free(sCumprod);
	free(bsubv);
	free(ssubv);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int     i, j, c, loop, nNodes, nCliques, ndomain, dims[2];
	double  *pClqs, *pr, *pt, *pSize;
	mxArray *pTemp, *pTemp1, *pStruct, *pCliques, *pBigpot, *pSmallpot;
	const char *field_names[] = {"domain", "T", "sizes"};

	nNodes = mxGetNumberOfElements(prhs[1]);
	pCliques = mxGetField(prhs[0], 0, "cliques");
	nCliques = mxGetNumberOfElements(pCliques);
	pTemp = mxGetField(prhs[0], 0, "eff_node_sizes");
	pSize = mxGetPr(pTemp);

	plhs[0] = mxCreateCellArray(1, &nCliques);
    for(i=0; i<nCliques; i++){
        pStruct = mxCreateStructMatrix(1, 1, 3, field_names);
		mxSetCell(plhs[0], i, pStruct);
		pTemp = mxGetCell(pCliques, i);
		ndomain = mxGetNumberOfElements(pTemp);
		pt = mxGetPr(pTemp);
		pTemp1 = mxDuplicateArray(pTemp);
		mxSetField(pStruct, 0, "domain", pTemp1);
		
		pTemp = mxCreateDoubleMatrix(1, ndomain, mxREAL);
		mxSetField(pStruct, 0, "sizes", pTemp);
		pr = mxGetPr(pTemp);
        for(j=0; j<ndomain; j++){
            pr[j] = pSize[(int)pt[j]-1];
        }
    }

	pClqs = mxGetPr(prhs[1]);
	for(loop=0; loop<nNodes; loop++){
		c = (int)pClqs[loop] - 1;
		pSmallpot = mxGetCell(prhs[2], loop);
		pTemp = mxGetField(pSmallpot, 0, "T");
		pBigpot = mxGetCell(plhs[0], c);
		pTemp1 = mxGetField(pBigpot, 0, "T");
		if(pTemp1){
			if(mxIsSparse(pTemp))
				multiply_spPot_by_spPot(pBigpot, pSmallpot);
			else multiply_spPot_by_fuPot(pBigpot, pSmallpot);
		}
		else{
			if(mxIsSparse(pTemp))
				multiply_null_by_spPot(pBigpot, pSmallpot);
			else multiply_null_by_fuPot(pBigpot, pSmallpot);
		}		
	}

	dims[0] = nCliques;
	dims[1] = nCliques;
	plhs[1] = mxCreateCellArray(2, dims);
}


