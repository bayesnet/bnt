/* C mex for distribute_evidence.c in @jtree_sparse_inf_engine directory*/
/* File enter_evidence.m in directory @jtree_sparse_inf_engine call it  */

/*********************************************/
/* distribute_evidence has 3 input & 2 output*/
/* engine                                    */
/* clpot                                     */
/* seppot                                    */
/*                                           */
/* clpot                                     */
/* seppot                                    */
/*********************************************/

#include "mex.h"

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

mxArray* convert_table_to_sparse(const double *bT, const int *index, const int nzCounts, const int N){
	mxArray  *spTable;
    int      i, *irs, *jcs;
    double   *sr;
    
	spTable = mxCreateSparse(N, 1, nzCounts, mxREAL);
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

mxArray* convert_ill_table_to_sparse(const double *Table, const int *sequence, const int nzCounts, const int N){
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

	if(sdim == 0){
		for(i=0; i<NZB; i++){
			bpr[i] *= *spr;
		}	
		return;
	}

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

void marginal_spPot_to_spPot(const mxArray *bigPot, mxArray *smallPot, const int maximize){
	int     i, j, count, bdim, sdim, NB, NS, NZB, position, bindex, sindex, nzCounts=0;
	int     *mask, *sequence, *result, *bir, *bjc, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double  *sTable, *pbDomain, *psDomain, *pbSize, *psSize, *bpr, *spr;
	mxArray *pTemp;

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	bdim = mxGetNumberOfElements(pTemp);
	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);

	pTemp = mxGetField(smallPot, 0, "domain");
	psDomain = mxGetPr(pTemp);
	sdim = mxGetNumberOfElements(pTemp);
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

	pTemp = mxGetField(bigPot, 0, "T");
	bpr = mxGetPr(pTemp);
	bir = mxGetIr(pTemp);
	bjc = mxGetJc(pTemp);
	NZB = bjc[1];

	if(sdim == 0){
		pTemp = mxGetField(smallPot, 0, "T");
		spr = mxGetPr(pTemp);
		*spr = 0;
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

	for(i=0; i<NZB; i++){
		sTable[i] = 0;
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
	
	pTemp = mxGetField(smallPot, 0, "T");
	if(pTemp)mxDestroyArray(pTemp);
	qsort(sequence, nzCounts, sizeof(int) * 2, compare);
	pTemp = convert_ill_table_to_sparse(sTable, sequence, nzCounts, NS);
	mxSetField(smallPot, 0, "T", pTemp);

	free(sTable);
	free(sequence);
	free(mask);
	free(bCumprod);
	free(sCumprod);
	free(bsubv);
	free(ssubv);
}

void divide_null_by_spPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, count1, match, temp, bdim, sdim, diffdim, NB, NS, ND, NZB, NZS, bindex, sindex;
	int     *samemask, *diffmask, *rir, *rjc, *sir, *sjc, *bCumprod, *sCumprod, *ssubv, *weight;
	double  *pbDomain, *psDomain, *pbSize, *psSize, *rpr, *spr, value;
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

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);
	sir = mxGetIr(pTemp);
	sjc = mxGetJc(pTemp);
	NZS = sjc[1];

	if(sdim == 0){
		pTemp = mxCreateSparse(NB, 1, NB, mxREAL);
		mxSetField(bigPot, 0, "T", pTemp);
		rpr = mxGetPr(pTemp);
		rir = mxGetIr(pTemp);
		rjc = mxGetJc(pTemp);
		rjc[0] = 0;
		rjc[1] = NB;
		value = *spr;
		if(value == 0) value = 1;
		for(i=0; i<NB; i++){
			rpr[i] = 1 / value;
			rir[i] = i;
		}	
		return;
	}

	NS = 1;
	for(i=0; i<sdim; i++){
		NS *= (int)psSize[i];
	}
	ND = NB / NS;


	pTemp = mxCreateSparse(NB, 1, NB, mxREAL);
	rpr = mxGetPr(pTemp);
	rir = mxGetIr(pTemp);
	rjc = mxGetJc(pTemp);
	rjc[0] = 0;
	rjc[1] = NB;
	for(i=0; i<NB; i++){
		rpr[i] = 1;
		rir[i] = i;
	}

	NZB = ND * NZS;

	diffdim = bdim - sdim;
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
			rpr[bindex] = 1 / (spr[i]);
		}
	}

	pTemp1 = mxGetField(bigPot, 0, "T");
	if(pTemp1)mxDestroyArray(pTemp1);
	mxSetField(bigPot, 0, "T", pTemp);

	free(samemask);
	free(diffmask);
	free(bCumprod);
	free(sCumprod);
	free(weight);
	free(ssubv);
}

void divide_spPot_by_spPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, bdim, sdim, NB, NZB, NZS, position, bindex, sindex;
	int     *mask, *result, *bir, *sir, *bjc, *sjc, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double  *pbDomain, *psDomain, *pbSize, *psSize, *bpr, *spr, value;
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

	pTemp1 = mxGetField(bigPot, 0, "T");
	bpr = mxGetPr(pTemp1);
	bir = mxGetIr(pTemp1);
	bjc = mxGetJc(pTemp1);
	NZB = bjc[1];

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);
	sir = mxGetIr(pTemp);
	sjc = mxGetJc(pTemp);
	NZS = sjc[1];

	if(sdim == 0){
		value = *spr;
		if(value == 0)value = 1;
		for(i=0; i<NZB; i++){
			bpr[i] /= value;
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
			bpr[i] /= spr[position];
		}
	}

	free(mask);
	free(bCumprod);
	free(sCumprod);
	free(bsubv);
	free(ssubv);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int     i, j, loop, loops, nCliques, temp, count, parent, child, maximize, *distribute_order;
	double  *pr, *pr1;
	mxArray *pTemp, *pPreCh, *pClpot, *pSeppot;

	pTemp = mxGetField(prhs[0], 0, "cliques");
	nCliques = mxGetNumberOfElements(pTemp);
	loops = nCliques - 1;
	pTemp = mxGetField(prhs[0], 0, "maximize");
	maximize = (int)mxGetScalar(pTemp);

	distribute_order = malloc(2 * loops * sizeof(int));
	pTemp = mxGetField(prhs[0], 0, "preorder");
	pr = mxGetPr(pTemp);
	pPreCh = mxGetField(prhs[0], 0, "preorder_children");
	count = 0;
	for(i=0; i<nCliques; i++){
		temp = (int)pr[i] - 1;
		pTemp = mxGetCell(pPreCh, temp);
		pr1 = mxGetPr(pTemp);
		loop = mxGetNumberOfElements(pTemp);
		for(j=0; j<loop; j++){
			distribute_order[count] = temp;
			distribute_order[count + loops] = (int)pr1[j] - 1;
			count++;
		}
	}

	plhs[0] = mxDuplicateArray(prhs[1]);
	plhs[1] = mxDuplicateArray(prhs[2]);

	for(loop=0; loop<loops; loop++){
		parent = distribute_order[loop];
		child  = distribute_order[loop+loops];
		i = nCliques * child + parent;
		pClpot = mxGetCell(plhs[0], child);
		pTemp = mxGetField(pClpot, 0, "T");
		pSeppot = mxGetCell(plhs[1], i);
		if(pTemp)
			divide_spPot_by_spPot(pClpot, pSeppot);
		else divide_null_by_spPot(pClpot, pSeppot);

		pClpot = mxGetCell(plhs[0], parent);
		marginal_spPot_to_spPot(pClpot, pSeppot, maximize);
		mxSetCell(plhs[1], i, pSeppot);

		pClpot = mxGetCell(plhs[0], child);
		multiply_spPot_by_spPot(pClpot, pSeppot); 
	}
	free(distribute_order);
}
