/* C mex for collect_evidence.c in @jtree_sparse_inf_engine directory */
/* File enter_evidence.m in directory @jtree_sparse_inf_engine call it*/

/******************************************/
/* collect_evidence has 3 input & 2 output*/
/* engine                                 */
/* clpot                                  */
/* seppot                                 */
/*                                        */
/* clpot                                  */
/* seppot                                 */
/******************************************/

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

void multiply_null_by_spPot(mxArray *bigPot, const mxArray *smallPot){
	int     i, j, count, count1, match, temp, bdim, sdim, diffdim, NB, NS, ND, NZB, NZS, bindex, sindex, nzCounts=0;
	int     *samemask, *diffmask, *sir, *sjc, *bCumprod, *sCumprod, *ssubv, *sequence, *weight;
	double  *bigTable, *pbDomain, *psDomain, *pbSize, *psSize, *spr, *bpr;
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

	pTemp = mxGetField(smallPot, 0, "T");
	spr = mxGetPr(pTemp);
	sir = mxGetIr(pTemp);
	sjc = mxGetJc(pTemp);
	NZS = sjc[1];

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
	}

	if(sdim == 0){
		pTemp = mxCreateSparse(NB, 1, NB, mxREAL);
		mxSetField(bigPot, 0, "T", pTemp);
		bpr = mxGetPr(pTemp);
		sir = mxGetIr(pTemp);
		sjc = mxGetJc(pTemp);
		sjc[0] = 0;
		sjc[1] = NB;
		for(i=0; i<NB; i++){
			bpr[i] = *spr;
			sir[i] = i;
		}	
		return;
	}

	NS = 1;
	for(i=0; i<sdim; i++){
		NS *= (int)psSize[i];
	}
	ND = NB / NS;

	if(ND == 1){
		pTemp1 = mxGetField(smallPot, 0, "T");
		pTemp = mxDuplicateArray(pTemp1);
		mxSetField(bigPot, 0, "T", pTemp);
		return;
	}


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

mxArray* marginal_null_to_spPot(const mxArray *bigPot, const mxArray *sDomain, const int maximize){
	int     i, j, count, bdim, sdim, NB, NS, ND;
	int     *mask, *sir, *sjc;
	double  *pbDomain, *psDomain, *pbSize, *psSize, *spr;
	mxArray *pTemp, *smallPot;
	const char *field_names[] = {"domain", "T", "sizes"};

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	bdim = mxGetNumberOfElements(pTemp);
	psDomain = mxGetPr(sDomain);
	sdim = mxGetNumberOfElements(sDomain);
	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);

	smallPot = mxCreateStructMatrix(1, 1, 3, field_names);
	pTemp = mxDuplicateArray(sDomain);
	mxSetField(smallPot, 0, "domain", pTemp);

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
	}

	if(sdim == 0){
		pTemp = mxCreateSparse(1, 1, 1, mxREAL);
		mxSetField(smallPot, 0, "T", pTemp);
		spr = mxGetPr(pTemp);
		sir = mxGetIr(pTemp);
		sjc = mxGetJc(pTemp);
		*spr = 0;
		*sir = 0;
		sjc[0] = 0;
		sjc[1] = 1;
		if(maximize) *spr = 1;
		else *spr = NB;

		pTemp = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(pTemp) = 1;
		mxSetField(smallPot, 0, "sizes", pTemp);
		return smallPot;
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
	pTemp = mxCreateDoubleMatrix(1, count, mxREAL);
	psSize = mxGetPr(pTemp);
	NS = 1;
	for(i=0; i<count; i++){
		psSize[i] = pbSize[mask[i]];
		NS *= (int)psSize[i];
	}
	mxSetField(smallPot, 0, "sizes", pTemp);

	ND = NB / NS;

	pTemp = mxCreateSparse(NS, 1, NS, mxREAL);
	mxSetField(smallPot, 0, "T", pTemp);
	spr = mxGetPr(pTemp);
	sir = mxGetIr(pTemp);
	sjc = mxGetJc(pTemp);
	if(maximize){
		for(i=0; i<NS; i++){
			spr[i] = 1;
			sir[i] = i;
		}
	}
	else{
		for(i=0; i<NS; i++){
			spr[i] = ND;
			sir[i] = i;
		}
	}
	sjc[0] = 0;
	sjc[1] = NS;

	free(mask);
	return smallPot;
}

mxArray* marginal_spPot_to_spPot(const mxArray *bigPot, const mxArray *sDomain, const int maximize){
	int     i, j, count, bdim, sdim, NB, NS, NZB, position, bindex, sindex, nzCounts=0;
	int     *mask, *sequence, *result, *bir, *bjc, *bCumprod, *sCumprod, *bsubv, *ssubv;
	double  *sTable, *pbDomain, *psDomain, *pbSize, *psSize, *bpr, *spr;
	mxArray *pTemp, *smallPot;
	const char *field_names[] = {"domain", "T", "sizes"};

	pTemp = mxGetField(bigPot, 0, "domain");
	pbDomain = mxGetPr(pTemp);
	bdim = mxGetNumberOfElements(pTemp);
	psDomain = mxGetPr(sDomain);
	sdim = mxGetNumberOfElements(sDomain);
	pTemp = mxGetField(bigPot, 0, "sizes");
	pbSize = mxGetPr(pTemp);

	pTemp = mxGetField(bigPot, 0, "T");
	bpr = mxGetPr(pTemp);
	bir = mxGetIr(pTemp);
	bjc = mxGetJc(pTemp);
	NZB = bjc[1];

	smallPot = mxCreateStructMatrix(1, 1, 3, field_names);
	pTemp = mxDuplicateArray(sDomain);
	mxSetField(smallPot, 0, "domain", pTemp);

	if(sdim == 0){
		pTemp = mxCreateSparse(1, 1, 1, mxREAL);
		mxSetField(smallPot, 0, "T", pTemp);
		spr = mxGetPr(pTemp);
		bir = mxGetIr(pTemp);
		bjc = mxGetJc(pTemp);
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

		pTemp = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(pTemp) = 1;
		mxSetField(smallPot, 0, "sizes", pTemp);
		return smallPot;
	}

	NB = 1;
	for(i=0; i<bdim; i++){
		NB *= (int)pbSize[i];
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
	pTemp = mxCreateDoubleMatrix(1, count, mxREAL);
	psSize = mxGetPr(pTemp);
	NS = 1;
	for(i=0; i<count; i++){
		psSize[i] = pbSize[mask[i]];
		NS *= (int)psSize[i];
	}
	mxSetField(smallPot, 0, "sizes", pTemp);


	sTable = malloc(NZB * sizeof(double));
	sequence = malloc(NZB * 2 * sizeof(double));
	bCumprod = malloc(bdim * sizeof(int));
	sCumprod = malloc(sdim * sizeof(int));
	bsubv = malloc(bdim * sizeof(int));
	ssubv = malloc(sdim * sizeof(int));

	for(i=0; i<NZB; i++)sTable[i] = 0;
	
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

	return smallPot;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int     i, n, p, np, pn, loop, loops, nCliques, temp, maximize;
	int     *collect_order;
	double  *pr, *pr1;
	mxArray *pTemp, *pTemp1, *pPostP, *pClpot, *pSeppot, *pSeparator;

	pTemp = mxGetField(prhs[0], 0, "cliques");
	nCliques = mxGetNumberOfElements(pTemp);
	loops = nCliques - 1;
	pTemp = mxGetField(prhs[0], 0, "maximize");
	maximize = (int)mxGetScalar(pTemp);
	pSeparator = mxGetField(prhs[0], 0, "separator");

	collect_order = malloc(2 * loops * sizeof(int));

	pTemp = mxGetField(prhs[0], 0, "postorder");
	pr = mxGetPr(pTemp);
	pPostP = mxGetField(prhs[0], 0, "postorder_parents");
	for(i=0; i<loops; i++){
		temp = (int)pr[i] - 1;
		pTemp = mxGetCell(pPostP, temp);
		pr1 = mxGetPr(pTemp);
		collect_order[i] = (int)pr1[0] - 1;
		collect_order[i+loops] = temp;
	}

	plhs[0] = mxDuplicateArray(prhs[1]);
	plhs[1] = mxDuplicateArray(prhs[2]);

	for(loop=0; loop<loops; loop++){
		p = collect_order[loop];
		n = collect_order[loop+loops];
		np = p * nCliques + n;
		pn = n * nCliques + p;
		pClpot = mxGetCell(plhs[0], n);
		pTemp1 = mxGetField(pClpot, 0, "T");
		pTemp = mxGetCell(pSeparator, pn);
		if(pTemp1)
			pSeppot = marginal_spPot_to_spPot(pClpot, pTemp, maximize);
		else pSeppot = marginal_null_to_spPot(pClpot, pTemp, maximize);
		mxSetCell(plhs[1], pn, pSeppot);

		pClpot = mxGetCell(plhs[0], p);
		pTemp1 = mxGetField(pClpot, 0, "T");
		if(pTemp1)
			multiply_spPot_by_spPot(pClpot, pSeppot);
		else multiply_null_by_spPot(pClpot, pSeppot);
	}
	free(collect_order);
}
	





