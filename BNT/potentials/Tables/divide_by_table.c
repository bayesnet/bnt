/* divide_by_table.c  ../potential/tables  */


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

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int     i, j, count, NB, NS, siz_b, siz_s, ndim, temp;
	int     *mask, *sx, *sy, *cpsy, *subs, *s, *cpsy2;
	double  *pbDomain, *psDomain, *sp, *zp, *bs, value;

	plhs[0] = mxDuplicateArray(prhs[0]);
	zp = mxGetPr(plhs[0]);

	siz_b = mxGetNumberOfElements(prhs[1]);
	siz_s = mxGetNumberOfElements(prhs[4]);
	pbDomain = mxGetPr(prhs[1]);
	psDomain = mxGetPr(prhs[4]);

	NB = mxGetNumberOfElements(prhs[0]);
	NS = mxGetNumberOfElements(prhs[3]);
	sp = mxGetPr(prhs[3]);

	bs = mxGetPr(prhs[2]);

	if(NS == 1){
		value = *sp;
		if(value == 0) value = 1;
		for(i=0; i<NB; i++){
			zp[i] /= value;
		}
		return;
	}

	if(NS == NB){
		for(i=0; i<NB; i++){
			value = sp[i];
			if(value == 0) value = 1;
			zp[i] /= value;
		}
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
		sx[i] = (int)bs[i];
		sy[i] = 1;
	}
	for(i=0; i<count; i++){
		temp = mask[i];
		sy[temp] = sx[temp];
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
		value = *sp;
		if(value == 0) value = 1;
		*zp++ /= value;
		for(i = 0; i < ndim; i++){
			if(subs[i] == s[i]){
				subs[i] = 0;
				if(sy[i])
					sp -= cpsy2[i];
			}
			else{
				subs[i]++;
				if(sy[i])
					sp += cpsy[i];
				break;
			}
		}
	}
	free(sx);
	free(sy);
	free(s);
	free(cpsy);
	free(subs);
	free(cpsy2);
    free(mask);
}

