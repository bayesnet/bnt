/* marg_table.c  ../potential/tables     */


/******************************************/
/* 5 input & 1 output                     */
/* Big table                              */
/* Big domain                             */
/* Big sizes                              */
/* onto                                   */
/* maximize, if missed, maximize=0        */
/*                                        */
/* small table                            */
/******************************************/

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int     i, j, count, NB, NS, siz_b, siz_s, ndim, temp, maximize;
	int     *mask, *sx, *sy, *cpsy, *subs, *s, *cpsy2, *ssize;
	double  *pb, *ps, *bp, *sp, *pbd;


	siz_b = mxGetNumberOfElements(prhs[1]);
	siz_s = mxGetNumberOfElements(prhs[3]);
	pb = mxGetPr(prhs[1]);
	ps = mxGetPr(prhs[3]);

	NB = mxGetNumberOfElements(prhs[0]);
	bp = mxGetPr(prhs[0]);

	pbd = mxGetPr(prhs[2]);

	if(nrhs < 5) maximize = 0;
	else maximize = (int)mxGetScalar(prhs[4]);

	if(siz_s == 0){
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		sp = mxGetPr(plhs[0]);
		if(maximize){
			for(i=0; i<NB; i++){
				*sp = (*sp < bp[i])? bp[i] : *sp;
			}
		}
		else{
			for(i=0; i<NB; i++){
				*sp += bp[i];
			}
		}
		return;
	}

	mask = malloc(siz_s * sizeof(int));
	ssize = malloc(siz_s * sizeof(int));
	count = 0;
	for(i=0; i<siz_s; i++){
		for(j=0; j<siz_b; j++){
			if(ps[i] == pb[j]){
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
		sx[i] = (int)pbd[i];
		sy[i] = 1;
	}
	for(i=0; i<siz_s; i++){
		temp = mask[i];
		sy[temp] = sx[temp];
		ssize[i] = sx[temp];
	}

	NS = 1;
	for(i=0; i<ndim; i++){
		NS *= sy[i];
	}

	plhs[0] = mxCreateNumericArray(siz_s, ssize, mxDOUBLE_CLASS, mxREAL);
	sp = mxGetPr(plhs[0]);

	if(NS == 1){
		if(maximize){
			for(i=0; i<NB; i++){
				*sp = (*sp < bp[i])? bp[i] : *sp;
			}
		}
		else{
			for(i=0; i<NB; i++){
				*sp += bp[i];
			}
		}
		free(mask);
		free(sx);
		free(sy);
		free(ssize);
		return;
	}

	if(NS == NB){
		for(i=0; i<NB; i++) *sp++ = *bp++;
		free(mask);
		free(sx);
		free(sy);
		free(ssize);
		return;
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

	if(maximize){
		for(j=0; j<NB; j++){
			*sp = (*sp < *bp)? *bp : *sp;
			bp++;
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
	}
	else{
		for(j=0; j<NB; j++){
			*sp += *bp++;
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
	}

	free(sx);
	free(sy);
	free(s);
	free(cpsy);
	free(subs);
	free(cpsy2);
    free(mask);
	free(ssize);
}

