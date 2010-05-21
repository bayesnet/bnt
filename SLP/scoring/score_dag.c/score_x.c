//
// File: score_x.c
//
//
// by Darima <darrimma@yahoo.com>, 27/12/2005 
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "score_x.h"

void score_dag_x(double *data,double *sz,double *dag,
    int ndata,int nsz,double* score)
//
// GOTO: score_dag_x.c
// by Darima <darrimma@yahoo.com>, 27/12/2005 
//    
{
	double *fdata,*fsz,*fscore;
    int fnsz,*family; 
	int i,j,k; 
        
	fdata = (double*)malloc(nsz*ndata*sizeof(double));
	fsz = (double*)malloc(nsz*sizeof(double));
    fscore = (double*)malloc(sizeof(double));
    family = (int*)malloc(nsz*sizeof(int));
	*score = 0;      
   
	for (j=0;j<nsz;j++) {
        
		// find parents of current node
		        
        fnsz = 0;
		for (i=0;i<nsz;i++)
			if (dag[i+j*nsz]==1) {family[fnsz] = i; fnsz++;}
        family[fnsz] = j; fnsz++;
        
        // initialize data for current family
        
		for (i=0;i<fnsz;i++)
            for (k=0;k<ndata;k++) {
                fdata[i+k*fnsz] = data[family[i]+k*nsz];
                fsz[i] = sz[family[i]];
            }
        
		// calculater score for the family
        
        score_family_x(fdata,fsz,fnsz,ndata,fscore);
        *score += *fscore;
	}

	free(fdata); free(fsz); free(fscore); free(family);
	return;
}

void score_family_x(double *data,double *sz,
    int nsz,int ndata,double* score)
//
// GOTO: score_family_x.c
// by Darima <darrimma@yahoo.com>, 27/12/2005 
//
{ 
    int self_sz, ps_sz;
    double prior1, prior2, *count;
    int i,j,k; double idx,tsz,N_ij;
        
    // self_sz - number of values of current node
    // ps_sz   - number of configurations of parents    
    
    self_sz = sz[nsz-1]; 
    ps_sz = 1; for (i=0;i<nsz-1;i++){ps_sz *= sz[i];}
       
    // calculate counts
    //
    // example: two binary parents for binary node
    //  idx |  p1 p2 |   0   1   <-- values of the node
    //  ---------------------
    //   0  |  0  0  |   1   0
    //   1  |  1  0  |   0  83
    //   2  |  0  1  |   3   0   <-- counts
    //   3  |  1  1  |  13   0
    //  ^------------------------ configurations of parents
    // count = [ps_sz x self_sz] array
        
    count = (double*)malloc(ps_sz*self_sz*sizeof(double));
    for (i=0;i<self_sz*ps_sz;i++) { count[i] = 0; }

    for (j=0;j<ndata;j++) {
        
        // for every data case:
        // 1) calculate parent configuration index: p1_val,p2_val -> idx
        // 2) increment corresponding count: count[idx,node_val]++
        //
        // example:
        //             case1  case2  etc
        // -----------------------------
        // p1   |      0      1      ...
        // p2   |      0      0      ...
        // node |      1      1      ...
        // 
        // notice: array[i,j,k] -> array[i+j*ni+k*ni*nj]
        //
        // example 1: data[i,j] -> data[i+j*nsz]
        //
        // example 2: consider family with 5 parents
        // - [p1,p2,..,p5] -> IDX5 
        //   IDX5 = p1 + p2*p1_sz + p3*p1_sz*p2_sz + ... +
        //          p5*p1_sz*p2_sz*p3_sz*p4_sz
        // - [p1,p2,..,p5,node_val] -> IDXFAM
        //   IDXFAM = IDX5 + node_val*ps_sz
        //   count[IDX5,node_val] ->count[IDXFAM]
        
        idx = data[0+j*nsz]; tsz = 1;
        for (i=1;i<nsz;i++) {tsz *= sz[i-1]; idx += data[i+j*nsz]*tsz;}
        count[(int)idx]++;
    }
       
    // calculate priors (BDeu,1)
    
    prior1 = 1/(double)(self_sz*ps_sz);
    prior2 = 1/(double)(ps_sz);
       
    // calculate score
    //
    // LL = log[  PROD_j gamma(alpha_ij)/gamma(alpha_ij + N_ij)  *
	//            PROD_k gamma(alpha_ijk + N_ijk)/gamma(alpha_ijk)  ] =
    //    = SUM_j {  gammaln(alpha_ij)-gammaln(alpha_ij + N_ij) +
    //         SUM_k [ gammaln(alpha_ijk + N_ijk)-gammaln(alpha_ijk) ]  }
    //
    // alpha_ijk  <- prior1
    // alpha_ij   <- prior2
    // N_ijk      <- count[j,k]

    *score = 0;
    for (j=0;j<ps_sz;j++) {    
        N_ij = .0; 
        for (k=0;k<self_sz;k++) {
            N_ij += count[j+k*ps_sz];
            *score += gammaln(prior1+count[j+k*ps_sz])-gammaln(prior1);
        }
        *score += gammaln(prior2)-gammaln(prior2 + N_ij);
    }

    free(count);
    return;
}

double gammaln(double xx)
//
// Returns the value ln[Gamma(xx)] for xx > 0.
// from "Numerical Recipes in C"
//
{
    double x,y,tmp,ser;
    static double cof[6]={
        76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
