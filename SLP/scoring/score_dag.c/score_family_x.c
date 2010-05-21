//
// File: score_family_x.c
//
// MATLAB:
// function score = score_family_x(data,sz)
//
// DESCRIPTION:
// Calculates Bayesian score for the family of tabular nodes
// (as in log_marg_prob_node.m):
//      data        [nsz x ndata] array
//                  IMPORTANT:                    <--- !!!!
//                  - parents go first                             
//                  - values of the nodes are in the range 0..sz[i] 
//      sz          sizes of the nodes in the family
//      ndata       # of observations
//      nsz         # of nodes in the family
// with default parameters:  
//      params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 }
//                   ...and 'dirichlet_type','BDeu'
//      scoring_fn = 'bayesian';
//      etc...
//
// EXAMPLE:
//   score_family_x([1 2 1; 1 2 2; 1 1 1],[2 2 2]);  // ans = -2.0794
//
// by Darima <darrimma@yahoo.com>, 27/12/2005 
//

#include "mex.h"
#include "score_x.h"

#define	IN_DATA	    prhs[0]
#define IN_SZ       prhs[1]
#define	OUT_SCORE	plhs[0]
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )   
{ 
    int nsz, ndata;
    double  *sz, *data, *score;
	int i,j;
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
        mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments."); 
    } 

    /* Assign input arguments */ 
    
    data   = mxGetPr(IN_DATA);
    sz     = mxGetPr(IN_SZ);
	ndata  = mxGetN(IN_DATA); 
    nsz    = MAX(mxGetM(IN_SZ),mxGetN(IN_SZ));   

    /* Create return argument */ 
    
    OUT_SCORE = mxCreateDoubleScalar(mxREAL);
    score = mxGetPr(OUT_SCORE);

    /* Do the actual computations in a subroutine */
    
    score_family_x(data,sz,nsz,ndata,score);
    return;
}

