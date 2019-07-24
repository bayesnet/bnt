/* triangulate.c written by Ilya Shpitser  */

#include <stdlib.h>

#ifdef UNIX
#include "matlab.h"
#endif

#include "matrix.h"
#include "mex.h"

#include "elim.h"
#include "map.h"
#include "misc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  int dims [2];
  int i, j, k, m, n;
  long index;
  double * G_pr;
  double * stage_pr;
  double * answer_G_pr, * fill_ins_pr;
  double * matlab_clique_pr;
  mxArray * matlab_clique;
  Elimination e;
  float ** adj_mat;
  int ** order = (int **) NULL;
  Iterator iter, iter2;
  word w, w2;
  int ** fill_ins;
  Map cliques;
  Map clique;
  mxArray * fill_ins_mat;
  int * nodes;
  mxArray * full;

// (original)  full = mlfFull((mxArray *) prhs[0]);
  full = (mxArray *) mlfFull((mxArray *) prhs[0]);  // added typecasting
  /* Obtain graph matrix information. */
  m = mxGetM(full);
  n = mxGetN(full);
  G_pr = mxGetPr(full);

  if(n < 1 || m < 1){
    return;
  }

  /* Allocate and populate the log weight adjacency matrix corresponding
     to the input graph. */
  adj_mat = (float **) malloc(sizeof(float *) * m);
  adj_mat[0] = (float *) malloc(sizeof(float) * m * n);
  for(i = 1; i < m; i++){
    adj_mat[i] = adj_mat[i - 1] + n;
  }
  /* We no longer have log weight info, but we have a (total) ordering on
     the nodes already, so we do not need this information. */
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      index = j * m + i;
      if(G_pr[index] > 0){
        adj_mat[i][j] = 1;
      } else {
        adj_mat[i][j] = 0;
      }
    }
  }

  /* Convert the total elimination ordering into a partial order argument
     for the elimination routine.  The elimination routine's purpose in this
     mode of operation is to return cliques and fill-in edges. */
  if(nrhs > 1){
    order = (int **) malloc(sizeof(int *) * m);
    order[0] = (int *) malloc(sizeof(int) * m * n);
    for(i = 1; i < m; i++){
      order[i] = order[i - 1] + n;
    }
    for(i = 0; i < m; i++){
      for(j = 0; j < n; j++){
        order[i][j] = 0;
      }
    }
    stage_pr = mxGetPr(prhs[1]);
    for(i = 0; i < mxGetN(prhs[1]) - 1; i++){
      order[(int) stage_pr[i] - 1][(int) stage_pr[i + 1] - 1] = 1;
    }
  }

  /* Find the elimination ordering. */
  e = find_elim(n, adj_mat, order, -1);

  /* Allocate memory for the answer, and set the answer. */
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  answer_G_pr = mxGetPr(plhs[0]);
  cliques = get_cliques(e);
/* 
  dims[0] = 1;
  dims[1] = get_size_Map(cliques);
  plhs[1] = mxCreateCellArray(2, (const int *) dims);*/
  plhs[1] = mxCreateCellMatrix(get_size_Map(cliques), 1);
  fill_ins = get_fill_ins(e);
  fill_ins_mat = mxCreateDoubleMatrix(m, n, mxREAL);
  fill_ins_pr = mxGetPr(fill_ins_mat);

  for(i = 0; i < n; i++){
    for(j = 0; j < m; j++){
      index = j * m + i;
      answer_G_pr[index] = G_pr[index];
      if(fill_ins[i][j] > 0){
        answer_G_pr[index] = 1;
        fill_ins_pr[index] = 1;
      }
    }
  }
  mxDestroyArray(full);
// (original)  plhs[2] = mlfSparse(fill_ins_mat, NULL, NULL, NULL, NULL, NULL);
  plhs[2] = (mxArray *) mlfSparse(fill_ins_mat, NULL, NULL, NULL, NULL, NULL); // added typecasting
  mxDestroyArray(fill_ins_mat);
  nodes = (int *) malloc(sizeof(int) * n);
  k = 0;
  iter = get_Iterator(cliques);
  while(!is_empty(iter)){
    w = next_key(iter);
    clique = (Map) w.v;
    matlab_clique = mxCreateDoubleMatrix(1, get_size_Map(clique), mxREAL);
    matlab_clique_pr = mxGetPr(matlab_clique);
    for(i = 0; i < n; i++){
      nodes[i] = 0;
    }
    iter2 = get_Iterator(clique);
    while(!is_empty(iter2)){
      w2 = next_key(iter2);
      nodes[w2.i] = w2.i + 1;
    }
    j = 0;
    for(i = 0; i < n; i++){
      if(nodes[i] > 0){
        matlab_clique_pr[j++] = nodes[i];
      }
    }
    mxSetCell(plhs[1], k++, matlab_clique);
  }
  free(nodes);

  /* Finally, free the allocated memory. */
  destroy_Elimination(e);
  if(adj_mat){
    if(adj_mat[0]){
      free(adj_mat[0]);
    }
    free(adj_mat);
  }
  if(order){
    if(order[0]){
      free(order[0]);
    }
    free(order);
  }
}
