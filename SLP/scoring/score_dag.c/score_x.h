//
// File: score_x.h
//
//
// by Darima <darrimma@yahoo.com>, 27/12/2005 
//

#ifndef score_x_h
#define mex_h

void score_family_x(double *data,double *sz,
    int ndata,int nsz,double* score);

void score_dag_x(double *data,double *sz,double *dag,
    int nsz,int ndata,double* score);

double gammaln(double xx);

#endif // score_x_h
