#ifndef __UTILS__
#define __UTILS__

#include <string>
#include <iostream>

/* forward declarations */
void copy(double *U1, double *U2, int N);
void nrerror(std::string msg);
void vecprd3(double *U1, double *U2, double *VP);
double ***d3tensor(long nrow, long ncol, long ndep);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_d3tensor(double ***t);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);

#endif

