
#include "common.h"


/*
 * This function computes the 2D scalar product of a matrix
 * and a vector.
 */

void plas_CalcMatVectScalarProduct_2D(double *value, double **m, double *a)
{
  value[0] = m[0][0]*a[0]+m[0][1]*a[1];
  value[1] = m[1][1]*a[0]+m[1][1]*a[1];
}

