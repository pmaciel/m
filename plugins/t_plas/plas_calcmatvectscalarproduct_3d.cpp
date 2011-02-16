
#include "common.h"


/*
 * This function computes the 3D scalar product of a matrix
 * and a vector.
 */

void plas_CalcMatVectScalarProduct_3D(double *value, double **m, double *a)
{
  value[0] = m[0][0]*a[0]+m[0][1]*a[1]+m[0][2]*a[2];
  value[1] = m[1][0]*a[0]+m[1][1]*a[1]+m[1][2]*a[2];
  value[2] = m[2][0]*a[0]+m[2][1]*a[1]+m[2][2]*a[2];
}

