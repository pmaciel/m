
#include "common.h"


/*
 * This function computes the 3D cross product of two vectors.
 */

void plas_CalcCrossProduct_3D(double *value, double *a, double *b)
{
  value[0] = a[1]*b[2] - a[2]*b[1];
  value[1] = a[2]*b[0] - a[0]*b[2];
  value[2] = a[0]*b[1] - a[1]*b[0];
}

