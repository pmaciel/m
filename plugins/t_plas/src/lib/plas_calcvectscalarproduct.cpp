
#include "common.h"


/*
 * This function computes the scalar product of two vectors.
 */

double plas_CalcVectScalarProduct(int numDim, double *a, double *b)
{
  if (numDim==2)
    return a[0]*b[0]+a[1]*b[1];
  else if (numDim==3)
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return 0.;
}

