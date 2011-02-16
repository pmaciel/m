
#include "common.h"


/*
 * This function computes the length of a vector.
 */

double plas_CalcVectorLength(int numDim, double *a)
{
  return sqrt(plas_CalcVectScalarProduct(numDim,a,a));
}

