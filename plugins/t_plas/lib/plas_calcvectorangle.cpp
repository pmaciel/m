
#include "common.h"


/*
 * This function computes the angle between two vectors.
 */

double plas_CalcVectorAngle(int numDim, double *a, double *b)
{
  return acos(plas_CalcVectScalarProduct(numDim,a,b)/(plas_CalcVectorLength(numDim,a)*plas_CalcVectorLength(numDim,b)));
}

