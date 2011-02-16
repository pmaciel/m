
#include "common.h"


/*
 * This function normalizes a vector to length one.
 */

void plas_NormalizeVector(int numDim, double *a)
{
  int idim;
  double length = plas_CalcVectorLength(numDim,a);

  for(idim=0; idim<numDim; idim++){
    a[idim] /= length;
  }
}

