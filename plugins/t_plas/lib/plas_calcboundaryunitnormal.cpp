
#include "common.h"


/*
 * This file includes all functionality concerning entities
 * interacting with boundaries of the computational domain,
 * such as wall bounces and entities leaving through outlets.
 *
 * This function calculates the unit normal vector of a
 * boundary face.
 */

void plas_CalcBoundaryUnitNormal(int numDim, int ibnd, int ifac, double *unitVec)
{
  int idim;
  double length;

  //***Get normal vector of a boundary face***//

  for(idim=0; idim<numDim; idim++){
    unitVec[idim] = plasinterface_getBndFaceNormComp(ibnd,ifac,idim);
  }

  //***Divide normal vector components by vector length***//

  length = plas_CalcVectorLength(numDim,unitVec);
  for(idim=0; idim<numDim; idim++){
    unitVec[idim] /= length;
  }
}

