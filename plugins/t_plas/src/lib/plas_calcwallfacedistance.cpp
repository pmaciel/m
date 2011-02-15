
#include "common.h"


/*
 * This file includes all functionality concerning entities
 * interacting with boundaries of the computational domain,
 * such as wall bounces and entities leaving through outlets.
 *
 * This function calculates the distance of an entity to a
 * wall face.
 */

double plas_CalcWallFaceDistance(int numDim, double *pos, int ibnd, int ifac)
{
  int idim;
  double posVec[numDim],unitVec[numDim];

  //***Calculate unit normal of boundary segment***//

  plas_CalcBoundaryUnitNormal(numDim,ibnd,ifac,unitVec);

  //***Calculate and return wall face distance***//

  for(idim=0; idim<numDim; idim++){
    posVec[idim] = pos[idim]-plasinterface_getBndFaceRefCoord(ibnd,ifac,idim);
  }

  return plas_CalcVectScalarProduct(numDim,posVec,unitVec);
}

