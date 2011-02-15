
#include "common.h"


/*
 * This file includes all functionality concerning entities
 * interacting with boundaries of the computational domain,
 * such as wall bounces and entities leaving through outlets.
 *
 * This function performs a wall bounce for an entity that
 * crossed a wall face of the domain boundary in the last
 * trajectory updata. Position and velocity are corrected.
 */

void plas_WallBounce(int numDim, double elasticity, LOCAL_ENTITY_VARIABLES *ent, int ibnd, int ifac)
{
  int idim;
  double unitVec[numDim];
  double wallDistancePositionVector,wallDistanceVelocityVector;

  //***Calculate unit normal of boundary segment***//

  plas_CalcBoundaryUnitNormal(numDim,ibnd,ifac,unitVec);

  //***Calculate normal wall distance of entity position***//

  wallDistancePositionVector = plas_CalcWallFaceDistance(numDim,ent->pos,ibnd,ifac);

  //***Adapt velocity due to elasticity factor***//

  for(idim=0; idim<numDim; idim++){
    ent->vel[idim] *= elasticity;
  }

  //***Add position vector to velocity vector***//

  for(idim=0; idim<numDim; idim++){
    ent->vel[idim] += ent->pos[idim];
  }

  //***Calculate normal wall distance of entity velocity vector***//

  wallDistanceVelocityVector = plas_CalcWallFaceDistance(numDim,ent->vel,ibnd,ifac);

  //***Mirror position and velocity vectors on wall segment***//

  for(idim=0; idim<numDim; idim++){
    ent->pos[idim] = ent->pos[idim]-2.0*wallDistancePositionVector*unitVec[idim];
    ent->vel[idim] = ent->vel[idim]-2.0*wallDistanceVelocityVector*unitVec[idim]-ent->pos[idim];
  }
}

