
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes the nodal dual cell volumes.
 */

void plasdriver_CalcNodalVolumes(DRIVER_GAMBIT_MESH *dmesh)
{
  int inod,ielm;

  dmesh->nodVolumes = (double*)calloc(dmesh->numNod,sizeof(double));

  for(inod=0; inod<dmesh->numNod; inod++){
    for(ielm=0; ielm<dmesh->numNodElms[inod]; ielm++){
      dmesh->nodVolumes[inod] += dmesh->elmVolumes[dmesh->nodElms[inod][ielm]]/dmesh->numNodElms[inod];
    }
  }
}

