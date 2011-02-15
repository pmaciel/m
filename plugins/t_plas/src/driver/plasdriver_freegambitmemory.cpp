
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function frees Gambit memory.
 */

void plasdriver_FreeGambitMemory(DRIVER_GAMBIT_MESH *dmesh)
{
  int inod,ielm,ibnd,ifac;

  free(dmesh->elmTypes);
  free(dmesh->numElmNodes);
  free(dmesh->numBndFaces);

  for(inod=0; inod<dmesh->numNod; inod++){
    free(dmesh->coords[inod]);
    free(dmesh->nodElms[inod]);
  }
  free(dmesh->coords);
  free(dmesh->nodElms);

  for(ielm=0; ielm<dmesh->numElm; ielm++){
    free(dmesh->elmNodes[ielm]);
    free(dmesh->elmNeighbs[ielm]);
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      free(dmesh->elmNorms[ielm][ifac]);
    }
    free(dmesh->elmNorms[ielm]);
  }
  free(dmesh->elmNodes);
  free(dmesh->elmNeighbs);
  free(dmesh->elmNorms);

  for(ibnd=0; ibnd<dmesh->numBnd; ibnd++){
    free(dmesh->bndNames[ibnd]);
    free(dmesh->bndFaces[ibnd]);
    free(dmesh->bndDomElms[ibnd]);
  }
  free(dmesh->bndNames);
  free(dmesh->bndFaces);
  free(dmesh->bndDomElms);
  free(dmesh->numNodElms);
  free(dmesh->numElmFaces);
  free(dmesh->elmVolumes);
  free(dmesh->nodVolumes);
  free(dmesh->bndTypes);
}

