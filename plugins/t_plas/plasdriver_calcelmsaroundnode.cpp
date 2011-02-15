
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes the elements around a node.
 */

#define MAXNODELMS 50
void plasdriver_CalcElmsAroundNode(DRIVER_GAMBIT_MESH *dmesh)
{
  int inod,ielm;

  dmesh->numNodElms = (int*)calloc(dmesh->numNod,sizeof(int));

  dmesh->nodElms = (int**)calloc(dmesh->numNod,sizeof(int*));
  for(inod=0; inod<dmesh->numNod; inod++){
    dmesh->nodElms[inod] = (int*)calloc(MAXNODELMS,sizeof(int));
  }

  for(ielm=0; ielm<dmesh->numElm; ielm++){
    for(inod=0; inod<dmesh->numElmNodes[ielm]; inod++){
      dmesh->nodElms[dmesh->elmNodes[ielm][inod]][dmesh->numNodElms[dmesh->elmNodes[ielm][inod]]] = ielm;
      dmesh->numNodElms[dmesh->elmNodes[ielm][inod]]++;
    }
  }
}
#undef MAXNODELMS

