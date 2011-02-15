
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes element neighbours.
 */

void plasdriver_CalcElementNeighbours(DRIVER_GAMBIT_MESH *dmesh)
{
  int neighbourFound,jnod,knod,lnod,ielm,jelm,kelm,lelm,ifac,faceNodes[4];

  dmesh->elmNeighbs = (int**)calloc(dmesh->numElm,sizeof(int*));
  for(ielm=0; ielm<dmesh->numElm; ielm++){
    dmesh->elmNeighbs[ielm] = (int*)calloc(dmesh->numElmFaces[ielm],sizeof(int));
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(dmesh,ielm,ifac,faceNodes);
      neighbourFound = 0;
      for(jnod=0; jnod<dmesh->numNodElms[faceNodes[0]]; jnod++){
        jelm = dmesh->nodElms[faceNodes[0]][jnod];
        for(knod=0; knod<dmesh->numNodElms[faceNodes[1]]; knod++){
          kelm = dmesh->nodElms[faceNodes[1]][knod];
          if(dmesh->numDim==2){
            if(jelm==kelm && jelm!=ielm){neighbourFound = 1;}
          }else if(dmesh->numDim==3){
            for(lnod=0; lnod<dmesh->numNodElms[faceNodes[2]]; lnod++){
              lelm = dmesh->nodElms[faceNodes[2]][lnod];
              if(jelm==kelm && jelm==lelm && jelm!=ielm){neighbourFound = 1;}
              if(neighbourFound==1){break;}
            }
          }
          if(neighbourFound==1){break;}
        }
        if(neighbourFound==1){break;}
      }
      if(neighbourFound==1){
        dmesh->elmNeighbs[ielm][ifac] = jelm;
      } else{
        dmesh->elmNeighbs[ielm][ifac] = -1;
      }
    }
  }
}

