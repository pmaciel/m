
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes element normals.
 */

void plasdriver_CalcElementNormals(DRIVER_GAMBIT_MESH *dmesh)
{
  int ielm,ifac,fnodes[4];

  for(ielm=0; ielm<dmesh->numElm; ielm++){
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(dmesh,ielm,ifac,fnodes);

      if(dmesh->elmTypes[ielm]==ELM_SIMPLEX && dmesh->numDim==2){

        dmesh->elmNorms[ielm][ifac][0] = dmesh->coords[fnodes[0]][1] - dmesh->coords[fnodes[1]][1];
        dmesh->elmNorms[ielm][ifac][1] = dmesh->coords[fnodes[1]][0] - dmesh->coords[fnodes[0]][0];

      } else if((dmesh->elmTypes[ielm]==ELM_SIMPLEX && dmesh->numDim==3)
                || (dmesh->elmTypes[ielm]==ELM_PRISM && ifac>2)
                || (dmesh->elmTypes[ielm]==ELM_PYRAMID && ifac>0)){

        dmesh->elmNorms[ielm][ifac][0] =
          0.5*((dmesh->coords[fnodes[2]][1]-dmesh->coords[fnodes[0]][1])
              *(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
              -(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
              *(dmesh->coords[fnodes[2]][2]-dmesh->coords[fnodes[0]][2]));
        dmesh->elmNorms[ielm][ifac][1] =
          0.5*((dmesh->coords[fnodes[2]][2]-dmesh->coords[fnodes[0]][2])
              *(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
              -(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
              *(dmesh->coords[fnodes[2]][0]-dmesh->coords[fnodes[0]][0]));
        dmesh->elmNorms[ielm][ifac][2] =
          0.5*((dmesh->coords[fnodes[2]][0]-dmesh->coords[fnodes[0]][0])
              *(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
              -(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
              *(dmesh->coords[fnodes[2]][1]-dmesh->coords[fnodes[0]][1]));

      } else if(dmesh->elmTypes[ielm]==ELM_QUAD){

        dmesh->elmNorms[ielm][ifac][0] = 2.0*(dmesh->coords[fnodes[0]][1] - dmesh->coords[fnodes[1]][1]);
        dmesh->elmNorms[ielm][ifac][1] = 2.0*(dmesh->coords[fnodes[1]][0] - dmesh->coords[fnodes[0]][0]);

      } else if(dmesh->elmTypes[ielm]==ELM_HEX
                || (dmesh->elmTypes[ielm]==ELM_PRISM && ifac<=2)
                || (dmesh->elmTypes[ielm]==ELM_PYRAMID && ifac==0)){

        dmesh->elmNorms[ielm][ifac][0] =
            ((dmesh->coords[fnodes[3]][1]-dmesh->coords[fnodes[0]][1])
            *(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
            -(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
            *(dmesh->coords[fnodes[3]][2]-dmesh->coords[fnodes[0]][2]));
        dmesh->elmNorms[ielm][ifac][1] =
            ((dmesh->coords[fnodes[3]][2]-dmesh->coords[fnodes[0]][2])
            *(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
            -(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
            *(dmesh->coords[fnodes[3]][0]-dmesh->coords[fnodes[0]][0]));
        dmesh->elmNorms[ielm][ifac][2] =
            ((dmesh->coords[fnodes[3]][0]-dmesh->coords[fnodes[0]][0])
            *(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
            -(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
            *(dmesh->coords[fnodes[3]][1]-dmesh->coords[fnodes[0]][1]));
      }
    }
  }
}

