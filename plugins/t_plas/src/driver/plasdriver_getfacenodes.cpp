
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function gets the nodes of a boundary face.
 */

void plasdriver_GetFaceNodes(DRIVER_GAMBIT_MESH *dmesh, int elm, int face, int *nodes)
{
  if(dmesh->elmTypes[elm]==ELM_SIMPLEX){
    if(dmesh->numDim==2){
      if(face==0){
        nodes[0] = dmesh->elmNodes[elm][0];
        nodes[1] = dmesh->elmNodes[elm][1];
        nodes[2] = -1;
        nodes[3] = -1;
      } else if(face==1){
        nodes[0] = dmesh->elmNodes[elm][1];
        nodes[1] = dmesh->elmNodes[elm][2];
        nodes[2] = -1;
        nodes[3] = -1;
      } else if(face==2){
        nodes[0] = dmesh->elmNodes[elm][2];
        nodes[1] = dmesh->elmNodes[elm][0];
        nodes[2] = -1;
        nodes[3] = -1;
      }
    } else if(dmesh->numDim==3){
      if(face==0){
        nodes[0] = dmesh->elmNodes[elm][1];
        nodes[1] = dmesh->elmNodes[elm][0];
        nodes[2] = dmesh->elmNodes[elm][2];
        nodes[3] = -1;
      } else if(face==1){
        nodes[0] = dmesh->elmNodes[elm][0];
        nodes[1] = dmesh->elmNodes[elm][1];
        nodes[2] = dmesh->elmNodes[elm][3];
        nodes[3] = -1;
      } else if(face==2){
        nodes[0] = dmesh->elmNodes[elm][1];
        nodes[1] = dmesh->elmNodes[elm][2];
        nodes[2] = dmesh->elmNodes[elm][3];
        nodes[3] = -1;
      } else if(face==3){
        nodes[0] = dmesh->elmNodes[elm][2];
        nodes[1] = dmesh->elmNodes[elm][0];
        nodes[2] = dmesh->elmNodes[elm][3];
        nodes[3] = -1;
      }
    }
  } else if(dmesh->elmTypes[elm]==ELM_QUAD){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][3];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = -1;
      nodes[3] = -1;
    }
  } else if(dmesh->elmTypes[elm]==ELM_HEX){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = dmesh->elmNodes[elm][5];
      nodes[3] = dmesh->elmNodes[elm][4];
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][3];
      nodes[2] = dmesh->elmNodes[elm][7];
      nodes[3] = dmesh->elmNodes[elm][5];
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][6];
      nodes[3] = dmesh->elmNodes[elm][7];
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = dmesh->elmNodes[elm][6];
    } else if(face==4){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][2];
      nodes[3] = dmesh->elmNodes[elm][3];
    } else if(face==5){
      nodes[0] = dmesh->elmNodes[elm][4];
      nodes[1] = dmesh->elmNodes[elm][5];
      nodes[2] = dmesh->elmNodes[elm][7];
      nodes[3] = dmesh->elmNodes[elm][6];
    }
  } else if(dmesh->elmTypes[elm]==ELM_PRISM){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = dmesh->elmNodes[elm][3];
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][5];
      nodes[3] = dmesh->elmNodes[elm][4];
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][3];
      nodes[3] = dmesh->elmNodes[elm][5];
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][1];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][4];
      nodes[2] = dmesh->elmNodes[elm][5];
      nodes[3] = -1;
    }
  } else if(dmesh->elmTypes[elm]==ELM_PYRAMID){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][3];
      nodes[3] = dmesh->elmNodes[elm][1];
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][3];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    }
  }
}

