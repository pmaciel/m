
#include "common.h"


/*
 * This function sets the faces and normal vectors of the
 * element assigned to a dispersed entity.
 */

void plas_SetElementFaces(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  int type = plasinterface_getElementType(ent->elm);
  int ifac,idim;

  //***Set number of faces accrding to element type***//

  if(type==ELM_SIMPLEX){
    ent->edata.numElmFaces = numDim+1;
  } else if(type==ELM_PRISM){
    ent->edata.numElmFaces = 5;
  } else if(type==ELM_QUAD){
    ent->edata.numElmFaces = 4;
  } else if(type==ELM_HEX){
    ent->edata.numElmFaces = 6;
  } else if(type==ELM_PYRAMID){
    ent->edata.numElmFaces = 5;
  }

  //***Get faces and normals from flow solver***//

  for(ifac=0; ifac<ent->edata.numElmFaces; ifac++){
    for(idim=0; idim<numDim; idim++){
      ent->edata.elmFaceVectors[ifac][idim] = plasinterface_getElmFaceMiddlePoint(ent->elm,ifac,idim);
      ent->edata.elmNorms[ifac][idim] = plasinterface_getElmNormComp(ent->elm,ifac,idim);
    }
  }
}

