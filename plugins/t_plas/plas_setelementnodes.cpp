
#include "common.h"


/*
 * This function sets the nodes of the element assigned to a
 * dispersed entity.
 */

void plas_SetElementNodes(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  int type = plasinterface_getElementType(ent->elm);
  int inod;

  //***Set number of nodes accrding to element type***//

  if(type==ELM_SIMPLEX){
    ent->edata.numElmNodes = numDim+1;
  } else if(type==ELM_PRISM){
    ent->edata.numElmNodes = 6;
  } else if(type==ELM_QUAD){
    ent->edata.numElmNodes = 4;
  } else if(type==ELM_HEX){
    ent->edata.numElmNodes = 8;
  } else if(type==ELM_PYRAMID){
    ent->edata.numElmNodes = 5;
  }

  //***Get node elements from flow solver***//

  for(inod=0; inod<ent->edata.numElmNodes; inod++){
    ent->edata.elmNodes[inod] = plasinterface_getElmNode(ent->elm,inod);
  }
}

