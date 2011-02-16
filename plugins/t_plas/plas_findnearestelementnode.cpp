
#include "common.h"


/*
 * This file includes all functinality to perform an element
 * search for a dispersed entity.
 *
 * This function finds the closest element node to an entity.
 */

int plas_FindNearestElementNode(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent)
{
  int idim,node=0;
  double impFac[ent->edata.numElmNodes],f_max=0.;

  plas_CalcNodeImpactFactors(data,ent,impFac);
  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    if(idim==0){
      node = ent->edata.elmNodes[idim];
      f_max = impFac[idim];
    } else{
      if(impFac[idim]>f_max){
        node = ent->edata.elmNodes[idim];
        f_max = impFac[idim];
      }
    }
  }

  if(data->fp.flowSolver==FLOWSOLVER_SFELES){node %= data->fp.numNod;}

  return node;
}

