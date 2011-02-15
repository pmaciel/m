
#include "common.h"


/*
 * This file includes all functinality to perform an element
 * search for a dispersed entity.
 *
 * This function computes the node impact factors of an
 * element according to the entity position.
 */

void plas_CalcNodeImpactFactors(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, double *imp)
{
  double sum = 0.0;
  int idim,jdim,inod;
  double distance[data->fp.numDim];

  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    inod = ent->edata.elmNodes[idim];
    if(data->fp.flowSolver==FLOWSOLVER_SFELES){inod %= data->fp.numNod;}
    for(jdim=0; jdim<data->fp.numDim; jdim++){
      distance[jdim] = plasinterface_getNodCoord(inod,jdim)-ent->pos[jdim];
    }
    imp[idim] = plas_CalcVectorLength(data->fp.numDim,distance);
    if(imp[idim]<data->rp.errTol){imp[idim] = data->rp.errTol;}
    sum += 1.0/imp[idim];
  }

  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    imp[idim] = 1.0/(imp[idim]*sum);
  }
}

