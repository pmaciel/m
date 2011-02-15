
#include "common.h"


/*
 * This function performs velocity interpolation from the
 * fluid flow solver.
 */

void plas_InterpolateVelocity(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim,jdim,kdim;
  double impFac[ent->edata.numElmNodes];

  //***Compute impact factors of element nodes***//

  plas_CalcNodeImpactFactors(data,ent,impFac);

  //***Flow velocity***//

  for(idim=0; idim<data->fp.numDim; idim++){
    flow->vel[idim] = 0.0;
    for(jdim=0; jdim<ent->edata.numElmNodes; jdim++){
      flow->vel[idim] += impFac[jdim]*((1.0-step)*plasinterface_getVelocityCompOld(ent->edata.elmNodes[jdim],idim)
                                           + step*plasinterface_getVelocityComp(ent->edata.elmNodes[jdim],idim));
    }
  }

  //***Flow velocity space derivatives***//

  for(idim=0; idim<data->fp.numDim; idim++){
    for(jdim=0; jdim<data->fp.numDim; jdim++){
      flow->velDx[idim][jdim] = 0.0;
      for(kdim=0; kdim<ent->edata.numElmNodes; kdim++){
        flow->velDx[idim][jdim] += impFac[kdim]*((1.0-step)*plasinterface_getVelocityDerivativeCompOld(ent->edata.elmNodes[kdim],idim,jdim)
                                                     + step*plasinterface_getVelocityDerivativeComp(ent->edata.elmNodes[kdim],idim,jdim));
      }
    }
  }

  //***Flow velocity time derivative***//

  for(idim=0; idim<data->fp.numDim; idim++){
    flow->velDt[idim] = 0.0;
    if (data->fp.dtEul>1.e-20) {
      for(jdim=0; jdim<ent->edata.numElmNodes; jdim++){
        flow->velDt[idim] += impFac[jdim]*(plasinterface_getVelocityComp(ent->edata.elmNodes[jdim],idim)
                                         - plasinterface_getVelocityCompOld(ent->edata.elmNodes[jdim],idim))/data->fp.dtEul;
      }
    }
  }
}

