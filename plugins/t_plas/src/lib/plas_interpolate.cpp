
#include "common.h"


/*
 * This function manages the interpolation of variables from
 * the fluid flow solver.
 */

void plas_Interpolate(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim;

  //***Interpolate velocity, temperature and pressure***//

  plas_InterpolateVelocity(data,ent,flow,step);
  plas_CalcVorticity(data->fp.numDim,flow);
  plas_InterpolateTemperature(data,ent,flow,step);
  plas_InterpolatePressure(data,ent,flow,step);

  //***Compute relative and normal velocity***//

  ent->normVel = 0.0;
  for(idim=0; idim<data->fp.numDim; idim++){
    ent->relVel[idim] = flow->vel[idim]-ent->vel[idim];
    ent->normVel += ent->relVel[idim]*ent->relVel[idim];
  }
  ent->normVel = sqrt(ent->normVel);

  //***Compute relative temperature***//

  ent->relTemp = flow->temp-ent->temp;
}

