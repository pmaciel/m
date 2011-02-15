
#include "common.h"


/*
 * This function computes the vorticity of the fluid flow.
 */

void plas_CalcVorticity(int numDim, LOCAL_FLOW_VARIABLES *flow)
{
  if(numDim==2){
    flow->vort[0] = 0.0;
    flow->vort[1] = 0.0;
  } else if(numDim==3){
    flow->vort[0] = flow->velDx[2][1]-flow->velDx[1][2];
    flow->vort[1] = flow->velDx[0][2]-flow->velDx[2][0];
    flow->vort[2] = flow->velDx[1][0]-flow->velDx[0][1];
  }
}

