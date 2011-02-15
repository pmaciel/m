
#include "common.h"


/*
 * This file includes allocation and deallocation routines
 * for PLaS data structures LOCAL_ENTITY_VARIABLES and
 * LOCAL_FLOW_VARIABLES.
 *
 * This function allocates data for an instance of the
 * structure LOCAL_FLOW_VARIABLES, which stores the flow
 * velocity and its derivatives at the entity position.
 */

void plas_AllocateLocalFlowVar(int numDim, LOCAL_FLOW_VARIABLES *flow)
{
  int idim;

  //***Allocate data***//

  flow->vel = (double*)calloc(numDim,sizeof(double));
  flow->velDt = (double*)calloc(numDim,sizeof(double));
  flow->velDx = (double**)calloc(numDim,sizeof(double*));
  for(idim=0; idim<numDim; idim++){
    flow->velDx[idim] = (double*)calloc(numDim,sizeof(double));
  }
  flow->vort = (double*)calloc(numDim,sizeof(double));
}

