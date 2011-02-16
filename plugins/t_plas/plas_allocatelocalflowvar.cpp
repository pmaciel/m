
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
  flow->vel   = new double[numDim];
  flow->velDt = new double[numDim];
  flow->velDx = new double*[numDim];
  for (int idim=0; idim<numDim; ++idim)
    flow->velDx[idim] = new double[numDim];
  flow->vort = new double[numDim];
}
