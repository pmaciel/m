
#include "common.h"


/*
 * This file includes allocation and deallocation routines
 * for PLaS data structures LOCAL_ENTITY_VARIABLES and
 * LOCAL_FLOW_VARIABLES.
 *
 * This function de-allocates data for an instance of the
 * structure LOCAL_FLOW_VARIABLES.
 */

void plas_DeallocateLocalFlowVar(int numDim, LOCAL_FLOW_VARIABLES *flow)
{
  delete[] flow->vel;
  delete[] flow->velDt;
  for (int idim=0; idim<numDim; ++idim)
    delete[] flow->velDx[idim];
  delete[] flow->velDx;
  delete[] flow->vort;
}
