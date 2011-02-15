
#include "common.h"

/*
 * This file includes allocation and deallocation routines
 * for PLaS data structures LOCAL_ENTITY_VARIABLES and
 * LOCAL_FLOW_VARIABLES.
 *
 * This function de-allocates data for an instance of the type
 * LOCAL_ENTITY_VARIABLES.
 */

void plas_DeallocateLocalEntityVar(LOCAL_ENTITY_VARIABLES *ent)
{
  int idim;

  //***De-allocate data***//

  free(ent->edata.elmNodes);
  for(idim=0; idim<MAXELMNORMS; idim++){
    free(ent->edata.elmNorms[idim]);
  }
  free(ent->edata.elmNorms);
  for(idim=0; idim<MAXELMFACES; idim++){
    free(ent->edata.elmFaceVectors[idim]);
  }
  free(ent->edata.elmFaceVectors);
  free(ent->vel);
  free(ent->velOld);
  free(ent->relVel);
  free(ent->pos);
  free(ent->posOld);
}

