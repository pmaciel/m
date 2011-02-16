
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
  delete[] ent->edata.elmNodes;
  for (int idim=0; idim<MAXELMNORMS; ++idim)
    delete[] ent->edata.elmNorms[idim];
  delete[] ent->edata.elmNorms;
  for (int idim=0; idim<MAXELMFACES; ++idim)
    delete[] ent->edata.elmFaceVectors[idim];
  delete[] ent->edata.elmFaceVectors;
  delete[] ent->vel;
  delete[] ent->velOld;
  delete[] ent->relVel;
  delete[] ent->pos;
  delete[] ent->posOld;
}
