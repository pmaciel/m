
#include "common.h"

/*
 * This file includes allocation and deallocation routines
 * for PLaS data structures LOCAL_ENTITY_VARIABLES and
 * LOCAL_FLOW_VARIABLES.
 *
 * This function allocates data for an instance of the
 * structure LOCAL_ENTITY_VARIABLES, which stores the position,
 * velocity and grid element information of an entity.
 */

void plas_AllocateLocalEntityVar(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  ent->edata.elmNodes = new int[MAXELMNODES];
  ent->edata.elmNorms = new double*[MAXELMNORMS];
  for (int idim=0; idim<MAXELMNORMS; ++idim)
    ent->edata.elmNorms[idim] = new double[numDim];
  ent->edata.elmFaceVectors = new double*[MAXELMFACES];
  for (int idim=0; idim<MAXELMFACES; ++idim)
    ent->edata.elmFaceVectors[idim] = new double[numDim];
  ent->vel    = new double[numDim];
  ent->velOld = new double[numDim];
  ent->relVel = new double[numDim];
  ent->pos    = new double[numDim];
  ent->posOld = new double[numDim];

  ent->elm  = -1;
  ent->node = -1;
}

