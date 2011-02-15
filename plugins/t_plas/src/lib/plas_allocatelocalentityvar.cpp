
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
  int idim;

  //***Allocate data***//

  ent->edata.elmNodes = (int*)calloc(MAXELMNODES,sizeof(int));
  ent->edata.elmNorms = (double**)calloc(MAXELMNORMS,sizeof(double*));
  for(idim=0; idim<MAXELMNORMS; idim++){
    ent->edata.elmNorms[idim] = (double*)calloc(numDim,sizeof(double));
  }
  ent->edata.elmFaceVectors = (double**)calloc(MAXELMFACES,sizeof(double*));
  for(idim=0; idim<MAXELMFACES; idim++){
    ent->edata.elmFaceVectors[idim] = (double*)calloc(numDim,sizeof(double));
  }
  ent->vel = (double*)calloc(numDim,sizeof(double));
  ent->velOld = (double*)calloc(numDim,sizeof(double));
  ent->relVel = (double*)calloc(numDim,sizeof(double));
  ent->pos = (double*)calloc(numDim,sizeof(double));
  ent->posOld = (double*)calloc(numDim,sizeof(double));

  ent->elm = -1;
  ent->node = -1;
}

