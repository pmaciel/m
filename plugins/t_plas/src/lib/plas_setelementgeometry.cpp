
#include "common.h"


/*
 * This function sets the geometry of an element assigned to
 * a dispersed entity.
 */

void plas_SetElementGeometry(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  plas_SetElementNodes(numDim,ent);
  plas_SetElementFaces(numDim,ent);
}

