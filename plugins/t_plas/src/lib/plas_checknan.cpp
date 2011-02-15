
#include "common.h"


/*
 * This file includes the functionality to perform  trajectory
 * integrations of dispersed entities.
 *
 * This function performs a not-a-number check for the entity
 * position and velocity.
 */

void plas_CheckNaN(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent)
{
  int idim;

  double checkNanPos = 0.0;
  double checkNanVel = 0.0;

  for(idim=0; idim<data->fp.numDim; idim++){
    checkNanPos += ent->pos[idim];
    checkNanVel += ent->vel[idim];
  }

  if(isnan(checkNanPos) || isnan(checkNanVel)){
    ent->flag = DFLAG_DISABLED;
    data->sd.lost++;
    plasinterface_screenWarning((char*) "Not-a-number detected.\n");
  }
}

