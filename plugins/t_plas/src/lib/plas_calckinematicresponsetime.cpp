
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the kinematic response time of the
 * dispersed entity.
 */

double plas_CalcKinematicResponseTime(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent)
{
  double tau = 0.;

  if(data->rp.flowType==FLOW_PARTIC || data->rp.flowType==FLOW_DROPLET){
    tau = 4.0*data->md.rhoDisp*ent->diam*ent->diam/(3.0*data->fp.muCont*ent->reynolds*ent->dragCoeff);
  } else if(data->rp.flowType==FLOW_BUBBLY){
    tau = 2.0*ent->diam*ent->diam/(3.0*data->fp.nuCont*ent->reynolds*ent->dragCoeff);
  }

  return tau;
}

