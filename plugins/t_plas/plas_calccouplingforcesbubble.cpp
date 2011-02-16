
#include "common.h"


/*
 * This file includes the computation of back-coupling terms
 * from the dispersed phase to the continuous phase.
 *
 * This function computes the momentum back coupling forces
 * for a bubble.
 */

void plas_CalcCouplingForcesBubble(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor)
{
  int idim,jdim;
  double densityRatio,bubbleVol,bubbleMass,vmCoeff;
  double uxw[data->fp.numDim],dudt[data->fp.numDim],dvdt[data->fp.numDim],entForce[data->fp.numDim];

  //***Compute momentum coupling forces***//

  if(data->ip.momentumCoupl){

    bubbleVol = PI*ent->diam*ent->diam*ent->diam/6.0;
    bubbleMass = data->md.rhoDisp*bubbleVol;
    densityRatio = (data->fp.rhoCont/data->md.rhoDisp);

    if(data->fp.numDim==2){
      uxw[0] = 0.0;
      uxw[1] = 0.0;
    } else if(data->fp.numDim==3){
      uxw[0] = ent->relVel[1]*flow->vort[2]-ent->relVel[2]*flow->vort[1];
      uxw[1] = ent->relVel[2]*flow->vort[0]-ent->relVel[0]*flow->vort[2];
      uxw[2] = ent->relVel[0]*flow->vort[1]-ent->relVel[1]*flow->vort[0];
    }

    vmCoeff = 0.5*(1.0+2.786*data->pd[ent->node].volFrac);

    for(idim=0; idim<data->fp.numDim; idim++){
      dudt[idim] = flow->velDt[idim];
      for(jdim=0; jdim<data->fp.numDim; jdim++){
        dudt[idim] += flow->vel[jdim]*flow->velDx[idim][jdim];
      }
      dvdt[idim] = (ent->vel[idim]-ent->velOld[idim])/data->fp.dtEul;
    }

    //**Compute surface force on particle (positive sign)***//

    for(idim=0; idim<data->fp.numDim; idim++){
      entForce[idim] = bubbleMass*(3.0/4.0)*(ent->dragCoeff/ent->diam)*densityRatio*ent->normVel*ent->relVel[idim]
                     + bubbleMass*ent->liftCoeff*densityRatio*uxw[idim]
                     + bubbleMass*vmCoeff*densityRatio*(dudt[idim]-dvdt[idim]);
    }
  }

  //***Compute back coupling terms***//

  plas_CalcBackCoupling(data,ent,flow,entForce,tFactor);

}

