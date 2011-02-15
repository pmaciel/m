
#include "common.h"


/*
 * This file includes the computation of back-coupling terms
 * from the dispersed phase to the continuous phase.
 *
 * This function calculates the momentum back coupling forces
 * for a particle or droplet.
 */

void plas_CalcCouplingForcesParticle(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor)
{
  int idim;
  double densityRatio,particleVol,particleMass;
  double uxw[data->fp.numDim],entForce[data->fp.numDim];

  //***Compute momentum coupling forces***//

  if(data->ip.momentumCoupl){

    particleVol = PI*ent->diam*ent->diam*ent->diam/6.0;
    particleMass = data->md.rhoDisp*particleVol;
    densityRatio = (data->fp.rhoCont/data->md.rhoDisp);

    if(data->fp.numDim==2){
      uxw[0] = 0.0;
      uxw[1] = 0.0;
    } else if(data->fp.numDim==3){
      uxw[0] = ent->relVel[1]*flow->vort[2]-ent->relVel[2]*flow->vort[1];
      uxw[1] = ent->relVel[2]*flow->vort[0]-ent->relVel[0]*flow->vort[2];
      uxw[2] = ent->relVel[0]*flow->vort[1]-ent->relVel[1]*flow->vort[0];
    }

    //**Compute surface force on particle (positive sign)***//

    for(idim=0; idim<data->fp.numDim; idim++){
      entForce[idim] = particleMass*(3.0/4.0)*(ent->dragCoeff/ent->diam)*densityRatio*ent->normVel*ent->relVel[idim]
                     + particleMass*ent->liftCoeff*densityRatio*uxw[idim];
    }
  }

  //***Compute back coupling terms***//

  plas_CalcBackCoupling(data,ent,flow,entForce,tFactor);

}

