
#include "common.h"


/*
 * This file includes the computation of back-coupling terms
 * from the dispersed phase to the continuous phase.
 *
 * This function computes the mass and momentum back coupling
 * terms for a dispersed entity.
 */

void plas_CalcBackCoupling(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double *force, double tFactor)
{
  double limitFrac = 0.95;
  int idim,jdim,inod;
  double contVol,iFrac,jFrac,rhsTerm,impFac[ent->edata.numElmNodes];

  //***Volume fraction trigger***//

  if(data->ip.volfracCoupl){
    if(data->pd[ent->node].volFrac>limitFrac){
      iFrac = limitFrac;
    } else{
      iFrac = data->pd[ent->node].volFrac;
    }
  } else{
    iFrac = 0.0;
  }

  if(data->pd[ent->node].volFrac>1.0){
    jFrac = 1.0/data->pd[ent->node].volFrac;
  } else{
    jFrac = 1.0;
  }

  //***Compute back-coupling force contribution (m/s^2) for fluid momentum equation (negative sign)***//

  if(data->ip.momentumCoupl==FORCE_PIC){

    inod = ent->node;
    contVol = plasinterface_getNodVol(ent->node);
    for(idim=0; idim<data->fp.numDim; idim++){
      data->pd[ent->node].dispForce[idim+1] -= tFactor*force[idim]*jFrac/((1.0-iFrac)*contVol*data->fp.rhoCont);
    }

  } else if(data->ip.momentumCoupl==FORCE_PROJ){

    plas_CalcNodeImpactFactors(data,ent,impFac);
    for(idim=0; idim<ent->edata.numElmNodes; idim++){
      inod = ent->edata.elmNodes[idim];
      if(data->fp.flowSolver==FLOWSOLVER_SFELES){inod %= data->fp.numNod;}
      contVol = plasinterface_getNodVol(inod);
      for(jdim=0; jdim<data->fp.numDim; jdim++){
        data->pd[inod].dispForce[jdim+1] -= tFactor*impFac[idim]*force[jdim]*jFrac/((1.0-iFrac)*contVol*data->fp.rhoCont);
      }
    }
  }

  //***Compute back-coupling term for fluid continuity equation***//

  if(data->ip.volfracCoupl){

    inod = ent->node;
    rhsTerm = data->pd[inod].volFracDt;
    for(idim=0; idim<data->fp.numDim; idim++){
      rhsTerm += flow->vel[idim]*data->pd[inod].volFracDx[idim];
    }

    data->pd[inod].dispForce[0] += tFactor*jFrac*rhsTerm/(1.0-iFrac);
    for(idim=1; idim<data->fp.numDim; idim++){
      data->pd[inod].dispForce[idim] += tFactor*jFrac*rhsTerm*flow->vel[idim]/(1.0-iFrac);
    }
  }
}

