
#include "common.h"


/*
 * This file includes the functionality to perform  trajectory
 * integrations of dispersed entities.
 *
 * This function composes the trajectory equation.
 */

void plas_CalcTrajectory(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double dtLagr)
{
  int idim,jdim;
  double
    vmCoeff,
    pvmTerm,
    convTerm,
    radTerm,
    massTerm,
    concTerm,
    solvec[data->fp.numDim+2],
    rhsvec[data->fp.numDim+2],
    eps       = 1e-9,
    theta     = 0.5,
    boltzmann = 5.670400e-8;


  double **mat = new double*[data->fp.numDim+2];
  for (int r=0; r<data->fp.numDim+2; ++r)
    mat[r] = new double[data->fp.numDim+2];


  // TODO: Concentration definition depends on saturation model
  double alpha = 1.25;
  double flow_concentration = flow->pressure*data->md.HeDisp*alpha;

  //***Initialize data structures (u,v,w,T,d)***//

  for(idim=0; idim<data->fp.numDim; idim++){
    solvec[idim] = ent->vel[idim];
  }
  solvec[data->fp.numDim] = ent->temp;
  solvec[data->fp.numDim+1] = ent->diam;

  for(idim=0; idim<data->fp.numDim+2; idim++){
    for(jdim=0; jdim<data->fp.numDim+2; jdim++){
      mat[idim][jdim] = 0.0;
    }
    rhsvec[idim] = 0.0;
  }

  //***Add time contribution to matrix***//

  for(idim=0; idim<data->fp.numDim+2; idim++){
    mat[idim][idim] += 1.0/dtLagr;
  }

  //***Drag force contribution to matrix and RHS***//

  for(idim=0; idim<data->fp.numDim; idim++){
    mat[idim][idim] += theta/ent->kinRespTime;
    rhsvec[idim] += ent->relVel[idim]/ent->kinRespTime;
  }

  //***Lift force contribution to matrix and RHS***//

  if(data->ip.liftForce && data->rp.flowType==FLOW_BUBBLY && data->fp.numDim==3){
    mat[0][1] += theta*2.0*ent->liftCoeff*flow->vort[2];
    mat[0][2] -= theta*2.0*ent->liftCoeff*flow->vort[1];
    mat[1][0] -= theta*2.0*ent->liftCoeff*flow->vort[2];
    mat[1][2] += theta*2.0*ent->liftCoeff*flow->vort[0];
    mat[2][0] += theta*2.0*ent->liftCoeff*flow->vort[1];
    mat[2][1] -= theta*2.0*ent->liftCoeff*flow->vort[0];
    rhsvec[0] += 2.0*ent->liftCoeff*(ent->relVel[1]*flow->vort[2]-ent->relVel[2]*flow->vort[1]);
    rhsvec[1] += 2.0*ent->liftCoeff*(ent->relVel[2]*flow->vort[0]-ent->relVel[0]*flow->vort[2]);
    rhsvec[2] += 2.0*ent->liftCoeff*(ent->relVel[0]*flow->vort[1]-ent->relVel[1]*flow->vort[0]);
  }

  //***Virtual mass & pressure force contribution to RHS***//

  if(data->rp.flowType==FLOW_BUBBLY){
    vmCoeff = 0.5*(1.0+2.786*data->pd[ent->node].volFrac);
    pvmTerm = 2.0*(1.0+0.5*(1.0+2.786*data->pd[ent->node].volFrac));
    for(idim=0; idim<data->fp.numDim; idim++){
      rhsvec[idim] += pvmTerm*flow->velDt[idim];
      for(jdim=0; jdim<data->fp.numDim; jdim++){
        rhsvec[idim] += pvmTerm*flow->vel[jdim]*flow->velDx[idim][jdim];
      }
    }
  }

  //***Gravitational force contribution to RHS***//

  for(idim=0; idim<data->fp.numDim; idim++){
    if(data->rp.flowType==FLOW_PARTIC || data->rp.flowType==FLOW_DROPLET){
      rhsvec[idim] += data->ip.gravVec[idim];
    } else if(data->rp.flowType==FLOW_BUBBLY){
      rhsvec[idim] -= 2.0*data->ip.gravVec[idim];
    }
  }

  //***Temperature equation contribution to matrix and RHS***//

  if(data->ip.energyCoupl && (data->rp.flowType==FLOW_PARTIC || data->rp.flowType==FLOW_DROPLET)){
    convTerm = ent->nusselt/(2.0*ent->thermRespTime);
    radTerm = (6.0*boltzmann*data->md.epsDisp)/(data->md.rhoDisp*data->md.cpDisp*ent->diam);
    if(data->ip.evapModel && data->rp.flowType==FLOW_DROPLET){
      massTerm = (3.0)*(ent->massTrCoeff/pow(ent->diam,2.0))*(data->md.latHeatDisp/data->md.cpDisp);
    } else{
      massTerm = 0.0;
    }
    mat[data->fp.numDim][data->fp.numDim] += theta*(convTerm-radTerm*4.0*pow(ent->temp,3.0));
    mat[data->fp.numDim][data->fp.numDim+1] -=
      theta*(-((1.0/(4.0*ent->thermRespTime*ent->diam))*(2.0+3.0*ent->nusselt)*(ent->relTemp))
   + ((radTerm/ent->diam)*(pow(ent->temp,4.0)-pow(flow->temp,4.0)))-((massTerm/ent->diam)*(1.5+1.0/ent->sherwood)));
    rhsvec[data->fp.numDim] += (convTerm*ent->relTemp)-(radTerm*(pow(ent->temp,4.0)-pow(flow->temp,4.0)))-massTerm;
  }

  //***Diameter equation contribution to matrix and RHS***//

  if(data->ip.evapModel && data->rp.flowType==FLOW_DROPLET ){
    mat[data->fp.numDim+1][data->fp.numDim+1] -= theta*(ent->massTrCoeff/pow(ent->diam,2.0))*(1.5-1.0/ent->sherwood);
    rhsvec[data->fp.numDim+1] += -ent->massTrCoeff/ent->diam;
  }

  //***Concentration Boundary Layer Model Payvar for bubble growth and Epstein & Plesset for Bubble shrink

  if(data->ip.saturModel && data->rp.flowType==FLOW_BUBBLY ){
   concTerm  = 4.0*data->md.massDiffCoeff*data->fp.rhoCont*(flow_concentration-ent->concInterf)/(ent->rhoBubble*(data->fp.rhoCont-ent->concInterf));

    if(flow_concentration > ent->concInterf){
      rhsvec[data->fp.numDim+1] += concTerm/ent->diam*(1.0+1.0/(pow(1.0+(flow_concentration-data->fp.rhoCont)/(ent->concInterf-flow_concentration)*(ent->rhoBubble/data->fp.rhoCont)*(1.0-(data->md.rhoDisp/ent->rhoBubble)*pow(2.0e-5/ent->diam,3.0)),0.5)-1.0));
      mat[data->fp.numDim+1][data->fp.numDim+1] -= concTerm*((1.0/(ent->diam+eps)+1.0/((ent->diam+eps)*pow(1.0+(flow_concentration-data->fp.rhoCont)/(ent->concInterf-flow_concentration)*(ent->rhoBubble/data->fp.rhoCont)*(1.0-(data->md.rhoDisp/ent->rhoBubble)*pow(2.0e-5/(ent->diam+eps),3.0)),0.5)-1.0))-(1.0/ent->diam+1.0/(ent->diam*pow(1.0+(flow_concentration-data->fp.rhoCont)/(ent->concInterf-flow_concentration)*(ent->rhoBubble/data->fp.rhoCont)*(1.0-(data->md.rhoDisp/ent->rhoBubble)*pow(2.0e-5/ent->diam,3.0)),0.5)-1.0)))/eps;

      }else{
      rhsvec[data->fp.numDim+1] += concTerm/ent->diam;
      mat[data->fp.numDim+1][data->fp.numDim+1] -= ((concTerm/(ent->diam+eps))-(concTerm/ent->diam))/eps;
        }
  }

  //***Solve linear system for velocity, temperature & diameter***//

  plas_SolveGaussSeidel(data->fp.numDim+2,mat,solvec,rhsvec);

  for(idim=0; idim<data->fp.numDim; idim++){
    ent->velOld[idim] = ent->vel[idim];
    ent->vel[idim] += solvec[idim];
  }
  ent->temp += solvec[data->fp.numDim];
  ent->diam += solvec[data->fp.numDim+1];

  //***Update entity position***//

  for(idim=0; idim<data->fp.numDim; idim++){
    ent->posOld[idim] = ent->pos[idim];
    ent->pos[idim] += 0.5*(ent->velOld[idim]+ent->vel[idim])*dtLagr;
  }


  for (int r=0; r<data->fp.numDim+2; ++r)
    delete[] mat[r];
  delete[] mat;
}

