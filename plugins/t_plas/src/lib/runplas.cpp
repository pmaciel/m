
#include "common.h"


/*
 * This is the main routine of the PLaS solver. It has to be
 * called at each time step of the driving flow solver.
 */

void runPLaS(PLAS_DATA *data)
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  char errMessage[100];
  int idx,ient,jent,inod,idim,iunk,ibnd,ifac,subIter,facFound,avgctr;
  double cellSize,minCellSize,entVel,dtLagr,dtRemaining,unitVec[data->fp.numDim];

  //***Flow solver parameters that have to be set at every time step***//

  plasinterface_screenOutput((char*) "Initialization...\n");
  plasinterface_setFlowSolverParamOnTimeStep(&(data->fp));

  //***Allocate dynamic memory***//

  plas_AllocateLocalEntityVar(data->fp.numDim,&ent);
  plas_AllocateLocalFlowVar(data->fp.numDim,&flow);

  //***Initializations***//

  data->sd.enabled = 0;
  data->sd.in = 0;
  data->sd.out = 0;
  data->sd.bounce = 0;
  data->sd.periodic = 0;
  data->sd.passed = 0;
  data->sd.lost = 0;
  data->sd.coll = 0;
  data->sd.coalesc = 0;
  data->sd.leftproc = 0;
  data->sd.dtLagrAvg = 0.0;
  data->sd.reynoldsAvg = 0.0;
  data->sd.nusseltAvg = 0.0;
  data->sd.subIterAvg = 0.0;
  avgctr = 0;
  minCellSize = pow(data->fp.minElmVolume,(1.0/data->fp.numDim));

  if(data->fp.iter==1){
    plas_CalcCellwiseData(data);
  }
  plas_MpiBarrier();

  for(inod=0; inod<data->fp.numNod; inod++){
    for(iunk=0; iunk<data->fp.numUnk; iunk++){
      data->pd[inod].dispForce[iunk] = 0.0;
    }
  }

  //***Sort and count active entities***//

  idx = data->ip.numMaxEnt-1;
  for(ient=0; ient<data->ip.numMaxEnt; ient++){
    if(data->ed[ient].flag==DFLAG_ENABLED || data->ed[ient].flag==DFLAG_CREATED){
      data->sd.enabled++;
      continue;
    }
    for(jent=idx; jent>ient; jent--){
      if(data->ed[jent].flag==DFLAG_DISABLED){
        idx--;
        continue;
      }
      data->ed[ient].element = data->ed[jent].element;
      data->ed[ient].node = data->ed[jent].node;
      for(idim=0; idim<data->fp.numDim; idim++){
        data->ed[ient].position[idim] = data->ed[jent].position[idim];
        data->ed[ient].velocity[idim] = data->ed[jent].velocity[idim];
      }
      data->ed[ient].diameter = data->ed[jent].diameter;
      data->ed[ient].flag = data->ed[jent].flag;
      data->ed[jent].flag = DFLAG_DISABLED;
      idx = jent-1;
      data->sd.enabled++;
      break;
    }
  }
  plas_MpiBarrier();

  //***Entity production***//

  if(data->ip.numProdDom>0){
    plasinterface_screenOutput((char*) "Imposing new entities in production domains...\n");
    plas_ImposeProductionDomains(data);
  }
  plas_MpiBarrier();

  if(plas_MpiAllSumInt(data->numExtEnt)>0){
    plasinterface_screenOutput((char*) "Imposing externally generated entities...\n");
    plas_ImposeExternal(data);
  }
  plas_MpiBarrier();

  //***Loop over dispersed entities to update trajectories***//

  plasinterface_screenOutput((char*) "Updating trajectories...\n");
  for(ient=0; ient<data->sd.enabled; ient++){

    //***Get entity information from global data structure***//

    ent.flag = data->ed[ient].flag;
    if(ent.flag!=DFLAG_ENABLED && ent.flag!=DFLAG_CREATED){
      sprintf(errMessage,"Entity %i had bad flag %d.",ient,ent.flag);
      plas_TerminateOnError(errMessage);
    }

    ent.elm = data->ed[ient].element;
    ent.node = data->ed[ient].node;
    ent.diam = data->ed[ient].diameter;
    ent.temp = data->ed[ient].temperature;

    for(idim=0; idim<data->fp.numDim; idim++){
      ent.pos[idim] = data->ed[ient].position[idim];
      ent.vel[idim] = data->ed[ient].velocity[idim];
      data->ed[ient].velocityOld[idim] = data->ed[ient].velocity[idim];
    }

    //***Subiterations of the trajectory increment***//

    subIter = 0;
    dtRemaining = data->fp.dtEul;

    do{

      //***Determine Lagrangian time step size***//

      entVel = plas_CalcVectorLength(data->fp.numDim,ent.vel);
      if(entVel<data->rp.errTol){
        entVel = data->rp.errTol;
      }

      cellSize = pow(plasinterface_getElmVol(ent.elm),(1.0/data->fp.numDim));
      if(cellSize<data->rp.errTol){
        cellSize = data->rp.errTol;
      }

      dtLagr = data->rp.lagrTimeFactor*(cellSize/entVel);
      if(dtLagr>dtRemaining){
        dtLagr = dtRemaining;
      }

      dtRemaining -= dtLagr;
      subIter++;

      //***Entity production***//

      if(ent.flag==DFLAG_CREATED){
        if(plas_RandomDouble()<=((data->fp.dtEul-dtRemaining)/data->fp.dtEul)){
          ent.flag = DFLAG_ENABLED;
        }
      }

      if(ent.flag!=DFLAG_ENABLED){
        continue;
      }

      //***Interpolate continuous phase velocity and pressure at entity position***//

      plas_SetElementGeometry(data->fp.numDim,&ent);
      plas_Interpolate(data,&ent,&flow,(data->fp.dtEul-dtRemaining)/data->fp.dtEul);
      plas_CalcEntityCoefficients(data,&ent,&flow);

      //***Solve Lagrangian trajectory equation***//

      plas_CalcTrajectory(data,&ent,&flow,dtLagr);
      plas_CheckNaN(data,&ent);

      //***Check for evaporated, burned or collapsed entity***//

      if(ent.diam<data->rp.errTol){
        ent.flag = DFLAG_DISABLED;
      }

      //***Re-calculate material data from database***//

      if(data->ip.energyCoupl){
        plas_CalcMaterialData(data,ent.temp,flow.pressure);
      }

      if(ent.flag!=DFLAG_DISABLED){

        //***Perform neighbour element search routine***//

        plas_SearchSuccessive(data,&ent);

        //***Treat entities that left the domain after unsiccessful element search***//

        if(ent.flag==DFLAG_LEFT){

          plas_FindExitFace(data->fp.numBnd,data->fp.numDim,&ent,&facFound,&ibnd,&ifac);
          if(facFound){

            if(plasinterface_getWallBndFlag(ibnd)){

              //***Perform wall bounce***//

              plas_WallBounce(data->fp.numDim,data->rp.wallElasticity,&ent,ibnd,ifac);
              data->sd.bounce++;
              ent.flag = DFLAG_ENABLED;


            } else if(data->ip.perBnd && plasinterface_getPerBndFlag(ibnd)){

              //***Pass entity through periodic boundary***//

              plas_CalcBoundaryUnitNormal(data->fp.numDim,ibnd,ifac,unitVec);
              for(idim=0; idim<data->fp.numDim; idim++){
                ent.pos[idim] += plasinterface_getPerBndOffset(ibnd,idim);
              }
              ent.flag = DFLAG_PASS;
              data->sd.leftproc++;
              data->sd.periodic++;


            } else{

              //***Entity left computational domain through outlet***//

              ent.flag = DFLAG_DISABLED;
              data->sd.out++;

            }

          } else{

            //***Entity left computational domain through inter-process boundary***//

            ent.flag = DFLAG_PASS;
            data->sd.leftproc++;
          }
        }

        //***Stochastic collision model***//

        if(data->fp.numDim==3 && ent.flag==DFLAG_ENABLED && data->ip.collisionModel && data->pd[ent.node].volFrac>1e-3){
          plas_CollisionModel(data,&ent,data->pd[ent.node].numDens,dtLagr);
          plas_CheckNaN(data,&ent);
        }
      }

      if(ent.flag==DFLAG_ENABLED){

        //***Interpolate continuous phase velocity and pressure to new entity position***//

        plas_SetElementGeometry(data->fp.numDim,&ent);
        plas_Interpolate(data,&ent,&flow,(data->fp.dtEul-dtRemaining)/data->fp.dtEul);
        plas_CalcEntityCoefficients(data,&ent,&flow);

        //***Compute back-coupling terms to flow solver***//

        if(data->rp.flowType==FLOW_PARTIC || data->rp.flowType==FLOW_DROPLET){
          plas_CalcCouplingForcesParticle(data,&ent,&flow,dtLagr/data->fp.dtEul);
        } else if (data->rp.flowType==FLOW_BUBBLY){
          plas_CalcCouplingForcesBubble(data,&ent,&flow,dtLagr/data->fp.dtEul);
        }
      }

    } while(dtRemaining>0.0); /***End subiteration loop***/

    //***Re-calculate and update statistics***//

    if(ent.flag==DFLAG_ENABLED){
      plas_SetElementGeometry(data->fp.numDim,&ent);
      plas_Interpolate(data,&ent,&flow,1.0);
      plas_CalcEntityCoefficients(data,&ent,&flow);
      data->sd.reynoldsAvg += ent.reynolds;
      data->sd.nusseltAvg += ent.nusselt;
      data->sd.dtLagrAvg += (data->fp.dtEul/subIter);
      data->sd.subIterAvg += (double)subIter;
      avgctr++;
    }

    //***Put back entity information to global data structure***//

    data->ed[ient].flag = ent.flag;
    data->ed[ient].diameter = ent.diam;
    data->ed[ient].temperature = ent.temp;

    for(idim=0; idim<data->fp.numDim; idim++){
      data->ed[ient].velocity[idim] = ent.vel[idim];
      data->ed[ient].position[idim] = ent.pos[idim];
    }

    if(ent.flag==DFLAG_ENABLED){
      data->ed[ient].element = ent.elm;
      data->ed[ient].node = ent.node;
    } else{
      data->ed[ient].element = -1;
      data->ed[ient].node = -1;
    }

  } //***End trajectory loop***//

  plas_MpiBarrier();

  //***Pass entities between parallel processes***//

  if(plas_MpiGetNumProc()>1){plasinterface_screenOutput((char*) "Inter-process communication...\n");}
  plas_PassEntities(data);
  plas_MpiBarrier();

  //***Update counters***//

  plasinterface_screenOutput((char*) "Post-processing...\n");

  data->sd.enabled = plas_MpiAllSumInt(data->sd.enabled);
  data->sd.in = plas_MpiAllSumInt(data->sd.in);
  data->sd.out = plas_MpiAllSumInt(data->sd.out);
  data->sd.bounce = plas_MpiAllSumInt(data->sd.bounce);
  data->sd.periodic = plas_MpiAllSumInt(data->sd.periodic);
  data->sd.passed = plas_MpiAllSumInt(data->sd.passed);
  data->sd.leftproc = plas_MpiAllSumInt(data->sd.leftproc);
  data->sd.coll = plas_MpiAllSumInt(data->sd.coll);
  data->sd.coalesc = plas_MpiAllSumInt(data->sd.coalesc);
  data->sd.lost = plas_MpiAllSumInt(data->sd.lost);
  data->sd.enabled -= data->sd.out;
  data->sd.lost += data->sd.leftproc-data->sd.passed;
  data->sd.enabled -= data->sd.lost;
  plas_MpiBarrier();

  //***Update statistics***//

  if(data->sd.enabled>0){
    data->sd.reynoldsAvg = plas_MpiAllSumDouble(data->sd.reynoldsAvg)/plas_MpiAllSumInt(avgctr);
    data->sd.nusseltAvg = plas_MpiAllSumDouble(data->sd.nusseltAvg)/plas_MpiAllSumInt(avgctr);
    data->sd.dtLagrAvg = plas_MpiAllSumDouble(data->sd.dtLagrAvg)/plas_MpiAllSumInt(avgctr);
    data->sd.subIterAvg = plas_MpiAllSumDouble(data->sd.subIterAvg)/plas_MpiAllSumInt(avgctr);
  } else{
    data->sd.reynoldsAvg = 0.0;
    data->sd.nusseltAvg = 0.0;
    data->sd.dtLagrAvg = 0.0;
    data->sd.subIterAvg = 0.0;
  }
  plas_MpiBarrier();

  //**Compute dispersed phase cellwise data***//

  plas_CalcCellwiseData(data);
  plas_MpiBarrier();

  //***Write PLaS output to files***//

  plas_WriteStatsFile(data,data->ip.writeStatsFilename,data->fp.iter,data->fp.time);
  if (data->fp.writeOutput) {
    plas_WriteTecplotFile(data,data->ip.writeTecplotFilename,data->fp.iter,data->fp.time);
  }
  plas_MpiBarrier();

  //***Free dynamic memory***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(data->fp.numDim,&flow);
  plas_MpiBarrier();
}

