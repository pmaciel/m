
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#ifdef MPI
#include <mpi.h>
#endif
#include "plas.h"


void plas::initialize()
{
  int ient,inod;
  FILE *inpFile;

  //***Flow solver parameters that have to be set only once***//

  setFlowSolverParamOnInit(&fp);

  //***Set partitioning data structure in case of parallel computation***//

  setPartitioningData(&(fp.part));

  //***Read parameters from input file and broadcast them***//

  if (!plas_MpiGetRank()) {
    plas_ReadParameters();
  }
  plas_MpiBarrier();
  plas_BroadcastParameters();
  plas_MpiBarrier();

  //***Allocate memory for dispersed phase data***//

  ed = new PLAS_ENTITY_DATA[ip.numMaxEnt];
  for (ient=0; ient<ip.numMaxEnt; ++ient) {
    ed[ient].flag        = DFLAG_DISABLED;
    ed[ient].position    = new double[fp.numDim];
    ed[ient].velocity    = new double[fp.numDim];
    ed[ient].velocityOld = new double[fp.numDim];
    ed[ient].node        = -1;
    ed[ient].element     = -1;
  }

  pd = new PLAS_PHASE_DATA[fp.numNod];
  for(inod=0; inod<fp.numNod; inod++){
    pd[inod].volFracDx = new double[fp.numDim];
    pd[inod].avgVel    = new double[fp.numDim];
    pd[inod].stdVel    = new double[fp.numDim];
    pd[inod].dispForce = new double[fp.numUnk];
  }

  //***Initialize fixed runtime parameters***//

  rp.lagrTimeFactor = 0.3;
  rp.errTol         = 1e-6;
  rp.wallElasticity = 1.;

  //***Initialize variable runtime parameters***//

  if (ip.numProdDom>0)
    rp.massResid = new double[ip.numProdDom];

  //***Initialize material data from database***//

  plas_CalcMaterialData(ip.iniTempDisp,getPressure(0));

  //***Initial entity distribution***//

  /* set output files */
  inpFile = fopen(ip.writeTecplotFilename,"r");
  if (inpFile==NULL || !ip.restart) {
    ip.restart = 0;
    plas_CreateTecplotFile(ip.writeTecplotFilename);
  }
  else
    fclose(inpFile);
  inpFile = fopen(ip.writeStatsFilename,"r");
  if (inpFile==NULL) {
    plas_CreateStatsFile(ip.writeStatsFilename);
  }
  else
    fclose(inpFile);

  /* loaded or random distribution */
  if (ip.restart)
    plas_LoadInitialDistribution(ip.writeTecplotFilename);
  else if (ip.numIniEnt>0)
    plas_RandomInitialDistribution();

  plas_MpiBarrier();
}


void plas::run()
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  char errMessage[100];
  int idx,ient,jent,inod,idim,iunk,ibnd,ifac,subIter,facFound,avgctr;
  double cellSize,minCellSize,entVel,dtLagr,dtRemaining,unitVec[fp.numDim];

  //***Flow solver parameters that have to be set at every time step***//

  screenOutput((char*) "Initialization...\n");
  setFlowSolverParamOnTimeStep(&(fp));

  //***Allocate dynamic memory***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);

  //***Initializations***//

  sd.enabled     = 0;
  sd.in          = 0;
  sd.out         = 0;
  sd.bounce      = 0;
  sd.periodic    = 0;
  sd.passed      = 0;
  sd.lost        = 0;
  sd.coll        = 0;
  sd.coalesc     = 0;
  sd.leftproc    = 0;
  sd.dtLagrAvg   = 0.;
  sd.reynoldsAvg = 0.;
  sd.nusseltAvg  = 0.;
  sd.subIterAvg  = 0.;
  avgctr         = 0;
  minCellSize = pow(fp.minElmVolume,(1.0/fp.numDim));

  if(fp.iter==1){
    plas_CalcCellwiseData();
  }
  plas_MpiBarrier();

  for(inod=0; inod<fp.numNod; inod++){
    for(iunk=0; iunk<fp.numUnk; iunk++){
      pd[inod].dispForce[iunk] = 0.0;
    }
  }

  //***Sort and count active entities***//

  idx = ip.numMaxEnt-1;
  for(ient=0; ient<ip.numMaxEnt; ient++){
    if(ed[ient].flag==DFLAG_ENABLED || ed[ient].flag==DFLAG_CREATED){
      sd.enabled++;
      continue;
    }
    for(jent=idx; jent>ient; jent--){
      if(ed[jent].flag==DFLAG_DISABLED){
        idx--;
        continue;
      }
      ed[ient].element = ed[jent].element;
      ed[ient].node    = ed[jent].node;
      for(idim=0; idim<fp.numDim; idim++){
        ed[ient].position[idim] = ed[jent].position[idim];
        ed[ient].velocity[idim] = ed[jent].velocity[idim];
      }
      ed[ient].diameter = ed[jent].diameter;
      ed[ient].flag = ed[jent].flag;
      ed[jent].flag = DFLAG_DISABLED;
      idx = jent-1;
      sd.enabled++;
      break;
    }
  }
  plas_MpiBarrier();

  //***Entity production***//

  if (ip.numProdDom>0){
    screenOutput((char*) "Imposing new entities in production domains...\n");
    plas_ImposeProductionDomains();
  }
  plas_MpiBarrier();

  if(plas_MpiAllSumInt(numExtEnt)>0){
    screenOutput((char*) "Imposing externally generated entities...\n");
    plas_ImposeExternal();
  }
  plas_MpiBarrier();

  //***Loop over dispersed entities to update trajectories***//

  screenOutput((char*) "Updating trajectories...\n");
  for(ient=0; ient<sd.enabled; ient++){

    //***Get entity information from global data structure***//

    ent.flag = ed[ient].flag;
    if(ent.flag!=DFLAG_ENABLED && ent.flag!=DFLAG_CREATED){
      sprintf(errMessage,"Entity %i had bad flag %d.",ient,ent.flag);
      plas_TerminateOnError(errMessage);
    }

    ent.elm  = ed[ient].element;
    ent.node = ed[ient].node;
    ent.diam = ed[ient].diameter;
    ent.temp = ed[ient].temperature;

    for(idim=0; idim<fp.numDim; idim++){
      ent.pos[idim] = ed[ient].position[idim];
      ent.vel[idim] = ed[ient].velocity[idim];
      ed[ient].velocityOld[idim] = ed[ient].velocity[idim];
    }

    //***Subiterations of the trajectory increment***//

    subIter = 0;
    dtRemaining = fp.dtEul;

    do{

      //***Determine Lagrangian time step size***//

      entVel = plas_CalcVectorLength(fp.numDim,ent.vel);
      if(entVel<rp.errTol){
        entVel = rp.errTol;
      }

      cellSize = pow(getElmVol(ent.elm),(1.0/fp.numDim));
      if(cellSize<rp.errTol){
        cellSize = rp.errTol;
      }

      dtLagr = rp.lagrTimeFactor*(cellSize/entVel);
      if(dtLagr>dtRemaining){
        dtLagr = dtRemaining;
      }

      dtRemaining -= dtLagr;
      subIter++;

      //***Entity production***//

      if(ent.flag==DFLAG_CREATED){
        if(plas_RandomDouble()<=((fp.dtEul-dtRemaining)/fp.dtEul)){
          ent.flag = DFLAG_ENABLED;
        }
      }

      if(ent.flag!=DFLAG_ENABLED){
        continue;
      }

      //***Interpolate continuous phase velocity and pressure at entity position***//

      plas_SetElementGeometry(fp.numDim,&ent);
      plas_Interpolate(&ent,&flow,(fp.dtEul-dtRemaining)/fp.dtEul);
      plas_CalcEntityCoefficients(&ent,&flow);

      //***Solve Lagrangian trajectory equation***//

      plas_CalcTrajectory(&ent,&flow,dtLagr);
      plas_CheckNaN(&ent);

      //***Check for evaporated, burned or collapsed entity***//

      if(ent.diam<rp.errTol){
        ent.flag = DFLAG_DISABLED;
      }

      //***Re-calculate material data from database***//

      if(ip.energyCoupl){
        plas_CalcMaterialData(ent.temp,flow.pressure);
      }

      if(ent.flag!=DFLAG_DISABLED){

        //***Perform neighbour element search routine***//

        plas_SearchSuccessive(&ent);

        //***Treat entities that left the domain after unsiccessful element search***//

        if(ent.flag==DFLAG_LEFT){

          plas_FindExitFace(fp.numBnd,fp.numDim,&ent,&facFound,&ibnd,&ifac);
          if(facFound){

            if(getWallBndFlag(ibnd)){

              //***Perform wall bounce***//

              plas_WallBounce(fp.numDim,rp.wallElasticity,&ent,ibnd,ifac);
              sd.bounce++;
              ent.flag = DFLAG_ENABLED;


            } else if(ip.perBnd && getPerBndFlag(ibnd)){

              //***Pass entity through periodic boundary***//

              plas_CalcBoundaryUnitNormal(fp.numDim,ibnd,ifac,unitVec);
              for(idim=0; idim<fp.numDim; idim++){
                ent.pos[idim] += getPerBndOffset(ibnd,idim);
              }
              ent.flag = DFLAG_PASS;
              sd.leftproc++;
              sd.periodic++;


            } else{

              //***Entity left computational domain through outlet***//

              ent.flag = DFLAG_DISABLED;
              sd.out++;

            }

          } else{

            //***Entity left computational domain through inter-process boundary***//

            ent.flag = DFLAG_PASS;
            sd.leftproc++;
          }
        }

        //***Stochastic collision model***//

        if(fp.numDim==3 && ent.flag==DFLAG_ENABLED && ip.collisionModel && pd[ent.node].volFrac>1e-3){
          plas_CollisionModel(&ent,pd[ent.node].numDens,dtLagr);
          plas_CheckNaN(&ent);
        }
      }

      if(ent.flag==DFLAG_ENABLED){

        //***Interpolate continuous phase velocity and pressure to new entity position***//

        plas_SetElementGeometry(fp.numDim,&ent);
        plas_Interpolate(&ent,&flow,(fp.dtEul-dtRemaining)/fp.dtEul);
        plas_CalcEntityCoefficients(&ent,&flow);

        //***Compute back-coupling terms to flow solver***//

        if(rp.flowType==FLOW_PARTIC || rp.flowType==FLOW_DROPLET){
          plas_CalcCouplingForcesParticle(&ent,&flow,dtLagr/fp.dtEul);
        } else if (rp.flowType==FLOW_BUBBLY){
          plas_CalcCouplingForcesBubble(&ent,&flow,dtLagr/fp.dtEul);
        }
      }

    } while(dtRemaining>0.0); /***End subiteration loop***/

    //***Re-calculate and update statistics***//

    if(ent.flag==DFLAG_ENABLED){
      plas_SetElementGeometry(fp.numDim,&ent);
      plas_Interpolate(&ent,&flow,1.0);
      plas_CalcEntityCoefficients(&ent,&flow);
      sd.reynoldsAvg += ent.reynolds;
      sd.nusseltAvg  += ent.nusselt;
      sd.dtLagrAvg   += (fp.dtEul/subIter);
      sd.subIterAvg  += (double)subIter;
      avgctr++;
    }

    //***Put back entity information to global data structure***//

    ed[ient].flag = ent.flag;
    ed[ient].diameter = ent.diam;
    ed[ient].temperature = ent.temp;

    for(idim=0; idim<fp.numDim; idim++){
      ed[ient].velocity[idim] = ent.vel[idim];
      ed[ient].position[idim] = ent.pos[idim];
    }

    if(ent.flag==DFLAG_ENABLED){
      ed[ient].element = ent.elm;
      ed[ient].node = ent.node;
    } else{
      ed[ient].element = -1;
      ed[ient].node = -1;
    }

  } //***End trajectory loop***//

  plas_MpiBarrier();

  //***Pass entities between parallel processes***//

  if(plas_MpiGetNumProc()>1){screenOutput((char*) "Inter-process communication...\n");}
  plas_PassEntities();
  plas_MpiBarrier();

  //***Update counters***//

  screenOutput((char*) "Post-processing...\n");

  sd.enabled   = plas_MpiAllSumInt(sd.enabled);
  sd.in        = plas_MpiAllSumInt(sd.in);
  sd.out       = plas_MpiAllSumInt(sd.out);
  sd.bounce    = plas_MpiAllSumInt(sd.bounce);
  sd.periodic  = plas_MpiAllSumInt(sd.periodic);
  sd.passed    = plas_MpiAllSumInt(sd.passed);
  sd.leftproc  = plas_MpiAllSumInt(sd.leftproc);
  sd.coll      = plas_MpiAllSumInt(sd.coll);
  sd.coalesc   = plas_MpiAllSumInt(sd.coalesc);
  sd.lost      = plas_MpiAllSumInt(sd.lost);
  sd.enabled  -= sd.out;
  sd.lost     += sd.leftproc-sd.passed;
  sd.enabled  -= sd.lost;
  plas_MpiBarrier();

  //***Update statistics***//

  if(sd.enabled>0){
    sd.reynoldsAvg = plas_MpiAllSumDouble(sd.reynoldsAvg)/plas_MpiAllSumInt(avgctr);
    sd.nusseltAvg  = plas_MpiAllSumDouble(sd.nusseltAvg)/plas_MpiAllSumInt(avgctr);
    sd.dtLagrAvg   = plas_MpiAllSumDouble(sd.dtLagrAvg)/plas_MpiAllSumInt(avgctr);
    sd.subIterAvg  = plas_MpiAllSumDouble(sd.subIterAvg)/plas_MpiAllSumInt(avgctr);
  } else{
    sd.reynoldsAvg = 0.;
    sd.nusseltAvg  = 0.;
    sd.dtLagrAvg   = 0.;
    sd.subIterAvg  = 0.;
  }
  plas_MpiBarrier();

  //**Compute dispersed phase cellwise data***//

  plas_CalcCellwiseData();
  plas_MpiBarrier();

  //***Write PLaS output to files***//

  plas_WriteStatsFile(ip.writeStatsFilename,fp.iter,fp.time);
  if (fp.writeOutput) {
    plas_WriteTecplotFile(ip.writeTecplotFilename,fp.iter,fp.time);
  }
  plas_MpiBarrier();

  //***Free dynamic memory***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(fp.numDim,&flow);
  plas_MpiBarrier();
}


plas::~plas()
{
  // free dynamic memory

  for (int i=0; i<ip.numMaxEnt; ++i) {
    delete[] ed[i].position;
    delete[] ed[i].velocity;
    delete[] ed[i].velocityOld;
  }
  delete[] ed;

  for (int i=0; i<fp.numNod; ++i) {
    delete[] pd[i].volFracDx;
    delete[] pd[i].avgVel;
    delete[] pd[i].stdVel;
    delete[] pd[i].dispForce;
  }
  delete[] pd;

  if (ip.numProdDom>0) {
    for (int i=0; i<ip.numProdDom; ++i)
      delete[] ip.prodParam[i];
    delete[] ip.prodParam;
    delete[] ip.prodDom;
    delete[] ip.massFluxes;
    delete[] rp.massResid;
  }

  if (fp.part.numNodPairs) {
    delete[] fp.part.sendRank;
    delete[] fp.part.sendNodeIdx;
    delete[] fp.part.recvRank;
    delete[] fp.part.recvNodeIdx;
  }
}


void plas::plas_AllocateLocalEntityVar(int numDim, LOCAL_ENTITY_VARIABLES *ent)
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


void plas::plas_AllocateLocalFlowVar(int numDim, LOCAL_FLOW_VARIABLES *flow)
{
  flow->vel   = new double[numDim];
  flow->velDt = new double[numDim];
  flow->velDx = new double*[numDim];
  for (int idim=0; idim<numDim; ++idim)
    flow->velDx[idim] = new double[numDim];
  flow->vort = new double[numDim];
}


void plas::plas_BroadcastParameters()
{
  int irank = plas_MpiGetRank();
  int i;

  //***Broadcast static data***//

  plas_MpiBroadcastInt(&ip.numMaxEnt,1,0);
  plas_MpiBroadcastInt(&ip.numIniEnt,1,0);
  plas_MpiBroadcastInt(&ip.numProdDom,1,0);
  plas_MpiBroadcastInt(&ip.iniDiamType,1,0);
  plas_MpiBroadcastDouble(&ip.iniDiam,1,0);
  plas_MpiBroadcastDouble(&ip.iniDiamStd,1,0);
  plas_MpiBroadcastDouble(ip.iniVel,3,0);
  plas_MpiBroadcastDouble(&ip.iniTempDisp,1,0);
  plas_MpiBroadcastInt(&ip.material,1,0);
  plas_MpiBroadcastInt(&rp.flowType,1,0);
  plas_MpiBroadcastInt(&ip.momentumCoupl,1,0);
  plas_MpiBroadcastInt(&ip.volfracCoupl,1,0);
  plas_MpiBroadcastInt(&ip.energyCoupl,1,0);
  plas_MpiBroadcastInt(&ip.collisionModel,1,0);
  plas_MpiBroadcastInt(&ip.liftForce,1,0);
  plas_MpiBroadcastInt(&ip.evapModel,1,0);
  plas_MpiBroadcastInt(&ip.perBnd,1,0);
  plas_MpiBroadcastDouble(ip.gravVec,3,0);
  plas_MpiBroadcastInt(&ip.restart,1,0);

  /*
   * ip.writeTecplotFilename and ip.writeStatsFilename don't need
   * broadcasting because only rank 0 writes them, and ip.confFilename
   * doesn't need reading (again)
   */

  //***Broadcast allocatable data***//

  if (ip.numProdDom>0) {
    if (irank!=0)
      ip.prodDom = new int[ip.numProdDom];
    plas_MpiBroadcastInt(ip.prodDom,ip.numProdDom,0);
    if (irank!=0) {
      ip.prodParam = new double*[ip.numProdDom];
      for (i=0; i<ip.numProdDom; ++i)
        ip.prodParam[i] = new double[6];
    }
    for (i=0; i<ip.numProdDom; ++i)
      plas_MpiBroadcastDouble(ip.prodParam[i],6,0);
    if (irank!=0)
      ip.massFluxes = new double[ip.numProdDom];
    plas_MpiBroadcastDouble(ip.massFluxes,ip.numProdDom,0);
  }
}


void plas::plas_CalcBackCoupling(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double *force, double tFactor)
{
  double limitFrac = 0.95;
  int idim,jdim,inod;
  double contVol,iFrac,jFrac,rhsTerm,impFac[ent->edata.numElmNodes];

  //***Volume fraction trigger***//

  if(ip.volfracCoupl){
    if(pd[ent->node].volFrac>limitFrac){
      iFrac = limitFrac;
    } else{
      iFrac = pd[ent->node].volFrac;
    }
  } else{
    iFrac = 0.0;
  }

  if(pd[ent->node].volFrac>1.0){
    jFrac = 1.0/pd[ent->node].volFrac;
  } else{
    jFrac = 1.0;
  }

  //***Compute back-coupling force contribution (m/s^2) for fluid momentum equation (negative sign)***//

  if(ip.momentumCoupl==FORCE_PIC){

    inod = ent->node;
    contVol = getNodVol(ent->node);
    for(idim=0; idim<fp.numDim; idim++){
      pd[ent->node].dispForce[idim+1] -= tFactor*force[idim]*jFrac/((1.0-iFrac)*contVol*fp.rhoCont);
    }

  } else if(ip.momentumCoupl==FORCE_PROJ){

    plas_CalcNodeImpactFactors(ent,impFac);
    for(idim=0; idim<ent->edata.numElmNodes; idim++){
      inod = ent->edata.elmNodes[idim];
      if(fp.flowSolver==FLOWSOLVER_SFELES){inod %= fp.numNod;}
      contVol = getNodVol(inod);
      for(jdim=0; jdim<fp.numDim; jdim++){
        pd[inod].dispForce[jdim+1] -= tFactor*impFac[idim]*force[jdim]*jFrac/((1.0-iFrac)*contVol*fp.rhoCont);
      }
    }
  }

  //***Compute back-coupling term for fluid continuity equation***//

  if(ip.volfracCoupl){

    inod = ent->node;
    rhsTerm = pd[inod].volFracDt;
    for(idim=0; idim<fp.numDim; idim++){
      rhsTerm += flow->vel[idim]*pd[inod].volFracDx[idim];
    }

    pd[inod].dispForce[0] += tFactor*jFrac*rhsTerm/(1.0-iFrac);
    for(idim=1; idim<fp.numDim; idim++){
      pd[inod].dispForce[idim] += tFactor*jFrac*rhsTerm*flow->vel[idim]/(1.0-iFrac);
    }
  }
}


void plas::plas_CalcBoundaryUnitNormal(int numDim, int ibnd, int ifac, double *unitVec)
{
  int idim;
  double length;

  //***Get normal vector of a boundary face***//

  for(idim=0; idim<numDim; idim++){
    unitVec[idim] = getBndFaceNormComp(ibnd,ifac,idim);
  }

  //***Divide normal vector components by vector length***//

  length = plas_CalcVectorLength(numDim,unitVec);
  for(idim=0; idim<numDim; idim++){
    unitVec[idim] /= length;
  }
}


void plas::plas_CalcCellwiseData()
{
  int irank = plas_MpiGetRank();
  int ient,inod,ielm,idim,jdim,itype,numElmNodes=0,correctNumDens[fp.part.numNodPairs];
  double ivol=0,xi,yi,xixi,xiyi,volFracOld[fp.numNod];
  double correctDiam[fp.part.numNodPairs],correctRespTime[fp.part.numNodPairs];
  double correctVel[fp.numDim*fp.part.numNodPairs];
  double correctStdDiam[fp.part.numNodPairs],correctStdVel[fp.numDim*fp.part.numNodPairs];
  double correctVolFracDx[fp.numDim*fp.part.numNodPairs];

  //***Initialize cellwise secondary phase data***//

  for(inod=0; inod<fp.numNod; inod++){
    pd[inod].numDens = 0;
    volFracOld[inod] = pd[inod].volFrac;
    pd[inod].volFrac = 0.0;
    pd[inod].volFracDt = 0.0;
    pd[inod].avgDiam = 0.0;
    pd[inod].stdDiam = 0.0;
    pd[inod].avgRespTime = 0.0;
    for(idim=0; idim<fp.numDim; idim++){
      pd[inod].volFracDx[idim] = 0.0;
      pd[inod].avgVel[idim] = 0.0;
      pd[inod].stdVel[idim] = 0.0;
    }
  }

  //***Assemble cellwise data by looping over all active entities***//

  for(ient=0; ient<ip.numMaxEnt; ient++){
    if(ed[ient].flag==DFLAG_ENABLED){
      inod = ed[ient].node;
      if(fp.flowSolver==FLOWSOLVER_SFELES){inod %= fp.numNod;}
      pd[inod].numDens++;
      if(fp.numDim==2){
        ivol = PI*ed[ient].diameter*ed[ient].diameter/4.0;
      } else if(fp.numDim==3){
        ivol = PI*ed[ient].diameter*ed[ient].diameter*ed[ient].diameter/6.0;
      }
      pd[inod].volFrac += ivol/getNodVol(inod);
      pd[inod].avgDiam += ed[ient].diameter;
      pd[inod].avgRespTime += (2.0*md.rhoDisp+fp.rhoCont)
                                   *ed[ient].diameter*ed[ient].diameter/(24.0*fp.muCont);
      for(idim=0; idim<fp.numDim; idim++){
        pd[inod].avgVel[idim] += ed[ient].velocity[idim];
      }
    }
  }

  //***Divide averaged quantities by number density***//

  for(inod=0; inod<fp.numNod; inod++){
    if(pd[inod].numDens>0){
      pd[inod].avgDiam /= pd[inod].numDens;
      pd[inod].avgRespTime /= pd[inod].numDens;
      for(idim=0; idim<fp.numDim; idim++){
        pd[inod].avgVel[idim] /= pd[inod].numDens;
      }
    }
  }

  //***Compute standard deviations of averaged quantities***//

  if(ip.collisionModel){

    for(ient=0; ient<ip.numMaxEnt; ient++){
      if(ed[ient].flag==DFLAG_ENABLED){
        inod = ed[ient].node;
        if(fp.flowSolver==FLOWSOLVER_SFELES){inod %= fp.numNod;}

        pd[inod].stdDiam += pow(ed[ient].diameter-pd[inod].avgDiam,2.0);
        for(idim=0; idim<fp.numDim; idim++){
          pd[inod].stdVel[idim] += pow(ed[ient].velocity[idim]-pd[inod].avgVel[idim],2.0);
        }
      }
    }

    for(inod=0; inod<fp.numNod; inod++){
      if(pd[inod].numDens>0){
        pd[inod].stdDiam = sqrt(pd[inod].stdDiam/pd[inod].numDens);
        for(idim=0; idim<fp.numDim; idim++){
          pd[inod].stdVel[idim] = sqrt(pd[inod].stdVel[idim]/pd[inod].numDens);
        }
      }
    }
  }

  //***Calculate volume fraction derivatrives in time and space***//

  if(ip.volfracCoupl){

    for(inod=0; inod<fp.numNod; inod++){
      pd[inod].volFracDt = (pd[inod].volFrac-volFracOld[inod])/fp.dtEul;
    }

    for(ielm=0; ielm<fp.numElm; ielm++){
      for(idim=0; idim<fp.numDim; idim++){

        itype = getElementType(ielm);
        if(itype==ELM_SIMPLEX){
          numElmNodes = fp.numDim+1;
        } else if(itype==ELM_PRISM){
          numElmNodes = 6;
        } else if(itype==ELM_QUAD){
          numElmNodes = 4;
        } else if(itype==ELM_HEX){
          numElmNodes = 8;
        } else if(itype==ELM_PYRAMID){
          numElmNodes = 5;
        }

        xi = yi = xixi = xiyi = 0.0;
        for(jdim=0; jdim<numElmNodes; jdim++){
          inod = getElmNode(ielm,jdim);
          xi += getNodCoord(inod,idim);
          yi += pd[inod].volFrac;
          xixi += getNodCoord(inod,idim)*getNodCoord(inod,idim);
          xiyi += getNodCoord(inod,idim)*pd[inod].volFrac;
        }
        pd[inod].volFracDx[idim] = (numElmNodes*xiyi-xi*yi)/(numElmNodes*xixi-xi*xi);
      }
    }
  }

  //***Gather phase data on partition boundaries***//

  for(inod=0; inod<fp.part.numNodPairs; inod++){

    correctNumDens[inod] = 0;
    correctDiam[inod] = 0.0;
    correctStdDiam[inod] = 0.0;
    correctRespTime[inod] = 0.0;
    for(idim=0; idim<fp.numDim; idim++){
      correctVel[fp.numDim*inod+idim] = 0.0;
      correctStdVel[fp.numDim*inod+idim] = 0.0;
      correctVolFracDx[fp.numDim*inod+idim] = 0.0;
    }

    if(fp.part.sendRank[inod]==irank){
      correctNumDens[inod] += pd[fp.part.sendNodeIdx[inod]].numDens;
      correctDiam[inod] += pd[fp.part.sendNodeIdx[inod]].avgDiam*correctNumDens[inod];
      correctStdDiam[inod] += pd[fp.part.sendNodeIdx[inod]].stdDiam*correctNumDens[inod];
      correctRespTime[inod] += pd[fp.part.sendNodeIdx[inod]].avgRespTime*correctNumDens[inod];
      for(idim=0; idim<fp.numDim; idim++){
        correctVel[fp.numDim*inod+idim] += pd[fp.part.sendNodeIdx[inod]].avgVel[idim]*correctNumDens[inod];
        correctStdVel[fp.numDim*inod+idim] += pd[fp.part.sendNodeIdx[inod]].stdVel[idim]*correctNumDens[inod];
        correctVolFracDx[fp.numDim*inod+idim] += pd[fp.part.sendNodeIdx[inod]].volFracDx[idim]*correctNumDens[inod];
      }
    }

    if(fp.part.recvRank[inod]==irank){
      correctNumDens[inod] += pd[fp.part.recvNodeIdx[inod]].numDens;
      correctDiam[inod] += pd[fp.part.recvNodeIdx[inod]].avgDiam*correctNumDens[inod];
      correctStdDiam[inod] += pd[fp.part.recvNodeIdx[inod]].stdDiam*correctNumDens[inod];
      correctRespTime[inod] += pd[fp.part.recvNodeIdx[inod]].avgRespTime*correctNumDens[inod];
      for(idim=0; idim<fp.numDim; idim++){
        correctVel[fp.numDim*inod+idim] += pd[fp.part.recvNodeIdx[inod]].avgVel[idim]*correctNumDens[inod];
        correctStdVel[fp.numDim*inod+idim] += pd[fp.part.recvNodeIdx[inod]].stdVel[idim]*correctNumDens[inod];
        correctVolFracDx[fp.numDim*inod+idim] += pd[fp.part.recvNodeIdx[inod]].volFracDx[idim]*correctNumDens[inod];
      }
    }
  }

  //***Sum phase data on partition boundaries***//

  if(fp.part.numNodPairs){
    plas_MpiAllSumIntArray(correctNumDens,fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctDiam,fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctStdDiam,fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctRespTime,fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctVel,fp.numDim*fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctStdVel,fp.numDim*fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctVolFracDx,fp.numDim*fp.part.numNodPairs);
  }

  //***Correct phase data on partition boundaries***//

  for(inod=0; inod<fp.part.numNodPairs; inod++){

    if(fp.part.sendRank[inod]==irank){
      pd[fp.part.sendNodeIdx[inod]].numDens = correctNumDens[inod];
      if(correctNumDens[inod]>0){
        pd[fp.part.sendNodeIdx[inod]].avgDiam = correctDiam[inod]/correctNumDens[inod];
        pd[fp.part.sendNodeIdx[inod]].stdDiam = correctStdDiam[inod]/correctNumDens[inod];
        pd[fp.part.sendNodeIdx[inod]].avgRespTime = correctRespTime[inod]/correctNumDens[inod];
        for(idim=0; idim<fp.numDim; idim++){
          pd[fp.part.sendNodeIdx[inod]].avgVel[idim] = correctVel[fp.numDim*inod+idim]/correctNumDens[inod];
          pd[fp.part.sendNodeIdx[inod]].stdVel[idim] = correctStdVel[fp.numDim*inod+idim]/correctNumDens[inod];
          pd[fp.part.sendNodeIdx[inod]].volFracDx[idim] = correctVolFracDx[fp.numDim*inod+idim]/correctNumDens[inod];
        }
      }
    }

    if(fp.part.recvRank[inod]==irank){
      pd[fp.part.recvNodeIdx[inod]].numDens = correctNumDens[inod];
      if(correctNumDens[inod]>0){
        pd[fp.part.recvNodeIdx[inod]].avgDiam = correctDiam[inod]/correctNumDens[inod];
        pd[fp.part.recvNodeIdx[inod]].stdDiam = correctStdDiam[inod]/correctNumDens[inod];
        pd[fp.part.recvNodeIdx[inod]].avgRespTime = correctRespTime[inod]/correctNumDens[inod];
        for(idim=0; idim<fp.numDim; idim++){
          pd[fp.part.recvNodeIdx[inod]].avgVel[idim] = correctVel[fp.numDim*inod+idim]/correctNumDens[inod];
          pd[fp.part.recvNodeIdx[inod]].stdVel[idim] = correctStdVel[fp.numDim*inod+idim]/correctNumDens[inod];
          pd[fp.part.recvNodeIdx[inod]].volFracDx[idim] = correctVolFracDx[fp.numDim*inod+idim]/correctNumDens[inod];
        }
      }
    }
  }
}


double plas::plas_CalcConcInterf(double pressBubble)
{
  return (pressBubble*md.HeDisp);
}


void plas::plas_CalcCouplingForcesBubble(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor)
{
  int idim,jdim;
  double densityRatio,bubbleVol,bubbleMass,vmCoeff;
  double uxw[fp.numDim],dudt[fp.numDim],dvdt[fp.numDim],entForce[fp.numDim];

  //***Compute momentum coupling forces***//

  if(ip.momentumCoupl){

    bubbleVol = PI*ent->diam*ent->diam*ent->diam/6.0;
    bubbleMass = md.rhoDisp*bubbleVol;
    densityRatio = (fp.rhoCont/md.rhoDisp);

    if(fp.numDim==2){
      uxw[0] = 0.0;
      uxw[1] = 0.0;
    } else if(fp.numDim==3){
      uxw[0] = ent->relVel[1]*flow->vort[2]-ent->relVel[2]*flow->vort[1];
      uxw[1] = ent->relVel[2]*flow->vort[0]-ent->relVel[0]*flow->vort[2];
      uxw[2] = ent->relVel[0]*flow->vort[1]-ent->relVel[1]*flow->vort[0];
    }

    vmCoeff = 0.5*(1.0+2.786*pd[ent->node].volFrac);

    for(idim=0; idim<fp.numDim; idim++){
      dudt[idim] = flow->velDt[idim];
      for(jdim=0; jdim<fp.numDim; jdim++){
        dudt[idim] += flow->vel[jdim]*flow->velDx[idim][jdim];
      }
      dvdt[idim] = (ent->vel[idim]-ent->velOld[idim])/fp.dtEul;
    }

    //**Compute surface force on particle (positive sign)***//

    for(idim=0; idim<fp.numDim; idim++){
      entForce[idim] = bubbleMass*(3.0/4.0)*(ent->dragCoeff/ent->diam)*densityRatio*ent->normVel*ent->relVel[idim]
                     + bubbleMass*ent->liftCoeff*densityRatio*uxw[idim]
                     + bubbleMass*vmCoeff*densityRatio*(dudt[idim]-dvdt[idim]);
    }
  }

  //***Compute back coupling terms***//

  plas_CalcBackCoupling(ent,flow,entForce,tFactor);

}


void plas::plas_CalcCouplingForcesParticle(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor)
{
  int idim;
  double densityRatio,particleVol,particleMass;
  double uxw[fp.numDim],entForce[fp.numDim];

  //***Compute momentum coupling forces***//

  if(ip.momentumCoupl){

    particleVol = PI*ent->diam*ent->diam*ent->diam/6.0;
    particleMass = md.rhoDisp*particleVol;
    densityRatio = (fp.rhoCont/md.rhoDisp);

    if(fp.numDim==2){
      uxw[0] = 0.0;
      uxw[1] = 0.0;
    } else if(fp.numDim==3){
      uxw[0] = ent->relVel[1]*flow->vort[2]-ent->relVel[2]*flow->vort[1];
      uxw[1] = ent->relVel[2]*flow->vort[0]-ent->relVel[0]*flow->vort[2];
      uxw[2] = ent->relVel[0]*flow->vort[1]-ent->relVel[1]*flow->vort[0];
    }

    //**Compute surface force on particle (positive sign)***//

    for(idim=0; idim<fp.numDim; idim++){
      entForce[idim] = particleMass*(3.0/4.0)*(ent->dragCoeff/ent->diam)*densityRatio*ent->normVel*ent->relVel[idim]
                     + particleMass*ent->liftCoeff*densityRatio*uxw[idim];
    }
  }

  //***Compute back coupling terms***//

  plas_CalcBackCoupling(ent,flow,entForce,tFactor);

}


void plas::plas_CalcCrossProduct_3D(double *value, double *a, double *b)
{
  value[0] = a[1]*b[2] - a[2]*b[1];
  value[1] = a[2]*b[0] - a[0]*b[2];
  value[2] = a[0]*b[1] - a[1]*b[0];
}


double plas::plas_CalcDispReynolds(double viscosity, double diameter, double normVel)
{
  double reynolds = diameter*normVel/viscosity;

  if(reynolds<1e-4){reynolds = 1e-4;}

  return reynolds;
}


double plas::plas_CalcDragCoeff(int flowType, double reynolds)
{
  double dragCoeff = 0.;

  if(flowType==FLOW_BUBBLY){

    //***Empirical drag coefficient correlations for a fluid sphere***//

    if(reynolds<1.5){
      dragCoeff = 16.0/reynolds;
    } else if(reynolds>=1.5 && reynolds<80.0){
      dragCoeff = 14.9/pow(reynolds,0.78);
    } else if(reynolds>=80.0 && reynolds<1530.0){
      dragCoeff = (49.9/reynolds)*(1.0-(2.21/pow(reynolds,0.5))) + (1.17e-5)*pow(reynolds,1.675);
    } else{
      dragCoeff = 2.61;
    }

  } else if(flowType==FLOW_PARTIC || flowType==FLOW_DROPLET){

    //Empirical drag coefficient correlations for a rigid sphere***//

    if(reynolds<0.5){
      dragCoeff = 24.0/reynolds;
    } else if(reynolds>=0.5 && reynolds<1000.0){
      dragCoeff = 24.0*(1.0+0.15*pow(reynolds,0.687))/reynolds;
    } else{
      dragCoeff = 0.44;
    }
  }

  return dragCoeff;
}


void plas::plas_CalcEntityCoefficients(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow)
{
  //***Compute flow coefficients***//

  ent->reynolds = plas_CalcDispReynolds(fp.nuCont,ent->diam,ent->normVel);
  ent->dragCoeff = plas_CalcDragCoeff(rp.flowType,ent->reynolds);
  ent->liftCoeff = plas_CalcLiftCoeff(rp.flowType);
  ent->kinRespTime = plas_CalcKinematicResponseTime(ent);
  ent->spalding = plas_CalcSpaldingNumber(flow->pressure);
  ent->prandtl = plas_CalcPrandtlNumber();
  ent->nusselt = plas_CalcNusseltNumber(ip.evapModel,ent->reynolds,ent->spalding,ent->prandtl);
  ent->thermRespTime = plas_CalcThermalResponseTime(ent->diam);
  ent->schmidt = plas_CalcSchmidtNumber();
  ent->sherwood = plas_CalcSherwoodNumber(ip.evapModel,ent->reynolds,ent->schmidt,ent->spalding);
  ent->massTrCoeff = plas_CalcMassTransferCoeff(ent->sherwood,ent->spalding);
  ent->pressBubble = plas_CalcPressBubble(ent->diam,flow->pressure);
  ent->rhoBubble = plas_CalcRhoBubble(ent->temp, ent->pressBubble);
  ent->concInterf = plas_CalcConcInterf( ent->pressBubble);
}


double plas::plas_CalcKinematicResponseTime(LOCAL_ENTITY_VARIABLES *ent)
{
  double tau = 0.;

  if(rp.flowType==FLOW_PARTIC || rp.flowType==FLOW_DROPLET){
    tau = 4.0*md.rhoDisp*ent->diam*ent->diam/(3.0*fp.muCont*ent->reynolds*ent->dragCoeff);
  } else if(rp.flowType==FLOW_BUBBLY){
    tau = 2.0*ent->diam*ent->diam/(3.0*fp.nuCont*ent->reynolds*ent->dragCoeff);
  }

  return tau;
}


double plas::plas_CalcLiftCoeff(int flowType)
{
  if (flowType==FLOW_PARTIC || flowType==FLOW_DROPLET)
    return 0.;
  else if (flowType==FLOW_BUBBLY)
    return 0.53;
  return 0.;
}


double plas::plas_CalcMassTransferCoeff(double sherwood, double spalding)
{
  double omega = log(1.0+spalding);

  return (2.0*(fp.rhoCont/md.rhoDisp)*md.binaryDiffCoeff*sherwood*omega);
}


void plas::plas_CalcMaterialData(double T, double p)
{
  if(ip.material==MAT_COPPER){

    //***Material data for copper***//

    md.rhoDisp = 8920.0;
    //md.muDisp = 0.0;
    md.cpDisp = 381.0;
    md.kDisp = 387.6;
    //md.sigDisp = 0.0;
    md.epsDisp = 0.03;
    //md.satPresDisp = 0.0;
    //md.vapPres = 0.0;
    //md.latHeatDisp = 0.0;
    //md.molarMassDisp = 0.0;
    //md.molarMassDispVap = 0.0;
    //md.binaryDiffCoeff = 0.0;
    //md.massDiffCoeff = 0.0;
    //md.HeDisp = 0.0;

  } else if(ip.material==MAT_POLY){

    //***Material data for C8H8***//

    md.rhoDisp = 1050.0;
    //md.muDisp = 0.0;
    md.cpDisp = 1300.0;
    md.kDisp = 0.08;
    //md.sigDisp = 0.0;
    //md.epsDisp = 0.0;
    //md.satPresDisp = 0.0;
    //md.vapPres = 0.0;
    //md.latHeatDisp = 0.0;
    //md.molarMassDisp = 0.0;
    //md.molarMassDispVap = 0.0;
    //md.binaryDiffCoeff = 0.0;
    //md.massDiffCoeff = 0.0;
    //md.HeDisp = 0.0;

  } else if(ip.material==MAT_WATER){

    //***Material data for water***//

    md.rhoDisp = 998.2;
    md.muDisp = 2.414e-5*pow(10.0,(247.8/(T-140.0)));
    md.cpDisp = (1.0/18.015)*(276376.0-2090.1*T+8.125*pow(T,2.0)-0.014116*pow(T,3.0)+9.3701e-6*pow(T,4.0));
    md.kDisp = (2.9041*(3.0+pow(1.0-(T/611.966),2.0/3.0)))/(3.0+20.0*pow(1.0-(273.0/611.966),2.0/3.0));
    md.sigDisp = 0.07275;
    md.epsDisp = 0.95;
    md.satPresDisp = exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T));
    //md.vapPres = 0.0;
    md.latHeatDisp = (5.2053e7/18.015)*pow((1.0-(T/611.966)),0.3199-0.212*(T/611.966)+0.25795*pow((T/611.966),2.0));;
    md.molarMassDisp = 18.015;
    md.molarMassDispVap = 29.0;
    md.binaryDiffCoeff = 0.220e-4;
    //md.massDiffCoeff = 0.0;
    //md.HeDisp = 0.0;

  } else if(ip.material==MAT_NHEPTANE){

    //***Material data for n-heptane***//

    md.rhoDisp = -941.03+19.9618*T-0.08612051*pow(T,2.0)+1.579494e-4*pow(T,3.0)-1.089345e-7*pow(T,4.0);
    md.muDisp = 0.386e-3;
    md.cpDisp = 799.3401062-126.5095282565*T+0.5279613848638*pow(T,2.0)-1664.890863e-8*pow(T,3.0)+644.6826474e-11*pow(T,4.0);
    md.kDisp = 0.25844890110-4.5549450549e-4*T;
    md.sigDisp = 0.059*pow(1.0-T/540.17,0.121);
    md.epsDisp = 1.0; // T and D dependancy still to be implemented...
    md.satPresDisp = 1.0e5*pow(10.0,4.02677-1258.34/(T-53.85));
    //md.vapPres = 0.0;
    md.latHeatDisp = 317.8e3*pow((540.17-T)/(540.17-371.4),0.38);
    md.molarMassDisp = 100.21;
    md.molarMassDispVap = 28.01;
    md.binaryDiffCoeff =
      ((1.8583e-7)/(p*pow(6.498,2.0)*(1.06036/(pow(T*399.3,0.15610)+0.19300/exp(0.47635*T*399.3))+1.76474/(exp(3.89411*T*399.3)))))
      *pow(pow(T,1.5)*(1.0/100.21+1.0/28.01),0.5);
    //md.massDiffCoeff = 0.0;
    //md.HeDisp = 0.0;

  } else if(ip.material==MAT_HYDROGEN){

    //***Material data for hydrogen***//

    md.rhoDisp = 0.0943;
    md.muDisp = 8.411e-6;
    md.cpDisp = 14283.0;
    md.kDisp = 0.1672;
    md.sigDisp = 0.07275;
    //md.epsDisp = 0.0;
    //md.satPresDisp = 0.0;
    md.vapPres = 2337.0;
    //md.latHeatDisp = 0.0;
    md.molarMassDisp = 2.0159;
    //md.molarMassDispVap = 0.0;
    //md.binaryDiffCoeff = 0.0;
    md.massDiffCoeff = 1.61e-9;
    md.HeDisp = 1.0/(7.099e4*101.325);

  } else if(ip.material==MAT_OXYGEN){

    //***Material data for oxygen***//

    md.rhoDisp = 1.2999;
    md.muDisp = 1.919e-5;
    md.cpDisp = 919.31;
    md.kDisp = 0.0246;
    //md.sigDisp = 0.0;
    //md.epsDisp = 0.0;
    //md.satPresDisp = 0.0;
    //md.vapPres = 0.0;
    //md.latHeatDisp = 0.0;
    md.molarMassDisp = 15.9994;
    //md.molarMassDispVap = 0.0;
    //md.binaryDiffCoeff = 0.0;
    //md.massDiffCoeff = 0.0;
    //md.HeDisp = 0.0;


  } else if(ip.material==MAT_AIR){

    //***Material data for air***//

    md.rhoDisp = (28.97e-3)/(Ru*T)*(p-exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T))+4.0*0.075/2e-5);//1.225;
    md.muDisp = 1.7894e-5;
    md.cpDisp = 1006.0;
    md.kDisp = 0.0242;
    md.sigDisp = 0.07275; //As bubbly flow air-water
    md.epsDisp = 0.0;
    md.satPresDisp = exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T));
    md.vapPres = 2337.0;
    md.latHeatDisp = 0.0;
    md.molarMassDisp = 28.97;
    md.molarMassDispVap = 0.0;
    md.binaryDiffCoeff = 0.0;
    md.massDiffCoeff = 2.0e-9;
    md.HeDisp = 2.1583e-7*exp(1342.0*(1.0/T-1.0/293.15));
  }
}


void plas::plas_CalcMatVectScalarProduct_2D(double *value, double **m, double *a)
{
  value[0] = m[0][0]*a[0]+m[0][1]*a[1];
  value[1] = m[1][1]*a[0]+m[1][1]*a[1];
}


void plas::plas_CalcMatVectScalarProduct_3D(double *value, double **m, double *a)
{
  value[0] = m[0][0]*a[0]+m[0][1]*a[1]+m[0][2]*a[2];
  value[1] = m[1][0]*a[0]+m[1][1]*a[1]+m[1][2]*a[2];
  value[2] = m[2][0]*a[0]+m[2][1]*a[1]+m[2][2]*a[2];
}


void plas::plas_CalcNodeImpactFactors(LOCAL_ENTITY_VARIABLES *ent, double *imp)
{
  double sum = 0.0;
  int idim,jdim,inod;
  double distance[fp.numDim];

  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    inod = ent->edata.elmNodes[idim];
    if(fp.flowSolver==FLOWSOLVER_SFELES){inod %= fp.numNod;}
    for(jdim=0; jdim<fp.numDim; jdim++){
      distance[jdim] = getNodCoord(inod,jdim)-ent->pos[jdim];
    }
    imp[idim] = plas_CalcVectorLength(fp.numDim,distance);
    if(imp[idim]<rp.errTol){imp[idim] = rp.errTol;}
    sum += 1.0/imp[idim];
  }

  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    imp[idim] = 1.0/(imp[idim]*sum);
  }
}


double plas::plas_CalcNusseltNumber(int evapModel, double reynolds, double spalding, double prandtl)
{
 double Nu = 2.0 + 0.6*sqrt(reynolds)*pow(prandtl,(1.0/3.0));

 if(evapModel){
   Nu *= (log(1.0+spalding)/spalding);
 }

 return Nu;
}


double plas::plas_CalcPrandtlNumber()
{
  return (fp.muCont*fp.cpCont/fp.kCont);
}


double plas::plas_CalcPressBubble(double diameter, double pressure)
{
  return (pressure-md.satPresDisp+4.0*md.sigDisp/diameter);
}


double plas::plas_CalcRhoBubble(double temperature, double pressBubble)
{
  return (md.molarMassDisp/(Ru*temperature)*pressBubble);
}


void plas::plas_CalcRotationMatrix_2D(double phi, double **m)
{
  m[0][0] = cos(phi);
  m[0][1] = -sin(phi);
  m[1][0] = sin(phi);
  m[1][1] = cos(phi);
}


void plas::plas_CalcRotationMatrix_3D(double phi, double **m, int axis)
{
  if(axis==0){
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = cos(phi);
    m[1][2] = -sin(phi);
    m[2][0] = 0.0;
    m[2][1] = sin(phi);
    m[2][2] = cos(phi);
  }else if(axis==1){
    m[0][0] = cos(phi);
    m[0][1] = 0.0;
    m[0][2] = -sin(phi);
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[2][0] = sin(phi);
    m[2][1] = 0.0;
    m[2][2] = cos(phi);
  }else if(axis==2){
    m[0][0] = cos(phi);
    m[0][1] = -sin(phi);
    m[0][2] = 0.0;
    m[1][0] = sin(phi);
    m[1][1] = cos(phi);
    m[1][2] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
  }
}


double plas::plas_CalcSchmidtNumber()
{
  return ((md.muDisp*fp.kCont)/(md.rhoDisp*fp.rhoCont*fp.cpCont));
}


double plas::plas_CalcSherwoodNumber(int evapModel, double reynolds, double schmidt, double spalding)
{
  double Sh = 2.0 + 0.6*pow(reynolds,(1.0/2.0))*pow(schmidt,(1.0/3.0));

  if(evapModel){
    double F_M = pow(1.0+spalding,0.7)*log(1.0+spalding)/spalding;
    Sh = 2.0+(Sh-2.0)/F_M;
  }

  return Sh;
}


double plas::plas_CalcSpaldingNumber(double pressure)
{
  if (md.satPresDisp<1.e-20 || md.molarMassDisp<1.e-20)
    return 0.;

  double Y_s = 1.0/(1.0+(pressure/md.satPresDisp-1.0)*(md.molarMassDispVap/md.molarMassDisp));
  return ((Y_s)/(1.0-Y_s));
}


double plas::plas_CalcThermalResponseTime(double diameter)
{
  return (1.0/(12.0*fp.kCont))*(md.rhoDisp*diameter*diameter*md.cpDisp);
}


void plas::plas_CalcTrajectory(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double dtLagr)
{
  int idim,jdim;
  double
    vmCoeff,
    pvmTerm,
    convTerm,
    radTerm,
    massTerm,
    concTerm,
    solvec[fp.numDim+2],
    rhsvec[fp.numDim+2],
    eps       = 1e-9,
    theta     = 0.5,
    boltzmann = 5.670400e-8;


  double **mat = new double*[fp.numDim+2];
  for (int r=0; r<fp.numDim+2; ++r)
    mat[r] = new double[fp.numDim+2];


  // TODO: Concentration definition depends on saturation model
  double alpha = 1.25;
  double flow_concentration = flow->pressure*md.HeDisp*alpha;

  //***Initialize data structures (u,v,w,T,d)***//

  for(idim=0; idim<fp.numDim; idim++){
    solvec[idim] = ent->vel[idim];
  }
  solvec[fp.numDim] = ent->temp;
  solvec[fp.numDim+1] = ent->diam;

  for(idim=0; idim<fp.numDim+2; idim++){
    for(jdim=0; jdim<fp.numDim+2; jdim++){
      mat[idim][jdim] = 0.0;
    }
    rhsvec[idim] = 0.0;
  }

  //***Add time contribution to matrix***//

  for(idim=0; idim<fp.numDim+2; idim++){
    mat[idim][idim] += 1.0/dtLagr;
  }

  //***Drag force contribution to matrix and RHS***//

  for(idim=0; idim<fp.numDim; idim++){
    mat[idim][idim] += theta/ent->kinRespTime;
    rhsvec[idim] += ent->relVel[idim]/ent->kinRespTime;
  }

  //***Lift force contribution to matrix and RHS***//

  if(ip.liftForce && rp.flowType==FLOW_BUBBLY && fp.numDim==3){
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

  if(rp.flowType==FLOW_BUBBLY){
    vmCoeff = 0.5*(1.0+2.786*pd[ent->node].volFrac);
    pvmTerm = 2.0*(1.0+0.5*(1.0+2.786*pd[ent->node].volFrac));
    for(idim=0; idim<fp.numDim; idim++){
      rhsvec[idim] += pvmTerm*flow->velDt[idim];
      for(jdim=0; jdim<fp.numDim; jdim++){
        rhsvec[idim] += pvmTerm*flow->vel[jdim]*flow->velDx[idim][jdim];
      }
    }
  }

  //***Gravitational force contribution to RHS***//

  for(idim=0; idim<fp.numDim; idim++){
    if(rp.flowType==FLOW_PARTIC || rp.flowType==FLOW_DROPLET){
      rhsvec[idim] += ip.gravVec[idim];
    } else if(rp.flowType==FLOW_BUBBLY){
      rhsvec[idim] -= 2.0*ip.gravVec[idim];
    }
  }

  //***Temperature equation contribution to matrix and RHS***//

  if(ip.energyCoupl && (rp.flowType==FLOW_PARTIC || rp.flowType==FLOW_DROPLET)){
    convTerm = ent->nusselt/(2.0*ent->thermRespTime);
    radTerm = (6.0*boltzmann*md.epsDisp)/(md.rhoDisp*md.cpDisp*ent->diam);
    if(ip.evapModel && rp.flowType==FLOW_DROPLET){
      massTerm = (3.0)*(ent->massTrCoeff/pow(ent->diam,2.0))*(md.latHeatDisp/md.cpDisp);
    } else{
      massTerm = 0.0;
    }
    mat[fp.numDim][fp.numDim] += theta*(convTerm-radTerm*4.0*pow(ent->temp,3.0));
    mat[fp.numDim][fp.numDim+1] -=
      theta*(-((1.0/(4.0*ent->thermRespTime*ent->diam))*(2.0+3.0*ent->nusselt)*(ent->relTemp))
   + ((radTerm/ent->diam)*(pow(ent->temp,4.0)-pow(flow->temp,4.0)))-((massTerm/ent->diam)*(1.5+1.0/ent->sherwood)));
    rhsvec[fp.numDim] += (convTerm*ent->relTemp)-(radTerm*(pow(ent->temp,4.0)-pow(flow->temp,4.0)))-massTerm;
  }

  //***Diameter equation contribution to matrix and RHS***//

  if(ip.evapModel && rp.flowType==FLOW_DROPLET ){
    mat[fp.numDim+1][fp.numDim+1] -= theta*(ent->massTrCoeff/pow(ent->diam,2.0))*(1.5-1.0/ent->sherwood);
    rhsvec[fp.numDim+1] += -ent->massTrCoeff/ent->diam;
  }

  //***Concentration Boundary Layer Model Payvar for bubble growth and Epstein & Plesset for Bubble shrink

  if(ip.saturModel && rp.flowType==FLOW_BUBBLY ){
   concTerm  = 4.0*md.massDiffCoeff*fp.rhoCont*(flow_concentration-ent->concInterf)/(ent->rhoBubble*(fp.rhoCont-ent->concInterf));

    if(flow_concentration > ent->concInterf){
      rhsvec[fp.numDim+1] += concTerm/ent->diam*(1.0+1.0/(pow(1.0+(flow_concentration-fp.rhoCont)/(ent->concInterf-flow_concentration)*(ent->rhoBubble/fp.rhoCont)*(1.0-(md.rhoDisp/ent->rhoBubble)*pow(2.0e-5/ent->diam,3.0)),0.5)-1.0));
      mat[fp.numDim+1][fp.numDim+1] -= concTerm*((1.0/(ent->diam+eps)+1.0/((ent->diam+eps)*pow(1.0+(flow_concentration-fp.rhoCont)/(ent->concInterf-flow_concentration)*(ent->rhoBubble/fp.rhoCont)*(1.0-(md.rhoDisp/ent->rhoBubble)*pow(2.0e-5/(ent->diam+eps),3.0)),0.5)-1.0))-(1.0/ent->diam+1.0/(ent->diam*pow(1.0+(flow_concentration-fp.rhoCont)/(ent->concInterf-flow_concentration)*(ent->rhoBubble/fp.rhoCont)*(1.0-(md.rhoDisp/ent->rhoBubble)*pow(2.0e-5/ent->diam,3.0)),0.5)-1.0)))/eps;

      }else{
      rhsvec[fp.numDim+1] += concTerm/ent->diam;
      mat[fp.numDim+1][fp.numDim+1] -= ((concTerm/(ent->diam+eps))-(concTerm/ent->diam))/eps;
        }
  }

  //***Solve linear system for velocity, temperature & diameter***//

  plas_SolveGaussSeidel(fp.numDim+2,mat,solvec,rhsvec);

  for(idim=0; idim<fp.numDim; idim++){
    ent->velOld[idim] = ent->vel[idim];
    ent->vel[idim] += solvec[idim];
  }
  ent->temp += solvec[fp.numDim];
  ent->diam += solvec[fp.numDim+1];

  //***Update entity position***//

  for(idim=0; idim<fp.numDim; idim++){
    ent->posOld[idim] = ent->pos[idim];
    ent->pos[idim] += 0.5*(ent->velOld[idim]+ent->vel[idim])*dtLagr;
  }


  for (int r=0; r<fp.numDim+2; ++r)
    delete[] mat[r];
  delete[] mat;
}


double plas::plas_CalcVectorAngle(int numDim, double *a, double *b)
{
  return acos(plas_CalcVectScalarProduct(numDim,a,b)/(plas_CalcVectorLength(numDim,a)*plas_CalcVectorLength(numDim,b)));
}


double plas::plas_CalcVectorLength(int numDim, double *a)
{
  return sqrt(plas_CalcVectScalarProduct(numDim,a,a));
}


void plas::plas_CalcVectorRotation_3D(double phi, double **m, double *a)
{
  m[0][0]=cos(phi)+(1-cos(phi))*a[0]*a[0];
  m[0][1]=(1-cos(phi))*a[0]*a[1]-sin(phi)*a[2];
  m[0][2]=(1-cos(phi))*a[0]*a[2]+sin(phi)*a[1];
  m[1][0]=(1-cos(phi))*a[0]*a[1]+sin(phi)*a[2];
  m[1][1]=cos(phi)+(1-cos(phi))*a[1]*a[1];
  m[1][2]=(1-cos(phi))*a[1]*a[2]-sin(phi)*a[0];
  m[2][0]=(1-cos(phi))*a[0]*a[2]-sin(phi)*a[1];
  m[2][1]=(1-cos(phi))*a[2]*a[1]+sin(phi)*a[0];
  m[2][2]=cos(phi)+(1-cos(phi))*a[2]*a[2];
}


double plas::plas_CalcVectScalarProduct(int numDim, double *a, double *b)
{
  if (numDim==2)
    return a[0]*b[0]+a[1]*b[1];
  else if (numDim==3)
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return 0.;
}


void plas::plas_CalcVorticity(int numDim, LOCAL_FLOW_VARIABLES *flow)
{
  if(numDim==2){
    flow->vort[0] = 0.0;
    flow->vort[1] = 0.0;
  } else if(numDim==3){
    flow->vort[0] = flow->velDx[2][1]-flow->velDx[1][2];
    flow->vort[1] = flow->velDx[0][2]-flow->velDx[2][0];
    flow->vort[2] = flow->velDx[1][0]-flow->velDx[0][1];
  }
}


double plas::plas_CalcWallFaceDistance(int numDim, double *pos, int ibnd, int ifac)
{
  int idim;
  double posVec[numDim],unitVec[numDim];

  //***Calculate unit normal of boundary segment***//

  plas_CalcBoundaryUnitNormal(numDim,ibnd,ifac,unitVec);

  //***Calculate and return wall face distance***//

  for(idim=0; idim<numDim; idim++){
    posVec[idim] = pos[idim]-getBndFaceRefCoord(ibnd,ifac,idim);
  }

  return plas_CalcVectScalarProduct(numDim,posVec,unitVec);
}


void plas::plas_CheckNaN(LOCAL_ENTITY_VARIABLES *ent)
{
  int idim;

  double checkNanPos = 0.0;
  double checkNanVel = 0.0;

  for(idim=0; idim<fp.numDim; idim++){
    checkNanPos += ent->pos[idim];
    checkNanVel += ent->vel[idim];
  }

  if(isnan(checkNanPos) || isnan(checkNanVel)){
    ent->flag = DFLAG_DISABLED;
    sd.lost++;
    screenWarning("Not-a-number detected");
  }
}


int plas::plas_Coalescence(double dj, double di, double *uijRelPrPr, double *x, double *y, double *z,
  double Mi, double Mj, double *uiPrPr, double *uiPrPrNew, double *uiNew, double *ui)
{
  double h0 = 1e-4;
  double hf = 1e-8;
  double Cc = 0.5;
  double sigma = 0.06;
  double Rij = 2.0*1/(2/dj+2/di);
  double T = sqrt(Rij*Rij*Rij*fp.rhoCont/(16*sigma))*log(h0/hf);
  double Tau_ij = Cc*Rij/uijRelPrPr[0];
  int idim,isCoal;
  double diStar;

  //***Coalescense model***//

  if(Tau_ij>T){

    isCoal = 1;
    sd.coalesc++;

    //***Compute new entity diameter***//

    diStar = pow((dj*dj*dj+di*di*di),1/3);

    if(fp.numDim==2){

      uiPrPrNew[0] = uiPrPr[0]+Mi/(Mi+Mj);
      uiPrPrNew[1] = uiPrPr[1];

      uiNew[0] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,x);
      uiNew[1] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,y);

    } else if(fp.numDim==3){

      uiPrPrNew[0] = uiPrPr[0]+Mi/(Mi+Mj);
      uiPrPrNew[1] = uiPrPr[1];
      uiPrPrNew[2] = uiPrPr[2];

      uiNew[0] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,x);
      uiNew[1] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,y);
      uiNew[2] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,z);
    }

    di = diStar;

    //***Update entity velocity***//

    for(idim=0; idim<fp.numDim; idim++){
      ui[idim] = uiNew[idim];
    }

  } else{
    isCoal = 0;
  }

  return isCoal;
}


void plas::plas_CollisionModel(LOCAL_ENTITY_VARIABLES *ent, int numDens, double dtLagr)
{
  int idim;
  double
    c_T = 0.3,
    e = 1.0,
    mu_static = 0.0,
    mu_dynamic = 0.0,
    ui[fp.numDim],
    uiFluct[fp.numDim],
    uj[fp.numDim],
    ujFluct[fp.numDim],
    uiNew[fp.numDim],
    uijRel[fp.numDim],
    Pcoll,
    FreqColl,
    conc,
    di,
    dj,
    tau_f,
    R,
    stokes,
    epsilon,
    rms,
    pos1[fp.numDim],
    pos1Pr[fp.numDim],
    pos2[fp.numDim],
    pos2Pr[fp.numDim],
    posRel[fp.numDim],
    x[fp.numDim],
    y[fp.numDim],
    z[fp.numDim],
    xPr[fp.numDim],
    yPr[fp.numDim],
    zPr[fp.numDim],
    xPrPr[fp.numDim],
    yPrPr[fp.numDim],
    zPrPr[fp.numDim],
    uiPrPr[fp.numDim],
    uiPrPrNew[fp.numDim],
    ujPrPr[fp.numDim],
    uijRelPrPr[fp.numDim],
    Phi,
    Psi,
    Mi,
    Mj,
    Jx,
    Jy,
    Jz,
    Y = 0.,
    Z = 0.,
    L = 0.,
    length;


  double **mat = new double*[fp.numDim];
  for (int r=0; r<fp.numDim; ++r)
    mat[r] = new double[fp.numDim];


  //***Set local entity properties (subscript i)***//

  for(idim=0; idim<fp.numDim; idim++){
    ui[idim] = ent->vel[idim];
  }
  di = ent->diam;

  //***Set fictitious collision partner properties (subscript j)***//

  for(idim=0; idim<fp.numDim; idim++){
    uj[idim] = plas_RandomGaussian(pd[ent->node].avgVel[idim],pd[ent->node].stdVel[idim]);
  }
  dj = plas_RandomGaussian(pd[ent->node].avgDiam,pd[ent->node].stdDiam);

  //***Correlation between the fluctuations of the real and fictitious particle***//

  if(ip.collisionModel==2){

    for(idim=0; idim<fp.numDim; idim++){
      uiFluct[idim] = ui[idim]-pd[ent->node].avgVel[idim];
    }

    tau_f = c_T*getEulerianTimeScale(ent->node);
    stokes = ent->kinRespTime/tau_f;
    epsilon = plas_RandomGaussian(0.0,1.0);
    R = exp(-0.55*pow(stokes,0.4));

    for(idim=0; idim<fp.numDim; idim++){
      rms = sqrt(pow(pd[ent->node].avgVel[idim],2.0)+pow(pd[ent->node].stdVel[idim],2.0));
      ujFluct[idim] = R*uiFluct[idim]+rms*sqrt(1-R*R)*epsilon;
      uj[idim] += ujFluct[idim];
    }
  }

  //***Calculate relative velocity***//

  for(idim=0; idim<fp.numDim; idim++){
    uijRel[idim] = ui[idim]-uj[idim];
  }

  //***Collision probability criterion***//

  conc = numDens/getNodVol(ent->node);
  FreqColl = PI/4.0*(di+dj)*(di+dj)*plas_CalcVectorLength(fp.numDim,uijRel)*conc;
  Pcoll = FreqColl*dtLagr;

  //***Perform collision***//

  if(plas_RandomDouble()<Pcoll){

    //***Transform into 2nd coordinate system***//

    x[0] = y[1] = z[2] = 1.0;
    x[1] = x[2] = y[0] = y[2] = z[0] = z[1] = 0.0;

    for(idim=0; idim<fp.numDim; idim++){
      xPr[idim] = uijRel[idim]/plas_CalcVectorLength(fp.numDim,uijRel);
    }

    plas_CalcCrossProduct_3D(yPr,xPr,x);
    plas_CalcCrossProduct_3D(zPr,xPr,yPr);

    pos1Pr[0] = plas_CalcVectScalarProduct(fp.numDim,ent->pos,xPr);
    pos1Pr[1] = plas_CalcVectScalarProduct(fp.numDim,ent->pos,yPr);
    pos1Pr[2] = plas_CalcVectScalarProduct(fp.numDim,ent->pos,zPr);

    while(L>1){
      Y = plas_RandomDouble();
      Z = plas_RandomDouble();
      L = sqrt(Y*Y+Z*Z);
    }

    Phi = asin(L);
    Psi = atan(Y/Z);

    plas_CalcRotationMatrix_3D(Phi,mat,2);
    plas_CalcMatVectScalarProduct_3D(pos2Pr,mat,xPr);
    plas_CalcVectorRotation_3D(Psi,mat,xPr);
    plas_CalcMatVectScalarProduct_3D(pos2Pr,mat,pos2Pr);

    for(idim=0; idim<fp.numDim; idim++){
      pos2Pr[idim] *= (dj+di)/2.0;
    }

    pos2[0] = plas_CalcVectScalarProduct(fp.numDim,pos2Pr,x);
    pos2[1] = plas_CalcVectScalarProduct(fp.numDim,pos2Pr,y);
    pos2[2] = plas_CalcVectScalarProduct(fp.numDim,pos2Pr,z);

    //***Transform into 3rd coordinate system***//

    for(idim=0; idim<fp.numDim; idim++){
      posRel[idim] = pos1[idim]-pos2[idim];
    }

    for(idim=0; idim<fp.numDim; idim++){
      xPrPr[idim] = (posRel[idim])/plas_CalcVectorLength(fp.numDim,posRel);
    }

    plas_CalcCrossProduct_3D(zPrPr,xPrPr,uijRel);

    for(idim=0; idim<fp.numDim; idim++){
      zPrPr[idim] /= plas_CalcVectorLength(fp.numDim,zPrPr);
    }

    plas_CalcCrossProduct_3D(yPrPr,zPrPr,xPrPr);

    uiPrPr[0] = plas_CalcVectScalarProduct(fp.numDim,ui,xPrPr);
    uiPrPr[1] = plas_CalcVectScalarProduct(fp.numDim,ui,yPrPr);
    uiPrPr[2] = plas_CalcVectScalarProduct(fp.numDim,ui,zPrPr);
    ujPrPr[0] = plas_CalcVectScalarProduct(fp.numDim,uj,xPrPr);
    ujPrPr[1] = plas_CalcVectScalarProduct(fp.numDim,uj,yPrPr);
    ujPrPr[2] = plas_CalcVectScalarProduct(fp.numDim,uj,zPrPr);
    uijRelPrPr[0] = plas_CalcVectScalarProduct(fp.numDim,uijRel,xPrPr);
    uijRelPrPr[1] = plas_CalcVectScalarProduct(fp.numDim,uijRel,yPrPr);
    uijRelPrPr[2] = plas_CalcVectScalarProduct(fp.numDim,uijRel,zPrPr);

    //***Sliding or non-sliding collision***//

    Mi = md.rhoDisp*PI*(dj/2.0)*(dj/2.0);
    Mj = md.rhoDisp*PI*(dj/2.0)*(dj/2.0);
    Jx = -(1.0-e)*uiPrPr[1]*(Mi*Mj)/(Mi+Mj);
    length = plas_CalcVectorLength(fp.numDim,uijRel);

    if(std::abs(length)<7.0/2.0*mu_static*(1.0+e)*std::abs(uiPrPr[0])){
      Jy = -2.0/7.0*uijRelPrPr[1]*Mi*Mj/(Mi+Mj);
      Jz = -2.0/7.0*uijRelPrPr[2]*Mi*Mj/(Mi+Mj);
    }else{
      Jy = -mu_dynamic*uijRelPrPr[1]/length*std::abs(Jx);
      Jz = -mu_dynamic*uijRelPrPr[2]/length*std::abs(Jx);
    }

    uiPrPrNew[0] = uiPrPr[0]+Jx/Mi;
    uiPrPrNew[1] = uiPrPr[1]+Jy/Mi;
    uiPrPrNew[2] = uiPrPr[2]+Jz/Mi;

    uiNew[0] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,x);
    uiNew[1] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,y);
    uiNew[2] = plas_CalcVectScalarProduct(fp.numDim,uiPrPrNew,z);

    sd.coll++;

    for(idim=0; idim<fp.numDim; idim++){
      ent->vel[idim] = uiNew[idim];
    }
  }

  for (int r=0; r<fp.numDim; ++r)
    delete[] mat[r];
  delete[] mat;
}


void plas::plas_CreateStatsFile(char *outpString)
{
  FILE *outpFile;
  if (!plas_MpiGetRank()) {
    outpFile = fopen(outpString,"w");
    if (outpFile==NULL)
      return;

    fprintf(outpFile,"Iter\tTime\tEnt\tIn\tOut\tBounc\tColl\tPer\tPass\tLost\tRe_disp\tNu_disp\tdt_Lagr\tsubit\n");
    fclose(outpFile);
  }
}


void plas::plas_CreateTecplotFile(char *outpString)
{
  if (!plas_MpiGetRank()) {

    FILE* outpFile;
    outpFile = fopen(outpString,"w");
    if (fp.numDim==2) {
      fprintf(outpFile,"VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"d\" \"T\" \"theta\"\n");
    } else if(fp.numDim==3) {
      fprintf(outpFile,"VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"d\" \"T\" \"theta\"\n");
    }
    fclose(outpFile);

  }
}


void plas::plas_DeallocateLocalEntityVar(LOCAL_ENTITY_VARIABLES *ent)
{
  delete[] ent->edata.elmNodes;
  for (int idim=0; idim<MAXELMNORMS; ++idim)
    delete[] ent->edata.elmNorms[idim];
  delete[] ent->edata.elmNorms;
  for (int idim=0; idim<MAXELMFACES; ++idim)
    delete[] ent->edata.elmFaceVectors[idim];
  delete[] ent->edata.elmFaceVectors;
  delete[] ent->vel;
  delete[] ent->velOld;
  delete[] ent->relVel;
  delete[] ent->pos;
  delete[] ent->posOld;
}


void plas::plas_DeallocateLocalFlowVar(int numDim, LOCAL_FLOW_VARIABLES *flow)
{
  delete[] flow->vel;
  delete[] flow->velDt;
  for (int idim=0; idim<numDim; ++idim)
    delete[] flow->velDx[idim];
  delete[] flow->velDx;
  delete[] flow->vort;
}


void plas::plas_FindExitFace(int numBnd, int numDim, LOCAL_ENTITY_VARIABLES *ent, int *f, int *i, int *j)
{
  int ibnd,ifac=0,found;
  double d;

  //***Loop over boundaries and faces to find the exit face***//

  found = 0;
  for(ibnd=0; ibnd<numBnd; ibnd++){
    for(ifac=0; ifac<getNumBndFaces(ibnd); ifac++){
      if(ent->elm==getBndDomElm(ibnd,ifac,ent->elm)){
        d = plas_CalcWallFaceDistance(numDim,ent->pos,ibnd,ifac);
        if(d<0){found = 1;}
      }
      if(found){break;}
    }
    if(found){break;}
  }

  //***Set the return parameters***//

  *f = found;
  *i = ibnd;
  *j = ifac;
}


void plas::plas_FindMinimumElementFaceDistance(int numDim, LOCAL_ENTITY_VARIABLES *ent, int *idx, double *dmin)
{
  int idim,jdim;
  double idist,posvec[numDim],normvec[numDim];

  for(idim=0; idim<ent->edata.numElmFaces; idim++){
    for(jdim=0; jdim<numDim; jdim++){
      posvec[jdim] = ent->pos[jdim]-ent->edata.elmFaceVectors[idim][jdim];
      normvec[jdim] = ent->edata.elmNorms[idim][jdim];
    }
    plas_NormalizeVector(numDim,normvec);
    idist = plas_CalcVectScalarProduct(numDim,posvec,normvec);
    if(idim==0){
      *idx = idim;
      *dmin = idist;
    } else{
      if(idist<*dmin){
        *idx = idim;
        *dmin = idist;
      }
    }
  }
}


int plas::plas_FindNearestElementNode(LOCAL_ENTITY_VARIABLES *ent)
{
  int idim,node=0;
  double impFac[ent->edata.numElmNodes],f_max=0.;

  plas_CalcNodeImpactFactors(ent,impFac);
  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    if(idim==0){
      node = ent->edata.elmNodes[idim];
      f_max = impFac[idim];
    } else{
      if(impFac[idim]>f_max){
        node = ent->edata.elmNodes[idim];
        f_max = impFac[idim];
      }
    }
  }

  if(fp.flowSolver==FLOWSOLVER_SFELES){node %= fp.numNod;}

  return node;
}


void plas::plas_ImposeExternal()
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  int numProc = plas_MpiGetNumProc();
  int *disps  = new int[numProc];
  int *recs   = new int[numProc];
  int ient,idim,totNewEnt;
  double *newPos,*newVel,*newDiam,*newTemp;

#ifdef MPI
  int iproc;
  int *dispsArr = new int[numProc];
  int *recsArr  = new int[numProc];
#endif

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);

  //***Broadcast entity data***//

  totNewEnt = plas_MpiAllSumInt(numExtEnt);
  newDiam = new double[totNewEnt];
  newTemp = new double[totNewEnt];
  newPos  = new double[totNewEnt*fp.numDim];
  newVel  = new double[totNewEnt*fp.numDim];

#ifdef MPI
  MPI_Allgather(&numExtEnt,1,MPI_INT,recs,1,MPI_INT,MPI_COMM_WORLD);

  disps[0] = 0;
  for(iproc=1; iproc<numProc; iproc++){
    disps[iproc] = disps[iproc-1] + recs[iproc-1];
  }

  for(iproc=0; iproc<numProc; iproc++){
    recsArr[iproc] = fp.numDim*recs[iproc];
    dispsArr[iproc] = fp.numDim*disps[iproc];
  }

  MPI_Allgatherv(extEntDiam,numExtEnt,MPI_DOUBLE,newDiam,recs,disps,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(extEntTemp,numExtEnt,MPI_DOUBLE,newTemp,recs,disps,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(extEntPos,numExtEnt*fp.numDim,MPI_DOUBLE,newPos,recsArr,dispsArr,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(extEntVel,numExtEnt*fp.numDim,MPI_DOUBLE,newVel,recsArr,dispsArr,MPI_DOUBLE,MPI_COMM_WORLD);
#else
  newDiam = extEntDiam;
  newTemp = extEntTemp;
  newPos = extEntPos;
  newVel = extEntVel;
#endif

  //***Loop over entities to generate***//

  for(ient=0; ient<totNewEnt; ient++){

    for(idim=0; idim<fp.numDim; idim++){
      ent.pos[idim] = newPos[fp.numDim*ient+idim];
    }
    for(idim=0; idim<fp.numDim; idim++){
      ent.vel[idim] = newVel[fp.numDim*ient+idim];
    }
    ent.diam = newDiam[ient];
    ent.temp = newTemp[ient];

    //***Element search***//

    plas_SearchDomainParallel(&ent);

    //***Initialize entity***//

    if(ent.flag==DFLAG_ENABLED && sd.enabled<ip.numMaxEnt){

      plas_SetElementGeometry(fp.numDim,&ent);
      plas_Interpolate(&ent,&flow,0.0);

      ed[sd.enabled].flag = DFLAG_CREATED;
      ed[sd.enabled].element = ent.elm;
      ed[sd.enabled].node = plas_FindNearestElementNode(&ent);
      ed[sd.enabled].diameter = ent.diam;
      ed[sd.enabled].temperature = ent.temp;
      for(idim=0; idim<fp.numDim; idim++){
        ed[sd.enabled].position[idim] = ent.pos[idim];
        ed[sd.enabled].velocity[idim] = flow.vel[idim]+ent.vel[idim];
      }
      sd.enabled++;
      sd.in++;
    }
  }

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(fp.numDim,&flow);

  delete[] newDiam;
  delete[] newTemp;
  delete[] newPos;
  delete[] newVel;
  delete[] disps;
  delete[] recs;
}


void plas::plas_ImposeProductionDomains()
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  int ient,ipd,idim,iterate,bCtr;
  double p,s,mass=0.,p1[3],p2[3];
  double *newPos,*newDiam;
  char msg[100];

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);

  newDiam = new double[ip.numMaxEnt];
  newPos  = new double[ip.numMaxEnt*fp.numDim];

  //***Creation of bubbles is performed only by the master process***//

  if(plas_MpiGetRank()==0){

    bCtr = 0;

    //***Loop over production domains***//

    for(ipd=0; ipd<ip.numProdDom; ipd++){

      if(ip.prodDom[ipd]==0){continue;}

      //***Accumulate produced secondary phase mass***//

      rp.massResid[ipd] += ip.massFluxes[ipd]*fp.dtEul;
      if(rp.massResid[ipd]<=0.0){continue;}

      //***Copy production parameters to local data structure***//

      for(idim=0; idim<3; idim++){
        p1[idim] = ip.prodParam[ipd][idim];
        p2[idim] = ip.prodParam[ipd][idim+3];
      }

      //***Iteratively generate entities according to mass flux***//

      iterate = 1;
      do{

        //***Check if maximum number of entities is not exceedes***//

        if(bCtr==ip.numMaxEnt){break;}

        //***Set diameter***//

        newDiam[bCtr] = plas_SetDiameter();

        //***Calculate mass of generated entity***//

        if(fp.numDim==2){
          mass = md.rhoDisp*PI*newDiam[bCtr]*newDiam[bCtr]/4.0;
        } else if(fp.numDim==3){
          mass = md.rhoDisp*PI*newDiam[bCtr]*newDiam[bCtr]*newDiam[bCtr]/6.0;
        }

        rp.massResid[ipd] -= mass;
        if(rp.massResid[ipd]<=0.0){iterate = 0;}

        //***Compute position***//

        if(ip.prodDom[ipd]==1){

          //***Line production domain***//

          s = plas_RandomDouble();
          for(idim=0; idim<fp.numDim; idim++){
            newPos[fp.numDim*bCtr+idim] = p1[idim]+s*(p2[idim]-p1[idim]);
          }

        } else if(ip.prodDom[ipd]==2){

          //***Rectangle production domain***//

          for(idim=0; idim<fp.numDim; idim++){
            s = plas_RandomDouble();
            newPos[fp.numDim*bCtr+idim] = p1[idim]+s*(p2[idim]-p1[idim]);
          }

        } else if(ip.prodDom[ipd]==3){

          //***Ellipse production domain***//

          do{
            p = 0.0;
            for(idim=0; idim<fp.numDim; idim++){
              s = 2.0*plas_RandomDouble()-1.0;
              newPos[fp.numDim*bCtr+idim] = p1[idim]+s*p2[idim];
              if(p2[idim]>rp.errTol){
                p += (newPos[fp.numDim*bCtr+idim]-p1[idim])*(newPos[fp.numDim*bCtr+idim]-p1[idim])/(p2[idim]*p2[idim]);
              }
            }
          } while(p>1.0);
        }

        bCtr++;

      } while(iterate);
    }
  }

  //***Broadcast entity data***//

  plas_MpiBroadcastInt(&bCtr,1,0);
  plas_MpiBroadcastDouble(newDiam,ip.numMaxEnt,0);
  plas_MpiBroadcastDouble(newPos,fp.numDim*ip.numMaxEnt,0);
  plas_MpiBroadcastDouble(rp.massResid,ip.numProdDom,0);

  //***Generate entities//

  for(ient=0; ient<bCtr; ient++){

    ent.diam = newDiam[ient];
    ent.temp = ip.iniTempDisp;
    for(idim=0; idim<fp.numDim; idim++){
      ent.pos[idim] = newPos[fp.numDim*ient+idim];
      ent.vel[idim] = ip.iniVel[idim];
    }

    //***Element search***//

    plas_SearchDomainParallel(&ent);

    //***Initialize entity***//

    if(ent.flag==DFLAG_ENABLED && sd.enabled<ip.numMaxEnt){

      plas_SetElementGeometry(fp.numDim,&ent);
      plas_Interpolate(&ent,&flow,0.0);

      ed[sd.enabled].flag = DFLAG_CREATED;
      ed[sd.enabled].element = ent.elm;
      ed[sd.enabled].node = plas_FindNearestElementNode(&ent);
      ed[sd.enabled].diameter = ent.diam;
      ed[sd.enabled].temperature = ent.temp;
      for(idim=0; idim<fp.numDim; idim++){
        ed[sd.enabled].position[idim] = ent.pos[idim];
        ed[sd.enabled].velocity[idim] = flow.vel[idim]+ent.vel[idim];
      }
      sd.enabled++;
      sd.in++;
    }
  }
  sprintf(msg,"imposed entities: %d\n",bCtr);
  screenOutput(msg);

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(fp.numDim,&flow);

  delete[] newDiam;
  delete[] newPos;
}


void plas::plas_Interpolate(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim;

  //***Interpolate velocity, temperature and pressure***//

  plas_InterpolateVelocity(ent,flow,step);
  plas_CalcVorticity(fp.numDim,flow);
  plas_InterpolateTemperature(ent,flow,step);
  plas_InterpolatePressure(ent,flow,step);

  //***Compute relative and normal velocity***//

  ent->normVel = 0.0;
  for(idim=0; idim<fp.numDim; idim++){
    ent->relVel[idim] = flow->vel[idim]-ent->vel[idim];
    ent->normVel += ent->relVel[idim]*ent->relVel[idim];
  }
  ent->normVel = sqrt(ent->normVel);

  //***Compute relative temperature***//

  ent->relTemp = flow->temp-ent->temp;
}


void plas::plas_InterpolatePressure(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim;
  double impFac[ent->edata.numElmNodes];

  //***Compute impact factors of element nodes***//

  plas_CalcNodeImpactFactors(ent,impFac);

  //***Flow pressure***//

  flow->pressure = 0.0;
  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    flow->pressure += impFac[idim]*((1.0-step)*getPressureOld(ent->edata.elmNodes[idim])
                                        + step*getPressure(ent->edata.elmNodes[idim]));
  }
}


void plas::plas_InterpolateTemperature(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim;
  double impFac[ent->edata.numElmNodes];

  //***Compute impact factors of element nodes***//

  plas_CalcNodeImpactFactors(ent,impFac);

  //***Flow temperature***//

  flow->temp = 0.0;
  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    flow->temp += impFac[idim]*((1.0-step)*getTemperatureOld(ent->edata.elmNodes[idim])
                                    + step*getTemperature(ent->edata.elmNodes[idim]));
  }
}


void plas::plas_InterpolateVelocity(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim,jdim,kdim;
  double impFac[ent->edata.numElmNodes];

  //***Compute impact factors of element nodes***//

  plas_CalcNodeImpactFactors(ent,impFac);

  //***Flow velocity***//

  for(idim=0; idim<fp.numDim; idim++){
    flow->vel[idim] = 0.0;
    for(jdim=0; jdim<ent->edata.numElmNodes; jdim++){
      flow->vel[idim] += impFac[jdim]*((1.0-step)*getVelocityCompOld(ent->edata.elmNodes[jdim],idim)
                                           + step*getVelocityComp(ent->edata.elmNodes[jdim],idim));
    }
  }

  //***Flow velocity space derivatives***//

  for(idim=0; idim<fp.numDim; idim++){
    for(jdim=0; jdim<fp.numDim; jdim++){
      flow->velDx[idim][jdim] = 0.0;
      for(kdim=0; kdim<ent->edata.numElmNodes; kdim++){
        flow->velDx[idim][jdim] += impFac[kdim]*((1.0-step)*getVelocityDerivativeCompOld(ent->edata.elmNodes[kdim],idim,jdim)
                                                     + step*getVelocityDerivativeComp(ent->edata.elmNodes[kdim],idim,jdim));
      }
    }
  }

  //***Flow velocity time derivative***//

  for(idim=0; idim<fp.numDim; idim++){
    flow->velDt[idim] = 0.0;
    if (fp.dtEul>1.e-20) {
      for(jdim=0; jdim<ent->edata.numElmNodes; jdim++){
        flow->velDt[idim] += impFac[jdim]*(getVelocityComp(ent->edata.elmNodes[jdim],idim)
                                         - getVelocityCompOld(ent->edata.elmNodes[jdim],idim))/fp.dtEul;
      }
    }
  }
}


void plas::plas_LoadInitialDistribution(char *inpString)
{
  LOCAL_ENTITY_VARIABLES ent;
  int ient,jent,idim,numIni;
  char buffer[200],*search;
  long fpos;
  FILE *inpFile;

  //***Read data file***//

  inpFile = fopen(inpString,"r");

  /* locate zone starting position ("^ZONE.+(dispersed)")) */
  /* set numIni (I) and fp.time (SOLUTIONTIME) */
  /* back to correct position */
  fpos = 0L;
  fp.time = 0.;
  while (fgets(buffer,200,inpFile)!=NULL) {
    if (strstr(buffer,"ZONE")==buffer && strstr(buffer,"(dispersed)")!=NULL) {
      if ((search = strstr(buffer," I="))!=NULL)
        sscanf(search+3,"%d",&numIni);
      if ((search = strstr(buffer," SOLUTIONTIME="))!=NULL)
        sscanf(search+14,"%lf",&fp.time);
      fpos = ftell(inpFile);
    }
  }
  numIni = (numIni<ip.numMaxEnt? numIni:ip.numMaxEnt);
  fseek(inpFile,fpos,SEEK_SET);

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);

  //***Loop over entities to generate***//

  jent = 0;
  for (ient=0; ient<numIni && fgets(buffer,200,inpFile)!=NULL; ++ient) {

    //***Read entity data from data file***//

    if (fp.numDim==2) {
      sscanf(buffer,"%lf %lf %lf %lf %lf %lf",
        &ent.pos[0],&ent.pos[1],
        &ent.vel[0],&ent.vel[1],
        &ent.temp,&ent.diam);
    }
    else if (fp.numDim==3) {
      sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",
        &ent.pos[0],&ent.pos[1],&ent.pos[2],
        &ent.vel[0],&ent.vel[1],&ent.vel[2],
        &ent.temp,&ent.diam);
    }

    //***Perform element search***//

    plas_SearchDomainParallel(&ent);

    //***Initialize entity***//

    if (ent.flag==DFLAG_ENABLED) {
      ed[jent].flag = DFLAG_ENABLED;
      ed[jent].element = ent.elm;
      ed[jent].node = plas_FindNearestElementNode(&ent);
      ed[jent].diameter = ent.diam;
      ed[jent].temperature = ent.temp;
      for(idim=0; idim<fp.numDim; idim++){
        ed[jent].position[idim] = ent.pos[idim];
        ed[jent].velocity[idim] = ent.vel[idim];
      }
      jent++;
    }
  }

  fclose(inpFile);

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
}


void plas::plas_MpiAllMaxIntArray(int *val, int size)
{
#ifdef MPI
  int i,sendval[size],recval[size];
  for(i=0; i<size; i++){sendval[i] = val[i];}
  MPI_Allreduce(sendval,recval,size,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  for(i=0; i<size; i++){val[i] = recval[i];}
#else
  return;
#endif
}


int plas::plas_MpiAllMaxInt(int val)
{
#ifdef MPI
  int recval;
  MPI_Allreduce(&val,&recval,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  return recval;
#else
  return val;
#endif
}


int plas::plas_MpiAllMinInt(int val)
{
#ifdef MPI
  int recval;
  MPI_Allreduce(&val,&recval,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  return recval;
#else
  return val;
#endif
}


void plas::plas_MpiAllSumDoubleArray(double *val, int size)
{
#ifdef MPI
  int i;
  double sendval[size],recval[size];
  for(i=0; i<size; i++){sendval[i] = val[i];}
  MPI_Allreduce(sendval,recval,size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=0; i<size; i++){val[i] = recval[i];}
#else
  return;
#endif
}


double plas::plas_MpiAllSumDouble(double val)
{
#ifdef MPI
  double recval;
  MPI_Allreduce(&val,&recval,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return recval;
#else
  return val;
#endif
}


void plas::plas_MpiAllSumIntArray(int *val, int size)
{
#ifdef MPI
  int i,sendval[size],recval[size];
  for(i=0; i<size; i++){sendval[i] = val[i];}
  MPI_Allreduce(sendval,recval,size,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  for(i=0; i<size; i++){val[i] = recval[i];}
#else
  return;
#endif
}


int plas::plas_MpiAllSumInt(int val)
{
#ifdef MPI
  int recval;
  MPI_Allreduce(&val,&recval,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return recval;
#else
  return val;
#endif
}


void plas::plas_MpiBarrier()
{
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#else
  return;
#endif
}


void plas::plas_MpiBroadcastDouble(double *variable, int size, int root)
{
#ifdef MPI
  MPI_Bcast(variable,size,MPI_DOUBLE,root,MPI_COMM_WORLD);
#else
  return;
#endif
}


void plas::plas_MpiBroadcastInt(int *variable, int size, int root)
{
#ifdef MPI
  MPI_Bcast(variable,size,MPI_INT,root,MPI_COMM_WORLD);
#else
  return;
#endif
}


int plas::plas_MpiGetNumProc()
{
#ifdef MPI
  int numProc;
  MPI_Comm_size(MPI_COMM_WORLD,&numProc);
  return numProc;
#else
  return 1;
#endif
}


int plas::plas_MpiGetRank()
{
#ifdef MPI
  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  return irank;
#else
  return 0;
#endif
}


void plas::plas_NormalizeVector(int numDim, double *a)
{
  int idim;
  double length = plas_CalcVectorLength(numDim,a);

  for(idim=0; idim<numDim; idim++){
    a[idim] /= length;
  }
}


void plas::plas_PassEntities()
{
  LOCAL_ENTITY_VARIABLES ent;
  int numProc        = plas_MpiGetNumProc();
  int irank          = plas_MpiGetRank();
  int *disps         = new int[numProc];
  int *recs          = new int[numProc];
  int leftCtrOneProc = sd.leftproc;
  int leftCtrAllProc = plas_MpiAllSumInt(sd.leftproc);
  int *foundProc     = new int[leftCtrAllProc];
  int *foundElm      = new int[leftCtrAllProc];
  int *foundNode     = new int[leftCtrAllProc];
  int idim,jdim,ient,jent,ielm,inod;
  char errMessage[100];

  double *leftDataOneProc = new double[(2*fp.numDim+1)*leftCtrOneProc];
#ifdef MPI
  double *leftDataAllProc = new double[(2*fp.numDim+1)*leftCtrAllProc];
#else
  double *leftDataAllProc = leftDataOneProc;
#endif

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);

  //***Collect all entities that left a process***//

  jdim = 0;
  for(ient=0; ient<ip.numMaxEnt; ient++){
    if(ed[ient].flag==DFLAG_PASS){
      ed[ient].flag = DFLAG_DISABLED;
      for(idim=0; idim<fp.numDim; idim++){
        leftDataOneProc[jdim] = ed[ient].position[idim];
        jdim++;
      }
      for(idim=0; idim<fp.numDim; idim++){
        leftDataOneProc[jdim] = ed[ient].velocity[idim];
        jdim++;
      }
      leftDataOneProc[jdim] = ed[ient].diameter;
      jdim++;
    }
  }

  //***Gather all entities that left any process***//

#ifdef MPI
  MPI_Allgather(&leftCtrOneProc,1,MPI_INT,recs,1,MPI_INT,MPI_COMM_WORLD);

  disps[0] = 0;
  for (int iproc=1; iproc<numProc; iproc++){
    disps[iproc] = disps[iproc-1] + recs[iproc-1];
  }

  for (int iproc=0; iproc<numProc; iproc++){
    recs[iproc] *= (2*fp.numDim+1);
    disps[iproc] *= (2*fp.numDim+1);
  }

  MPI_Allgatherv(leftDataOneProc,(2*fp.numDim+1)*leftCtrOneProc,MPI_DOUBLE,leftDataAllProc,recs,disps,MPI_DOUBLE,MPI_COMM_WORLD);
#endif

  //***Search for all free entities on all processes***//

  idim = 0;
  for(ient=0; ient<leftCtrAllProc; ient++){

    for(jdim=0; jdim<fp.numDim; jdim++){
      ent.pos[jdim] = leftDataAllProc[idim+jdim];
    }

    //***Element search***//

    plas_SearchDomainParallel(&ent);

    if(ent.flag==DFLAG_ENABLED){
      foundProc[ient] = irank;
      foundElm[ient] = ent.elm;
      foundNode[ient] = ent.node;
    } else{
      foundProc[ient] = -1;
    }
    idim += 2*fp.numDim+1;
  }

  //***Determine a definite process on which an entity is found***//

  plas_MpiAllMaxIntArray(foundProc,leftCtrAllProc);

  //***Broadcast information of found entities***//

  idim = 0;
  jent = 0;

  for(ient=0; ient<leftCtrAllProc; ient++){

    if(foundProc[ient]>-1){
      ielm = foundElm[ient];
      plas_MpiBroadcastInt(&ielm,1,foundProc[ient]);
      inod = foundNode[ient];
      plas_MpiBroadcastInt(&inod,1,foundProc[ient]);
    }

    if(foundProc[ient]!=irank){
      idim += 2*fp.numDim+1;
      continue;
    }

    while(ed[jent].flag==DFLAG_ENABLED){
      jent++;
      if(jent==ip.numMaxEnt){
        sprintf(errMessage,"Memory exceeded when passing trajectory information.");
        plas_TerminateOnError(errMessage);
      }
    }

    //***Generate entity on the new process***//

    sd.passed++;
    ed[jent].flag = DFLAG_ENABLED;
    ed[jent].element = ielm;
    ed[jent].node = inod;

    for(jdim=0; jdim<fp.numDim; jdim++){
      ed[jent].position[jdim] = leftDataAllProc[idim];
      idim++;
    }

    for(jdim=0; jdim<fp.numDim; jdim++){
      ed[jent].velocity[jdim] = leftDataAllProc[idim];
      idim++;
    }

    ed[jent].diameter = leftDataAllProc[idim];
    idim++;
  }

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);

  delete[] foundProc;
  delete[] foundElm;
  delete[] foundNode;
  delete[] leftDataOneProc;
#ifdef MPI
  delete[] leftDataAllProc;
#endif
  delete[] disps;
  delete[] recs;
}


double plas::plas_RandomDouble()
{
  int r = plas_RandomInteger(0,(int)1e5);

  return (double)(r/1e5);
}


double plas::plas_RandomGaussian(float m, float s)
{
  double xx1, xx2, w, yy1;
  static double yy2;
  static int use_last = 0;

  if (use_last){
    yy1 = yy2;
    use_last = 0;
  } else{
    do {
      xx1 = 2.0 * plas_RandomDouble() - 1.0;
      xx2 = 2.0 * plas_RandomDouble() - 1.0;
      w = xx1 * xx1 + xx2 * xx2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    yy1 = xx1 * w;
    yy2 = xx2 * w;
    use_last = 1;
  }

  return( m + yy1 * s );
}


void plas::plas_RandomInitialDistribution()
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  int ient,idim,numIni;
  double p,q,rand1,rand2,rand3,rand4,setDiam,partVolume,totalVolume,elmVolume;

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);

  //***Initializations***//

  partVolume = fp.domainVolume;
  totalVolume = plas_MpiAllSumDouble(partVolume);
  numIni = (int)(ip.numIniEnt*partVolume/totalVolume);
  if(numIni>ip.numMaxEnt){numIni = ip.numMaxEnt;}

  //***Loop over entities to generate***//

  for(ient=0; ient<numIni; ient++){

    //***Generate entity in a random element***//

    do{
      ent.elm = plas_RandomInteger(0,fp.numElm-1);
      plas_SetElementGeometry(fp.numDim,&ent);
      elmVolume = getElmVol(ent.elm);
      p = plas_RandomDouble();
      q = elmVolume/fp.maxElmVolume;
    } while(p>q);

    //***Set diameter***//

    setDiam = plas_SetDiameter();

    //***Random position according to element type***//

    if(getElementType(ent.elm)==ELM_SIMPLEX){

      rand1 = plas_RandomDouble();
      rand2 = (1.0-rand1)*plas_RandomDouble();
      rand3 = (1.0-rand1-rand2)*plas_RandomDouble();
      for(idim=0; idim<fp.numDim; idim++){
        ent.pos[idim] = getNodCoord(ent.edata.elmNodes[0],idim)
          + rand1*(getNodCoord(ent.edata.elmNodes[1],idim)-getNodCoord(ent.edata.elmNodes[0],idim))
          + rand2*(getNodCoord(ent.edata.elmNodes[2],idim)-getNodCoord(ent.edata.elmNodes[0],idim));
        if(fp.numDim==3){
          ent.pos[idim] +=
            rand3*(getNodCoord(ent.edata.elmNodes[3],idim)-getNodCoord(ent.edata.elmNodes[0],idim));
        }
      }

    } else if(getElementType(ent.elm)==ELM_PRISM){

      rand1 = plas_RandomDouble();
      rand2 = (1.0-rand1)*plas_RandomDouble();
      rand3 = plas_RandomDouble();

      for(idim=0; idim<fp.numDim; idim++){
        ent.pos[idim] = rand3*getNodCoord(ent.edata.elmNodes[0],idim)
          + rand3*rand1*(getNodCoord(ent.edata.elmNodes[1],idim)-getNodCoord(ent.edata.elmNodes[0],idim))
          + rand3*rand2*(getNodCoord(ent.edata.elmNodes[2],idim)-getNodCoord(ent.edata.elmNodes[0],idim))
          + (1.0-rand3)*getNodCoord(ent.edata.elmNodes[3],idim)
          + (1.0-rand3)*rand1*(getNodCoord(ent.edata.elmNodes[4],idim)-getNodCoord(ent.edata.elmNodes[3],idim))
          + (1.0-rand3)*rand2*(getNodCoord(ent.edata.elmNodes[5],idim)-getNodCoord(ent.edata.elmNodes[3],idim));
      }

    } else if(getElementType(ent.elm)==ELM_QUAD){

      rand1 = plas_RandomDouble();
      rand2 = plas_RandomDouble();
      rand3 = plas_RandomDouble();

      for(idim=0; idim<fp.numDim; idim++){
        ent.pos[idim] = rand3*(getNodCoord(ent.edata.elmNodes[0],idim)
                               + rand1*(getNodCoord(ent.edata.elmNodes[1],idim)
                                        - getNodCoord(ent.edata.elmNodes[0],idim)))
          + (1.0-rand3)*(getNodCoord(ent.edata.elmNodes[2],idim)
                         + rand2*(getNodCoord(ent.edata.elmNodes[3],idim)
                                  - getNodCoord(ent.edata.elmNodes[2],idim)));
      }

    } else if(getElementType(ent.elm)==ELM_HEX){

      rand1 = plas_RandomDouble();
      rand2 = plas_RandomDouble();
      rand3 = plas_RandomDouble();
      rand4 = plas_RandomDouble();

      for(idim=0; idim<fp.numDim; idim++){
        ent.pos[idim] = rand4*(rand3*(getNodCoord(ent.edata.elmNodes[0],idim)
                                      + rand1*(getNodCoord(ent.edata.elmNodes[4],idim)
                                               - getNodCoord(ent.edata.elmNodes[0],idim)))
                               + (1.0-rand3)*(getNodCoord(ent.edata.elmNodes[1],idim)
                                              + rand2*(getNodCoord(ent.edata.elmNodes[5],idim)
                                                       - getNodCoord(ent.edata.elmNodes[1],idim))))
          +(1.0-rand4)*(rand3*(getNodCoord(ent.edata.elmNodes[2],idim)
                              + rand1*(getNodCoord(ent.edata.elmNodes[6],idim)
                                       - getNodCoord(ent.edata.elmNodes[2],idim)))
                       + (1.0-rand3)*(getNodCoord(ent.edata.elmNodes[3],idim)
                         + rand2*(getNodCoord(ent.edata.elmNodes[7],idim)
                                  - getNodCoord(ent.edata.elmNodes[3],idim))));
      }

    } else if(getElementType(ent.elm)==ELM_PYRAMID){

      rand1 = plas_RandomDouble();
      rand2 = plas_RandomDouble();
      rand3 = plas_RandomDouble();
      rand4 = plas_RandomDouble();

      for(idim=0; idim<fp.numDim; idim++){
        ent.pos[idim] = rand4*getNodCoord(ent.edata.elmNodes[4],idim)
          + (1.0-rand4)*(rand3*(getNodCoord(ent.edata.elmNodes[0],idim)
                                + rand1*(getNodCoord(ent.edata.elmNodes[1],idim)
                                         - getNodCoord(ent.edata.elmNodes[0],idim)))
                         + (1.0-rand3)*(getNodCoord(ent.edata.elmNodes[2],idim)
                                        + rand2*(getNodCoord(ent.edata.elmNodes[3],idim)
                                                 - getNodCoord(ent.edata.elmNodes[2],idim))));
      }
    }

    //***Initialize entity***//

    plas_Interpolate(&ent,&flow,0.0);

    ed[ient].flag = DFLAG_CREATED;
    ed[ient].element = ent.elm;
    ed[ient].node = plas_FindNearestElementNode(&ent);
    ed[ient].diameter = setDiam;
    ed[ient].temperature = ip.iniTempDisp;
    for(idim=0; idim<fp.numDim; idim++){
      ed[ient].position[idim] = ent.pos[idim];
      ed[ient].velocity[idim] = flow.vel[idim]+ip.iniVel[idim];
    }
    sd.in++;
  }

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(fp.numDim,&flow);
}


int plas::plas_RandomInteger(int min, int max)
{
  return min + (int)( ((double)(max - min) + 1.) * rand()/(RAND_MAX + 1.));
}


double plas::plas_ReadDoubleParam(FILE *inpFile)
{
  double d;
  int ignore_i;
  char text[100],*ignore_cp;

  ignore_cp = fgets(text,100,inpFile);
  ignore_i = fscanf(inpFile,"%lf",&d);
  ignore_cp = fgets(text,100,inpFile);

  return d;
}


int plas::plas_ReadIntParam(FILE *inpFile)
{
  int i,ignore_i;
  char text[100],*ignore_cp;

  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&i);
  ignore_cp = fgets(text,100,inpFile);

  return i;
}


void plas::plas_ReadParameters()
{
  int i,j,ignore_i;
  char inpText[200],fileString[200],errMessage[200],*ignore_cp;
  FILE *inpFile;

  //***Open input file***//

  inpFile = NULL;
  if (ip.confFilename!=NULL) {
    if (strlen(ip.confFilename)) {
      inpFile = fopen(ip.confFilename,"r");
      if (inpFile==NULL) {
        sprintf(errMessage,"Configuration file: \"%s\" was not found.\n",ip.confFilename);
        screenOutput(errMessage);
      }
    }
  }
  if (inpFile==NULL) {
    sprintf(fileString,"./plas.conf");
    inpFile = fopen(fileString,"r");
    if (inpFile==NULL) {
      sprintf(errMessage,"Configuration file: \"%s\" was not found.",fileString);
      plas_TerminateOnError(errMessage);
    }
  }

  //***Read maximum number of entities***//

  ip.numMaxEnt = plas_ReadIntParam(inpFile);
  if(ip.numMaxEnt<0){
    sprintf(errMessage,"Bad value for maximum number of dispersed entities.");
    plas_TerminateOnError(errMessage);
  }

  //***Read number of initially distributed entities***//

  ip.numIniEnt = plas_ReadIntParam(inpFile);
  if(ip.numIniEnt<0){
    sprintf(errMessage,"Bad value for number of initially distributed dispersed entities.");
    plas_TerminateOnError(errMessage);
  }

  //***Read production domains and mass fluxes***//

  ip.numProdDom = plas_ReadIntParam(inpFile);
  if(ip.numProdDom<0){
    sprintf(errMessage,"Invalid number of production domains.");
    plas_TerminateOnError(errMessage);
  }

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);

  if(ip.numProdDom>0){

    ip.prodDom = new int[ip.numProdDom];
    if (ip.numProdDom>0) {
      ip.prodParam  = new double*[ip.numProdDom];
      ip.massFluxes = new double [ip.numProdDom];
    }

    for(i=0; i<ip.numProdDom; i++){
      ignore_i = fscanf(inpFile,"%d",&ip.prodDom[i]);
      if(ip.prodDom[i]<1 || ip.prodDom[i]>3){
        sprintf(errMessage,"Bad value for production domain type #%d.",i);
        plas_TerminateOnError(errMessage);
      }
      ip.prodParam[i] = new double[6];
      for(j=0; j<6; j++){
        ignore_i = fscanf(inpFile,"%lf",&ip.prodParam[i][j]);
      }
      ignore_i = fscanf(inpFile,"%lf",&ip.massFluxes[i]);
      if(ip.massFluxes[i]<0.0){
        sprintf(errMessage,"Bad value for mass flux #%d.",i);
        plas_TerminateOnError(errMessage);
      }
      ignore_cp = fgets(inpText,200,inpFile);
    }
  }

  //***Read diameter spectrum of dispersed entities***//

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&ip.iniDiamType);
  if(ip.iniDiamType<0 || ip.iniDiamType>2){
    sprintf(errMessage,"Bad value for initial diameter distribution type.");
    plas_TerminateOnError(errMessage);
  }
  ignore_i = fscanf(inpFile,"%lf",&ip.iniDiam);
  if(ip.iniDiam<1e-6){
    sprintf(errMessage,"Bad value for entity diameter.");
    plas_TerminateOnError(errMessage);
  }
  ignore_i = fscanf(inpFile,"%lf",&ip.iniDiamStd);
  if(ip.iniDiamStd<0.0){
    sprintf(errMessage,"Negative value for entity diameter standard deviation.");
    plas_TerminateOnError(errMessage);
  }
  ignore_cp = fgets(inpText,200,inpFile);

  //***Read initial velocity of dispersed entities***//

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);
  for(i=0; i<3; i++){
    ignore_i = fscanf(inpFile,"%lf",&ip.iniVel[i]);
  }
  ignore_cp = fgets(inpText,200,inpFile);

  //***Read temperature of dispersed entities***//

  ip.iniTempDisp = plas_ReadDoubleParam(inpFile);
  if(ip.iniTempDisp<0.0){
    sprintf(errMessage,"Bad value for dispersed phase temperature.");
    plas_TerminateOnError(errMessage);
  }

  //***Read entity material***//

  ip.material = plas_ReadIntParam(inpFile);

  if(ip.material==MAT_COPPER
     || ip.material==MAT_POLY){rp.flowType = FLOW_PARTIC;}
  else if(ip.material==MAT_WATER
     || ip.material==MAT_NHEPTANE){rp.flowType = FLOW_DROPLET;}
  else if(ip.material==MAT_HYDROGEN
     || ip.material==MAT_OXYGEN
     || ip.material==MAT_AIR){rp.flowType = FLOW_BUBBLY;}
  else{
    sprintf(errMessage,"No match in the entity material database.");
    plas_TerminateOnError(errMessage);
  }

  //***Read momentum back-coupling option***//

  ip.momentumCoupl = plas_ReadIntParam(inpFile);
  if(ip.momentumCoupl<0 && ip.momentumCoupl>2){
    sprintf(errMessage,"Bad value for momentum coupling parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read volume fraction back-coupling option***//

  ip.volfracCoupl = plas_ReadIntParam(inpFile);
  if(ip.volfracCoupl!=0 && ip.volfracCoupl!=1){
    sprintf(errMessage,"Bad value for mass coupling parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read energy coupling option***//

  ip.energyCoupl = plas_ReadIntParam(inpFile);
  if(ip.energyCoupl!=0 && ip.energyCoupl!=1){
    sprintf(errMessage,"Bad value for energy coupling parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read collision model option***//

  ip.collisionModel = plas_ReadIntParam(inpFile);
  if(ip.collisionModel<0 || ip.collisionModel>2){
    sprintf(errMessage,"Bad value for collision model parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read slip-shear lift force option***//

  ip.liftForce = plas_ReadIntParam(inpFile);
  if(ip.liftForce!=0 && ip.liftForce!=1){
    sprintf(errMessage,"Bad value for lift force parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read evaporation model option***//

  ip.evapModel = plas_ReadIntParam(inpFile);
  if(ip.evapModel!=0 && ip.evapModel!=1){
    sprintf(errMessage,"Bad value for evaporation model parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read saturation model option***//

  ip.saturModel = plas_ReadIntParam(inpFile);
  if(ip.saturModel!=0 && ip.saturModel!=1){
    sprintf(errMessage,"Bad value for saturation model parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read periodic boundaries option***//

  ip.perBnd = plas_ReadIntParam(inpFile);
  if(ip.perBnd!=0 && ip.perBnd!=1){
    sprintf(errMessage,"Bad value for periodic boundary parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read gravity vector***//

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);
  for(i=0; i<3; i++){
    ignore_i = fscanf(inpFile,"%lf",&ip.gravVec[i]);
  }
  ignore_cp = fgets(inpText,200,inpFile);

  //***Read restart flag***//

  ip.restart = plas_ReadIntParam(inpFile);
  if(ip.restart!=0 && ip.restart!=1){
    sprintf(errMessage,"Bad value for restart flag.");
    plas_TerminateOnError(errMessage);
  }

  //***Read the output statistics filename***//

  ip.writeStatsFilename = new char[200];
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%s",inpText);
  strncpy(ip.writeStatsFilename,inpText,200);
  ignore_cp = fgets(inpText,200,inpFile);
  {
    sprintf(errMessage,"Output statistics filename: \"%s\"\n",ip.writeStatsFilename);
    screenOutput(errMessage);
  }

  //***Read the output tecplot filename***//

  ip.writeTecplotFilename = new char[200];
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%s",inpText);
  strncpy(ip.writeTecplotFilename,inpText,200);
  ignore_cp = fgets(inpText,200,inpFile);
  {
    sprintf(errMessage,"Output tecplot filename: \"%s\"\n",ip.writeTecplotFilename);
    screenOutput(errMessage);
  }

  //***Read the configuration filename***//

  ip.confFilename = new char[200];
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%s",inpText);
  strncpy(ip.confFilename,inpText,200);
  ignore_cp = fgets(inpText,200,inpFile);
  {
    sprintf(errMessage,"Configuration filename: \"%s\"\n",ip.confFilename);
    screenOutput(errMessage);
  }

  fclose(inpFile);
}


void plas::plas_SearchDomainParallel(LOCAL_ENTITY_VARIABLES *ent)
{
  int bfCtr = 0;
  int elmFoundAnyProc = 0;
  int idx,elmFound,procFound;
  double dist;

  //***Perform a successive neighbour search (find by chance)***//

  do{
    ent->elm = plas_RandomInteger(0,fp.numElm-1);
    bfCtr++;
    plas_SearchSuccessive(ent);
    if(ent->flag==DFLAG_ENABLED){
      elmFound = 1;
      procFound = plas_MpiGetRank();
    }else{
      elmFound = 0;
      procFound = plas_MpiGetNumProc();
    }

  } while(!elmFound && bfCtr<5);

  //***Determins if an element has been found on any process***//

  elmFoundAnyProc = plas_MpiAllMaxInt(elmFound);

  //***Brute force search (relevant elements given by flow solver)***//

  if(!elmFoundAnyProc){
    ent->elm = StartElementSearch(ent->pos);
    do{
      plas_SetElementGeometry(fp.numDim,ent);
      plas_FindMinimumElementFaceDistance(fp.numDim,ent,&idx,&dist);
      if(dist>(0.0-rp.errTol)){
        elmFound = 1;
        procFound = plas_MpiGetRank();
      }else{
        elmFound = 0;
        procFound = plas_MpiGetNumProc();
        ent->elm++;
      }
    } while(!elmFound && ent->elm<=EndElementSearch(ent->pos));
  }

  //***In case an element is found on more than opne process, take the minimum one***//

  procFound = plas_MpiAllMinInt(procFound);

  //***Disable entity in case of not found***//

  if(elmFound && plas_MpiGetRank()==procFound){
    ent->flag = DFLAG_ENABLED;
    ent->node = plas_FindNearestElementNode(ent);
  } else{
    ent->flag = DFLAG_DISABLED;
    ent->elm = -1;
    ent->node = -1;
  }
}


void plas::plas_SearchSuccessive(LOCAL_ENTITY_VARIABLES *ent)
{
  int elmFound = 0;
  int leftDomain = 0;
  int lastElm = ent->elm;
  int lastNode = ent->node;
  int neighbourElm,idx;
  double dist;

  //***Successive neighbour search***//

  do{

    //***Set geometry of current element***//

    plas_SetElementGeometry(fp.numDim,ent);

    //***Find the face with the minimum distance to the entity***//

    plas_FindMinimumElementFaceDistance(fp.numDim,ent,&idx,&dist);

    if(dist>(0.0-rp.errTol)){

      //***In case the minimum distance is positive, the element is found***//

      elmFound = 1;
    } else{

      //***Search the neighbour element in directin of the minimum face distance***//

      neighbourElm = getElmNeighbour(ent->elm,idx);

      if(neighbourElm==-1){
        leftDomain = 1;
      } else{
        ent->elm = neighbourElm;
        if(lastElm!=-1){
          lastElm = neighbourElm;
        }
      }
    }
  } while(!elmFound && !leftDomain);

  //***Disable entity in case of not found***//

  if(elmFound){
    ent->flag = DFLAG_ENABLED;
    ent->node = plas_FindNearestElementNode(ent);
  } else{
    ent->flag = DFLAG_LEFT;
    ent->elm = lastElm;
    ent->node = lastNode;
  }
}


double plas::plas_SetDiameter()
{
  double muN = ip.iniDiam;
  double sigN = ip.iniDiamStd;
  double diameter=0.,muLn,sigLn;

  //***Compute diameter***//

  do{

    if(ip.iniDiamType==0){

      //***Constant diameter***//

      diameter = muN;

    } else if(ip.iniDiamType==1){

      //***Gaussian normal distribution***//

      diameter = plas_RandomGaussian(muN,sigN);

    } else if(ip.iniDiamType==2){

      //***Log-normal distribution***//

      sigLn = sqrt(log((sigN*sigN/(muN*muN))+1.0));
      muLn = log(muN)-0.5*sigLn*sigLn;
      diameter = exp(plas_RandomGaussian(muLn,sigLn));
    }

  } while(diameter<rp.errTol);

  return diameter;
}


void plas::plas_SetElementFaces(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  int type = getElementType(ent->elm);
  int ifac,idim;

  //***Set number of faces accrding to element type***//

  if(type==ELM_SIMPLEX){
    ent->edata.numElmFaces = numDim+1;
  } else if(type==ELM_PRISM){
    ent->edata.numElmFaces = 5;
  } else if(type==ELM_QUAD){
    ent->edata.numElmFaces = 4;
  } else if(type==ELM_HEX){
    ent->edata.numElmFaces = 6;
  } else if(type==ELM_PYRAMID){
    ent->edata.numElmFaces = 5;
  }

  //***Get faces and normals from flow solver***//

  for(ifac=0; ifac<ent->edata.numElmFaces; ifac++){
    for(idim=0; idim<numDim; idim++){
      ent->edata.elmFaceVectors[ifac][idim] = getElmFaceMiddlePoint(ent->elm,ifac,idim);
      ent->edata.elmNorms[ifac][idim] = getElmNormComp(ent->elm,ifac,idim);
    }
  }
}


void plas::plas_SetElementGeometry(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  plas_SetElementNodes(numDim,ent);
  plas_SetElementFaces(numDim,ent);
}


void plas::plas_SetElementNodes(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  int type = getElementType(ent->elm);
  int inod;

  //***Set number of nodes accrding to element type***//

  if(type==ELM_SIMPLEX){
    ent->edata.numElmNodes = numDim+1;
  } else if(type==ELM_PRISM){
    ent->edata.numElmNodes = 6;
  } else if(type==ELM_QUAD){
    ent->edata.numElmNodes = 4;
  } else if(type==ELM_HEX){
    ent->edata.numElmNodes = 8;
  } else if(type==ELM_PYRAMID){
    ent->edata.numElmNodes = 5;
  }

  //***Get node elements from flow solver***//

  for(inod=0; inod<ent->edata.numElmNodes; inod++){
    ent->edata.elmNodes[inod] = getElmNode(ent->elm,inod);
  }
}


void plas::plas_SolveGaussSeidel(int numDim, double **mat, double *s, double *rhs)
{
  int idim,jdim;
  double errLoc,errMax,sOld[numDim];

  for(idim=0; idim<numDim; idim++){
    sOld[idim] = s[idim];
  }

  do{
    errMax = 0.0;
    for(idim=0; idim<numDim; idim++){
      s[idim] = 0.0;
      for(jdim=0; jdim<numDim; jdim++){
        if(jdim<idim){
          s[idim] += mat[idim][jdim]*s[jdim];
        }
        if(jdim>idim){
          s[idim] += mat[idim][jdim]*sOld[jdim];
        }
      }
      s[idim] = (rhs[idim]-s[idim])/mat[idim][idim];
      errLoc = std::abs(s[idim]-sOld[idim]);
      if(errLoc>errMax){errMax = errLoc;}
      sOld[idim] = s[idim];
    }
  }while(errMax>1e-6);
}


void plas::plas_TerminateOnError(const std::string& errMessage)
{
  screenWarning(errMessage);
  screenWarning("FATAL ERROR!");
  throw 42;
}


void plas::plas_WallBounce(int numDim, double elasticity, LOCAL_ENTITY_VARIABLES *ent, int ibnd, int ifac)
{
  int idim;
  double unitVec[numDim];
  double wallDistancePositionVector,wallDistanceVelocityVector;

  //***Calculate unit normal of boundary segment***//

  plas_CalcBoundaryUnitNormal(numDim,ibnd,ifac,unitVec);

  //***Calculate normal wall distance of entity position***//

  wallDistancePositionVector = plas_CalcWallFaceDistance(numDim,ent->pos,ibnd,ifac);

  //***Adapt velocity due to elasticity factor***//

  for(idim=0; idim<numDim; idim++){
    ent->vel[idim] *= elasticity;
  }

  //***Add position vector to velocity vector***//

  for(idim=0; idim<numDim; idim++){
    ent->vel[idim] += ent->pos[idim];
  }

  //***Calculate normal wall distance of entity velocity vector***//

  wallDistanceVelocityVector = plas_CalcWallFaceDistance(numDim,ent->vel,ibnd,ifac);

  //***Mirror position and velocity vectors on wall segment***//

  for(idim=0; idim<numDim; idim++){
    ent->pos[idim] = ent->pos[idim]-2.0*wallDistancePositionVector*unitVec[idim];
    ent->vel[idim] = ent->vel[idim]-2.0*wallDistanceVelocityVector*unitVec[idim]-ent->pos[idim];
  }
}


void plas::plas_WriteStatsFile(char *outpString, int iter, double time)
{
  FILE *outpFile;
  if (!plas_MpiGetRank()) {
    outpFile = fopen(outpString,"a");
    if (outpFile==NULL)
      return;

    fprintf(outpFile,"%5d\t%6.2f\t%6d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%11.4e\t%11.4e\t%11.4e\t%6.2f\n",
    iter,time,
    sd.enabled,sd.in,sd.out,sd.bounce,sd.coll,
    sd.periodic,sd.passed,sd.lost,sd.reynoldsAvg,sd.nusseltAvg,
    sd.dtLagrAvg,sd.subIterAvg);
    fclose(outpFile);
  }
}


void plas::plas_WriteTecplotFile(char *outpString, int iter, double time)
{
  FILE *outpFile;
  int ient, iproc, idim;
  const int numProc = plas_MpiGetNumProc();
  const int irank   = plas_MpiGetRank();

  for (iproc=0; iproc<numProc; ++iproc) {
    if (iproc==irank) {
      outpFile = fopen(outpString,"a");
      if (irank==0) {

        // write header (if no entities are active, write dummy)
        fprintf(outpFile,"ZONE T=\"Entities (dispersed)\" I=%d AUXDATA ITER=\"%d\" SOLUTIONTIME=%12.6f DATAPACKING=POINT\n",sd.enabled? sd.enabled:1,iter,time);
        if (!sd.enabled) {
          for (idim=0; idim<fp.numDim*2+3; ++idim)
            fprintf(outpFile,"0 ");
          fprintf(outpFile,"\n");
        }

      }
      for (ient=0; ient<ip.numMaxEnt && sd.enabled; ++ient) {
        if (ed[ient].flag==DFLAG_ENABLED || ed[ient].flag==DFLAG_CREATED) {

          PLAS_ENTITY_DATA *e = &(ed[ient]);
          for (idim=0; idim<fp.numDim; ++idim)
            fprintf(outpFile,"%.12f ",e->position[idim]);
          for (idim=0; idim<fp.numDim; ++idim)
            fprintf(outpFile,"%.12f ",e->velocity[idim]);
          fprintf(outpFile,"%.12f %.12f %.12f\n",e->diameter,e->temperature,0.);

        }
      }
      fclose(outpFile);
    }
    plas_MpiBarrier();
  }
}

