
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "plas.h"


void plas::initialize(const XMLNode& x)
{
  int ient,inod;


  //***Flow solver parameters that have to be set only once***//
  setFlowSolverParamOnInit(&fp);


  //***Read parameters from input file***//
  plas_ReadParameters(x);


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
  FILE *inpFile = fopen(ip.writeTecplotFilename.c_str(),"r");
  if (ip.writeTecplotFilename.length() && (inpFile==NULL || !ip.restart)) {
    ip.restart = 0;
    plas_CreateTecplotFile(ip.writeTecplotFilename);
  }
  if (inpFile)
    fclose(inpFile);

  if (ip.writeStatsFilename.length())
    plas_CreateStatsFile(ip.writeStatsFilename);

  /* loaded or random distribution */
  if (ip.restart)
    plas_LoadInitialDistribution(ip.writeTecplotFilename);
  else if (ip.numIniEnt>0)
    plas_RandomInitialDistribution();
}


void plas::run()
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  char errMessage[100];
  int idx,ient,jent,inod,idim,iunk,ibnd,ifac,subIter,facFound,avgctr;
  double cellSize,minCellSize,entVel,dtLagr,dtRemaining;


  //***Flow solver parameters that have to be set at every time step***//
  screenOutput("Initialization...");
  setFlowSolverParamOnTimeStep(&(fp));


  //***Allocate dynamic memory***//
  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);

  //***Initializations***//

  sd.enabled     = 0;
  sd.in          = 0;
  sd.out         = 0;
  sd.bounce      = 0;
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

  //***Entity production***//

  if (ip.numProdDom>0){
    screenOutput("Imposing new entities in production domains...");
    plas_ImposeProductionDomains();
  }

  if (numExtEnt>0){
    screenOutput("Imposing externally generated entities...");
    plas_ImposeExternal();
  }

  //***Loop over dispersed entities to update trajectories***//

  screenOutput("Updating trajectories...");
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



//***Pass entities between parallel processes***//
  plas_PassEntities();


//***Update counters***//
  screenOutput("Post-processing...");
  sd.enabled  -= sd.out;
  sd.lost     += sd.leftproc-sd.passed;
  sd.enabled  -= sd.lost;


  //***Update statistics***//
  if(sd.enabled>0){
    sd.reynoldsAvg /= avgctr;
    sd.nusseltAvg  /= avgctr;
    sd.dtLagrAvg   /= avgctr;
    sd.subIterAvg  /= avgctr;
  } else{
    sd.reynoldsAvg = 0.;
    sd.nusseltAvg  = 0.;
    sd.dtLagrAvg   = 0.;
    sd.subIterAvg  = 0.;
  }

  //**Compute dispersed phase cellwise data***//

  plas_CalcCellwiseData();

  //***Write PLaS output to files***//

  if (ip.writeStatsFilename.length())
    plas_WriteStatsFile(ip.writeStatsFilename,fp.iter,fp.time);
  if (ip.writeTecplotFilename.length())
    plas_WriteTecplotFile(ip.writeTecplotFilename,fp.iter,fp.time);

  //***Free dynamic memory***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(fp.numDim,&flow);
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
  // get normal vector of a boundary face
  for (int idim=0; idim<numDim; ++idim)
    unitVec[idim] = getBndFaceNormComp(ibnd,ifac,idim);

  // divide normal vector components by vector length
  const double length = plas_CalcVectorLength(numDim,unitVec);
  for (int idim=0; idim<numDim; ++idim)
    unitVec[idim] /= length;
}


void plas::plas_CalcCellwiseData()
{
  int ient,inod,ielm,idim,jdim,itype,numElmNodes=0;
  double ivol=0,xi,yi,xixi,xiyi,volFracOld[fp.numNod];


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
        numElmNodes = (itype==ELM_SIMPLEX? fp.numDim+1 :
                      (itype==ELM_PRISM?   6 :
                      (itype==ELM_QUAD?    4 :
                      (itype==ELM_HEX?     8 :
                      (itype==ELM_PYRAMID? 5 :
                                           0 )))));

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


void plas::plas_CreateStatsFile(const std::string &outpString)
{
  FILE *outpFile;
  outpFile = fopen(outpString.c_str(),"w");
  if (outpFile==NULL)
    return;

  fprintf(outpFile,"Iter\tTime\tEnt\tIn\tOut\tBounc\tColl\tPass\tLost\tRe_disp\tNu_disp\tdt_Lagr\tsubit\n");
  fclose(outpFile);
}


void plas::plas_CreateTecplotFile(const std::string &outpString)
{
  FILE* outpFile;
  outpFile = fopen(outpString.c_str(),"w");
  if (fp.numDim==2) {
    fprintf(outpFile,"VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"d\" \"T\" \"theta\"\n");
  } else if(fp.numDim==3) {
    fprintf(outpFile,"VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"d\" \"T\" \"theta\"\n");
  }
  fclose(outpFile);
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
      if(ent->elm==getBndDomElm(ibnd,ifac)){
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

  return node;
}


void plas::plas_ImposeExternal()
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  int ient,idim;
  double *newPos,*newVel,*newDiam,*newTemp;


  //***Allocation of local data structure***//
  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);


  //***Broadcast entity data***//
  newDiam = new double[numExtEnt];
  newTemp = new double[numExtEnt];
  newPos  = new double[numExtEnt*fp.numDim];
  newVel  = new double[numExtEnt*fp.numDim];

  newDiam = extEntDiam;
  newTemp = extEntTemp;
  newPos = extEntPos;
  newVel = extEntVel;

  //***Loop over entities to generate***//

  for(ient=0; ient<numExtEnt; ient++){

    for(idim=0; idim<fp.numDim; idim++){
      ent.pos[idim] = newPos[fp.numDim*ient+idim];
    }
    for(idim=0; idim<fp.numDim; idim++){
      ent.vel[idim] = newVel[fp.numDim*ient+idim];
    }
    ent.diam = newDiam[ient];
    ent.temp = newTemp[ient];

    //***Element search***//

    plas_SearchSuccessive(&ent);

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


  //***Creation of bubbles***//
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


  //***Generate entities***//
  for(ient=0; ient<bCtr; ient++){

    ent.elm = 0;
    ent.diam = newDiam[ient];
    ent.temp = ip.iniTempDisp;
    for(idim=0; idim<fp.numDim; idim++){
      ent.pos[idim] = newPos[fp.numDim*ient+idim];
      ent.vel[idim] = ip.iniVel[idim];
    }

    //***Element search***//

    plas_SearchSuccessive(&ent);

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
  sprintf(msg,"imposed entities: %d",bCtr);
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
  double impFac[ent->edata.numElmNodes];

  // compute impact factors of element nodes
  plas_CalcNodeImpactFactors(ent,impFac);

  // flow temperature
  flow->temp = 0.;
  for (int idim=0; idim<ent->edata.numElmNodes; ++idim) {
    flow->temp += impFac[idim]*((1.0-step)*getTemperatureOld(ent->edata.elmNodes[idim])
                                    + step*getTemperature(ent->edata.elmNodes[idim]));
  }
}


void plas::plas_InterpolateVelocity(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim,jdim,kdim;
  double impFac[ent->edata.numElmNodes];

  // compute impact factors of element nodes
  plas_CalcNodeImpactFactors(ent,impFac);

  // flow velocity
  for(idim=0; idim<fp.numDim; idim++){
    flow->vel[idim] = 0.0;
    for(jdim=0; jdim<ent->edata.numElmNodes; jdim++){
      flow->vel[idim] += impFac[jdim] * (
        (1. - step) * getVelocityCompOld(ent->edata.elmNodes[jdim],idim)
            + step  * getVelocityComp(ent->edata.elmNodes[jdim],idim) );
    }
  }

  // flow velocity space derivatives
  for(idim=0; idim<fp.numDim; idim++){
    for(jdim=0; jdim<fp.numDim; jdim++){
      flow->velDx[idim][jdim] = 0.0;
      for(kdim=0; kdim<ent->edata.numElmNodes; kdim++){
        flow->velDx[idim][jdim] += impFac[kdim]*((1.0-step)*getVelocityDerivativeCompOld(ent->edata.elmNodes[kdim],idim,jdim)
                                                     + step*getVelocityDerivativeComp(ent->edata.elmNodes[kdim],idim,jdim));
      }
    }
  }

  // flow velocity time derivative
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


void plas::plas_LoadInitialDistribution(const std::string &inpString)
{
  LOCAL_ENTITY_VARIABLES ent;
  int ient,jent,idim,numIni;
  char buffer[200],*search;
  long fpos;
  FILE *inpFile;

  //***Read data file***//

  inpFile = fopen(inpString.c_str(),"r");

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

    plas_SearchSuccessive(&ent);

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
  int *disps         = new int[1];
  int *recs          = new int[1];
  int leftCtrOneProc = sd.leftproc;
  int *foundProc     = new int[leftCtrOneProc];
  int *foundElm      = new int[leftCtrOneProc];
  int *foundNode     = new int[leftCtrOneProc];
  int idim,jdim,ient,jent,ielm,inod;

  double *leftDataOneProc = new double[(2*fp.numDim+1)*leftCtrOneProc];


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


  //***Search for all free entities on all processes***//
  idim = 0;
  for(ient=0; ient<leftCtrOneProc; ient++){

    for(jdim=0; jdim<fp.numDim; jdim++){
      ent.pos[jdim] = leftDataOneProc[idim+jdim];
    }


    //***Element search***//
    plas_SearchSuccessive(&ent);

    if(ent.flag==DFLAG_ENABLED){
      foundProc[ient] = 0;
      foundElm[ient] = ent.elm;
      foundNode[ient] = ent.node;
    } else{
      foundProc[ient] = -1;
    }
    idim += 2*fp.numDim+1;
  }


  //***Broadcast information of found entities***//
  idim = 0;
  jent = 0;

  for(ient=0; ient<leftCtrOneProc; ient++){

    if(foundProc[ient]>-1){
      ielm = foundElm[ient];
      inod = foundNode[ient];
    }

    if(foundProc[ient]!=0){
      idim += 2*fp.numDim+1;
      continue;
    }

    while(ed[jent].flag==DFLAG_ENABLED){
      jent++;
      if(jent>=ip.numMaxEnt)
        plas_TerminateOnError("memory exceeded when passing trajectory information");
    }


    //***Generate entity on the new process***//
    sd.passed++;
    ed[jent].flag = DFLAG_ENABLED;
    ed[jent].element = ielm;
    ed[jent].node = inod;

    for(jdim=0; jdim<fp.numDim; jdim++){
      ed[jent].position[jdim] = leftDataOneProc[idim];
      idim++;
    }

    for(jdim=0; jdim<fp.numDim; jdim++){
      ed[jent].velocity[jdim] = leftDataOneProc[idim];
      idim++;
    }

    ed[jent].diameter = leftDataOneProc[idim];
    idim++;
  }


  //***De-allocation of local data structure***//
  plas_DeallocateLocalEntityVar(&ent);

  delete[] foundProc;
  delete[] foundElm;
  delete[] foundNode;
  delete[] leftDataOneProc;
  delete[] disps;
  delete[] recs;
}


double plas::plas_RandomDouble()
{
  int r = plas_RandomInteger(0,100000);
  return (double)(r/100000.);
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
  int ient,idim;
  double p,q,rand1,rand2,rand3,rand4,setDiam,elmVolume;

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(fp.numDim,&ent);
  plas_AllocateLocalFlowVar(fp.numDim,&flow);

  //***Initializations***//
  ip.numIniEnt = (ip.numIniEnt>ip.numMaxEnt? ip.numMaxEnt : ip.numIniEnt);


  //***Loop over entities to generate***//

  for(ient=0; ient<ip.numIniEnt; ient++){

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


void plas::plas_ReadParameters(const XMLNode& x)
{
  // read parameters
  ip.numMaxEnt  = std::max(0, x.getAttribute< int    >("entities.max",1000));
  ip.numIniEnt  = std::max(0, x.getAttribute< int    >("entities.init",0));
  ip.iniDiam    = std::max(0.,x.getAttribute< double >("entities.production.diameter.mean",1.));
  ip.iniDiamStd = std::max(0.,x.getAttribute< double >("entities.production.diameter.std", 0.));

  ip.iniVel[0] = x.getAttribute< double >("entities.production.velocity.x",0.);
  ip.iniVel[1] = x.getAttribute< double >("entities.production.velocity.y",0.);
  ip.iniVel[2] = x.getAttribute< double >("entities.production.velocity.z",0.);
  ip.iniTempDisp = std::max(0.,x.getAttribute< double >("entities.production.temperature",293.15));

  const std::string
    e_material  = x.getAttribute< std::string >("entities.material","air"),
    e_ddist     = x.getAttribute< std::string >("entities.production.diameter.distribution","constant"),
    c_momentum  = x.getAttribute< std::string >("couple.momentum","no"),
    c_energy    = x.getAttribute< std::string >("couple.energy",  "no"),
    m_collision = x.getAttribute< std::string >("model.collision","no");

  ip.volfracCoupl = (x.getAttribute< std::string >("couple.volumefraction",             "no")=="yes");
  ip.evapModel    = (x.getAttribute< std::string >("model.droplets.thinfilmevaporation","no")=="yes");
  ip.saturModel   = (x.getAttribute< std::string >("model.bubbles.saturation",          "no")=="yes");
  ip.liftForce    = (x.getAttribute< std::string >("model.bubbles.slip-shearliftforce", "no")=="yes");

  ip.gravVec[0] = x.getAttribute< double >("gravity.x",0.);
  ip.gravVec[1] = x.getAttribute< double >("gravity.y",0.);
  ip.gravVec[2] = x.getAttribute< double >("gravity.z",0.);

  ip.restart = (x.getAttribute< std::string >("restart","no")=="yes");
  ip.writeStatsFilename   = x.getAttribute< std::string >("output.statistics","plas.txt");
  ip.writeTecplotFilename = x.getAttribute< std::string >("output.results","plas.plt");

  const std::string label = x.getAttribute< std::string >("label","output");
  bool
    write_stats   = !(x.getAttribute< std::string >("output.statistics","yes")=="no"),
    write_tecplot = !(x.getAttribute< std::string >("output.results",   "yes")=="no");
  ip.writeStatsFilename   = (write_stats?
                            (label.length()? label + ".plas.txt" : std::string("output.plas.txt"))
                          :  std::string());
  ip.writeTecplotFilename = (write_tecplot?
                            (label.length()? label + ".plas.plt" : std::string("output.plas.plt"))
                          :  std::string());


  // apply corrections and set additional paramters
  ip.material = (e_material=="Cu" || e_material=="copper"?    MAT_COPPER   :
                (e_material=="C8H8"?                          MAT_POLY     :
                (e_material=="H2O"|| e_material=="water"?     MAT_WATER    :
                (                    e_material=="n-Heptane"? MAT_NHEPTANE :
                (e_material=="H2" || e_material=="hydrogen"?  MAT_HYDROGEN :
                (e_material=="O2" || e_material=="oxygen"?    MAT_OXYGEN   :
                (                    e_material=="air"?       MAT_AIR      :
                                                              MAT_AIR )))))));
  rp.flowType = (ip.material==MAT_COPPER   || ip.material==MAT_POLY?                           FLOW_PARTIC  :
                (ip.material==MAT_WATER    || ip.material==MAT_NHEPTANE?                       FLOW_DROPLET :
                (ip.material==MAT_HYDROGEN || ip.material==MAT_OXYGEN || ip.material==MAT_AIR? FLOW_BUBBLY  :
                  FLOW_BUBBLY )));
  ip.iniDiamType = (e_ddist=="constant"?   0 :
                   (e_ddist=="normal"?     1 :
                   (e_ddist=="log-normal"? 2 : 0 )));

  ip.momentumCoupl  = (c_momentum=="PIC"?        1 :
                      (c_momentum=="projection"? 2 : 0 ));
  ip.energyCoupl    = (c_energy=="yes" || c_energy=="one-way");
  ip.collisionModel = (m_collision=="uncorrelated"? 1 :
                      (m_collision=="correlated"?   2 : 0 ));

  if ((ip.numProdDom = x.nChildNode("production"))) {
    ip.prodDom    = new int    [ip.numProdDom];
    ip.massFluxes = new double [ip.numProdDom];
    ip.prodParam  = new double*[ip.numProdDom];
    for (int i=0; i<ip.numProdDom; ++i)
      ip.prodParam[i] = new double[6];
  }
  for (int i=0; i<ip.numProdDom; ++i) {
    XMLNode p = x.getChildNode("production",i);
    ip.prodDom[i] = (p.getAttribute< std::string >("type")=="line"?      1 :
                    (p.getAttribute< std::string >("type")=="rectangle"? 2 :
                    (p.getAttribute< std::string >("type")=="ellipse"?   3 : 0)));
    ip.massFluxes[i] = std::max(p.getAttribute< double >("massflux",0.),0.);
    ip.prodParam[i][0] = p.getAttribute< double >("x0",0.);
    ip.prodParam[i][1] = p.getAttribute< double >("y0",0.);
    ip.prodParam[i][2] = p.getAttribute< double >("z0",0.);
    ip.prodParam[i][3] = p.getAttribute< double >("x1",0.);
    ip.prodParam[i][4] = p.getAttribute< double >("y1",0.);
    ip.prodParam[i][5] = p.getAttribute< double >("z1",0.);
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
    ent->elm  = lastElm;
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
  // set number of faces according to element type
  const int type = getElementType(ent->elm);
  ent->edata.numElmFaces = (type==ELM_SIMPLEX? numDim+1 :
                           (type==ELM_PRISM?   5 :
                           (type==ELM_QUAD?    4 :
                           (type==ELM_HEX?     6 :
                           (type==ELM_PYRAMID? 5 : 0 )))));

  // get faces and normals
  for (int ifac=0; ifac<ent->edata.numElmFaces; ++ifac) {
    for (int idim=0; idim<numDim; ++idim) {
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
  // set number of nodes according to element type
  const int type = getElementType(ent->elm);
  ent->edata.numElmNodes = (type==ELM_SIMPLEX? numDim+1 :
                           (type==ELM_PRISM?   6 :
                           (type==ELM_QUAD?    4 :
                           (type==ELM_HEX?     8 :
                           (type==ELM_PYRAMID? 5 : 0 )))));

  // get node elements
  for (int inod=0; inod<ent->edata.numElmNodes; ++inod)
    ent->edata.elmNodes[inod] = getElmNode(ent->elm,inod);
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


void plas::plas_WriteStatsFile(const std::string &outpString, int iter, double time)
{
  FILE *outpFile;
    outpFile = fopen(outpString.c_str(),"a");
    if (outpFile==NULL)
      return;

    fprintf(outpFile,"%5d\t%6.2f\t%6d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%11.4e\t%11.4e\t%11.4e\t%6.2f\n",
    iter,time,
    sd.enabled,sd.in,sd.out,sd.bounce,sd.coll,
    sd.passed,sd.lost,sd.reynoldsAvg,sd.nusseltAvg,
    sd.dtLagrAvg,sd.subIterAvg);
    fclose(outpFile);
}


void plas::plas_WriteTecplotFile(const std::string &outpString, int iter, double time)
{
  FILE *outpFile;
  int ient, idim;

  outpFile = fopen(outpString.c_str(),"a");

  // write header (if no entities are active, write dummy)
  fprintf(outpFile,"ZONE T=\"Entities (dispersed)\" I=%d AUXDATA ITER=\"%d\" SOLUTIONTIME=%12.6f DATAPACKING=POINT\n",sd.enabled? sd.enabled:1,iter,time);
  if (!sd.enabled) {
    for (idim=0; idim<fp.numDim*2+3; ++idim)
      fprintf(outpFile,"0 ");
    fprintf(outpFile,"\n");
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
