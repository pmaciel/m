
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <set>
#include <list>

#include "boost/progress.hpp"
#include "mfactory.h"
#include "pillaz.h"


/// Forward declaration: Pillaz materials database initialisation
void pillaz_material_database();


void pillaz::initialize(const XMLNode& x)
{
  // utility variables
  long unsigned count = 0;
  boost::progress_display pbar(count);


  screenOutput("initializing Pillaz...");


  screenOutput("initialize materials database...");
  pillaz_material_database();
  screenOutput("initialize materials database.");


  screenOutput("set flow solver parameters...");
  setFlowSolverParamOnInit(&fp);
  screenOutput("set flow solver parameters.");


  screenOutput("recreating node-to-element connectivity...");
  pillaz_ReadParameters(x);
  screenOutput("recreating node-to-element connectivity.");


  screenOutput("recreating map of nodes to (inner) elements...");

  // allocate
  m_nodetoelem.resize(fp.numNod);

  // push back all the elements to their composing node indices
  for (size_t iz=0; iz<fp.nInnerElements.size(); ++iz) {
    for (size_t ie=0; ie<fp.nInnerElements[iz]; ++ie) {
      std::vector< int > en( pillaz_getElmNNodes(getElmType(iz,ie)) );
      getElmNodes(iz,ie,&en[0]);
      for (std::vector< int >::const_iterator n=en.begin(); n!=en.end(); ++n)
        m_nodetoelem[ *n ].push_back(pillaz_elementaddress(iz,ie));
    }
  }
  screenOutput("recreating map of nodes to (inner) elements.");


  screenOutput("recreating map of (inner) elements to elements (sharing a face)...");

  // allocate
  count = 0;
  m_elemtoelem.resize(fp.nInnerElements.size());
  for (size_t iz=0; iz<fp.nInnerElements.size(); ++iz) {
    if (fp.nInnerElements[iz])
      m_elemtoelem[iz].resize(fp.nInnerElements[iz]);
    for (size_t ie=0; ie<fp.nInnerElements[iz]; ++ie) {
      count += pillaz_getElmNFaces(getElmType(iz,ie));
      m_elemtoelem[iz][ie].resize( pillaz_getElmNFaces(getElmType(iz,ie)) );
    }
  }

  // hold information about elements on a boundary
  std::list< pillaz_elementaddress > list_boundelems;

  // perform the (inner) element to element search
  pbar.restart(count);
  for (size_t iz=0; iz<m_elemtoelem.size(); ++iz) {
    for (size_t ie=0; ie<m_elemtoelem[iz].size(); ++ie) {
      for (size_t iface=0; iface<m_elemtoelem[iz][ie].size(); ++iface, ++pbar) {

        // get element face nodes
        std::vector< int > fn( pillaz_getElmNNodes(pillaz_getElmFaceType(getElmType(iz,ie),iface)), -1 );
        pillaz_getElmFaceNodes(iz,ie,iface,&fn[0]);
        std::sort(fn.begin(),fn.end());

        // create set of searcheable elements (those sharing a node)
        std::set< std::pair< int,int > > neighelems;
        for (std::vector< int >::const_iterator n=fn.begin(); n!=fn.end(); ++n)
          for (std::vector< pillaz_elementaddress >::const_iterator e=m_nodetoelem[*n].begin(); e!=m_nodetoelem[*n].end(); ++e)
            neighelems.insert( neighelems.end(), std::pair< int,int >(e->izone,e->ielem));
        neighelems.erase(std::pair< int,int >(iz,ie));  // avoid self-search

        // find the element with a face with the same nodes
        bool isneighborfacefound = false;
        for (std::set< std::pair< int,int > >::const_iterator ne=neighelems.begin(); ne!=neighelems.end() && !isneighborfacefound; ++ne) {
          for (int jface=0; jface<pillaz_getElmNFaces(getElmType(ne->first,ne->second)) && !isneighborfacefound; ++jface) {

            // get neighbor element face nodes
            std::vector< int > nfn( pillaz_getElmNNodes(pillaz_getElmFaceType(getElmType(ne->first,ne->second),jface)), -1 );
            pillaz_getElmFaceNodes(ne->first,ne->second,jface,&nfn[0]);
            std::sort(nfn.begin(),nfn.end());

            // if neighbor face matches, it must be that one
            if (std::equal(fn.begin(),fn.end(),nfn.begin())) {
              m_elemtoelem[iz][ie][iface] = pillaz_elementaddress(ne->first,ne->second);
              isneighborfacefound = true;
            }

          }
        }

        // if neighbor face is not found, it just might touch a boundary :)
        if (!isneighborfacefound)
          list_boundelems.insert(list_boundelems.end(),pillaz_elementaddress(iz,ie));

      }
    }
  }
  screenOutput("recreating map of (inner) elements to elements (sharing a face).");


  screenOutput("recreating map of boundary elements to (inner) elements...");

  // allocate
  count = 0;
  m_boundtoinner.resize(fp.nBoundElements.size());
  for (size_t iz=0; iz<fp.nBoundElements.size(); ++iz) {
    count += fp.nBoundElements[iz];
    if (fp.nBoundElements[iz])
      m_boundtoinner[iz].resize(fp.nBoundElements[iz]);
  }

  // perform the boundary elements to (inner) elements search
  pbar.restart(count);
  for (size_t iz=0; iz<m_boundtoinner.size(); ++iz) {
    for (size_t ie=0; ie<m_boundtoinner[iz].size(); ++ie, ++pbar) {

      // get boundary element nodes (sorted for comparison later)
      std::vector< int > ben( pillaz_getElmNNodes(getBndElmType(iz,ie)), -1 );
      getBndElmNodes(iz,ie,&ben[0]);
      std::sort(ben.begin(),ben.end());

      // find the element with a face with the same nodes
      bool isinnerfacefound = false;
      for (std::list< pillaz_elementaddress >::iterator be=list_boundelems.begin(); be!=list_boundelems.end() && !isinnerfacefound; ++be) {
        for (int jface=0; jface<pillaz_getElmNFaces(getElmType(be->izone,be->ielem)) && !isinnerfacefound; ++jface) {

          // get inner element face nodes
          std::vector< int > befn( pillaz_getElmNNodes(pillaz_getElmFaceType(getElmType(be->izone,be->ielem),jface)), -1 );
          pillaz_getElmFaceNodes(be->izone,be->ielem,jface,&befn[0]);
          std::sort(befn.begin(),befn.end());

          // if neighbor face matches, it must be that one
          if (std::equal(ben.begin(),ben.end(),befn.begin())) {
            m_boundtoinner[iz][ie]       = *be;
            m_boundtoinner[iz][ie].iface = jface;
            list_boundelems.erase(be);
            isinnerfacefound = true;
          }

        }
      }

    }
  }
  screenOutput("recreating map of boundary elements to (inner) elements.");


  screenOutput("calculating inner element normals...");

  // allocate
  count = 0;
  m_innerelem_normals.resize(fp.nInnerElements.size());
  for (size_t iz=0; iz<fp.nInnerElements.size(); ++iz) {
    if (fp.nInnerElements[iz])
      m_innerelem_normals[iz].resize(fp.nInnerElements[iz]);
    for (size_t ie=0; ie<fp.nInnerElements[iz]; ++ie) {
      m_innerelem_normals[iz][ie].assign( pillaz_getElmNFaces(getElmType(iz,ie)), std::vector< double >(fp.numDim,0.) );
      count += m_innerelem_normals[iz][ie].size();
    }
  }

  // calculate zone/element/face normals
  pbar.restart(count);
  for (size_t iz=0; iz<m_innerelem_normals.size(); ++iz) {
    for (size_t ie=0; ie<m_innerelem_normals[iz].size(); ++ie)
      for (size_t f=0; f<m_innerelem_normals[iz][ie].size(); ++f, ++pbar)
        pillaz_CalcElmFaceNormal( iz,ie,f, &(m_innerelem_normals[iz][ie][f])[0] );
  }
  screenOutput("calculating inner element normals.");


  screenOutput("calculating element volumes...");
  volumeElm.resize(fp.nInnerElements.size());
  for (size_t iz=0; iz<fp.nInnerElements.size(); ++iz) {
    if (fp.nInnerElements[iz])
      volumeElm[iz].assign(fp.nInnerElements[iz],0.);
    for (size_t ie=0; ie<(size_t) fp.nInnerElements[iz]; ++ie)
      volumeElm[iz][ie] = pillaz_CalcElmSize(iz,ie);
  }
  screenOutput("calculating element volumes.");


  screenOutput("calculating nodal volumes (dual cell)...");
  volumeNod.assign(fp.numNod,0.);
  for (int in=0; in<fp.numNod; ++in)
    for (std::vector< pillaz_elementaddress >::const_iterator e=m_nodetoelem[in].begin(); e!=m_nodetoelem[in].end(); ++e)
      volumeNod[in] += volumeElm[e->izone][e->ielem] / double(pillaz_getElmNNodes(getElmType(e->izone,e->ielem)));
  screenOutput("calculating nodal volumes (dual cell).");


  // allocate memory for dispersed phase data
  ed = new PILLAZ_ENTITY_DATA[ip.numMaxEnt];
  for (int ient=0; ient<ip.numMaxEnt; ++ient) {
    ed[ient].flag        = DFLAG_DISABLED;
    ed[ient].position    = new double[fp.numDim];
    ed[ient].velocity    = new double[fp.numDim];
    ed[ient].velocityOld = new double[fp.numDim];
    ed[ient].node        = -1;
  }

  pd = new PILLAZ_PHASE_DATA[fp.numNod];
  for (int inod=0; inod<fp.numNod; ++inod) {
    pd[inod].volFracDx = new double[fp.numDim];
    pd[inod].avgVel    = new double[fp.numDim];
    pd[inod].stdVel    = new double[fp.numDim];
    pd[inod].dispForce = new double[fp.numUnk];
  }


  //***Initialize variable runtime parameters***//
  if (ip.numProdDom>0)
    massResid = new double[ip.numProdDom];


  //***Initialize material data from database***//
  mdd->update(ip.iniTempDisp,getQuantityScalar(PRESSURE,0));


  // set output tecplot file (and initial entity distribution)
  bool TecplotFileExists( std::ifstream(ip.writeTecplotFilename.c_str()) );
  TecplotFileExists = false;   //FIXME!!!
  if (ip.writeTecplotFilename.length() && !TecplotFileExists)
    pillaz_CreateTecplotFile(ip.writeTecplotFilename);

  if (TecplotFileExists)
    pillaz_LoadInitialDistribution(ip.writeTecplotFilename);
  else if (ip.numIniEnt>0)
    pillaz_RandomInitialDistribution();


  // set output statistics file
  if (ip.writeStatsFilename.length())
    pillaz_CreateStatsFile(ip.writeStatsFilename);


  screenOutput("initializing Pillaz.");
}


void pillaz::run()
{
  screenOutput("perform Pillaz iteration...");


  LOCAL_ENTITY_VARIABLES ent(fp.numDim);
  LOCAL_FLOW_VARIABLES flow(fp.numDim);
  char errMessage[100];
  int idx,ient,jent,inod,idim,iunk,ibnd,ifac,subIter,avgctr;
  double cellSize,minCellSize,entVel,dtLagr,dtRemaining;


  //***Flow solver parameters that have to be set at every time step***//
  fp.time += fp.dtEul;
  ++fp.iter;


  //***Initializations***//
  sd.enabled     = 0;
  sd.in          = 0;
  sd.out         = 0;
  sd.bounce      = 0;
  sd.lost        = 0;
  sd.coll        = 0;
  sd.coalesc     = 0;
  sd.dtLagrAvg   = 0.;
  sd.reynoldsAvg = 0.;
  sd.nusseltAvg  = 0.;
  sd.subIterAvg  = 0.;
  avgctr         = 0;

  // find minimum cell size
  {
    double min_volumeElm = 1.e99;
    for (size_t iz=0; iz<volumeElm.size(); ++iz)
      for  (size_t ie=0; ie<volumeElm[iz].size(); ++ie)
        min_volumeElm = std::min(min_volumeElm,volumeElm[iz][ie]);
    minCellSize = pow(min_volumeElm,(1./fp.numDim));
  }

  if(fp.iter==1){
    pillaz_CalcCellwiseData();
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
      ed[ient].eaddress = ed[jent].eaddress;
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
    screenOutput("imposing new entities in production domains...");
    pillaz_ImposeProductionDomains();
  }

  if (numExtEnt>0){
    screenOutput("imposing externally generated entities...");
    pillaz_ImposeExternal();
  }


  //***Loop over dispersed entities to update trajectories***//
  screenOutput("updating trajectories...");
  for(ient=0; ient<sd.enabled; ient++){


    //***Get entity information from global data structure***//
    ent.flag = ed[ient].flag;
    if (ent.flag!=DFLAG_ENABLED && ent.flag!=DFLAG_CREATED) {
      sprintf(errMessage,"entity %i had bad flag %d.",ient,ent.flag);
      pillaz_TerminateOnError(errMessage);
    }

    ent.eaddress = ed[ient].eaddress;
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
      entVel   = std::max(ip.errTol, pillaz_CalcVectorLength(fp.numDim,&ent.vel[0]) );
      cellSize = std::max(ip.errTol, pow( volumeElm[ent.eaddress.izone][ent.eaddress.ielem], 1./fp.numDim) );
      dtLagr   = std::min(dtRemaining, ip.lagrTimeFactor*(cellSize/entVel) );
      dtRemaining -= dtLagr;
      subIter++;


      //***Entity production***//
      if ((ent.flag==DFLAG_CREATED) && (pillaz_RandomDouble()<=((fp.dtEul-dtRemaining)/fp.dtEul)))
        ent.flag = DFLAG_ENABLED;

      if(ent.flag!=DFLAG_ENABLED){
        continue;
      }


      //***Interpolate continuous phase velocity and pressure at entity position***//
      pillaz_SetElementGeometry(fp.numDim,&ent);
      pillaz_Interpolate(&ent,&flow,(fp.dtEul-dtRemaining)/fp.dtEul);
      pillaz_CalcEntityCoefficients(&ent,&flow);


      //***Solve Lagrangian trajectory equation***//
      pillaz_CalcTrajectory(&ent,&flow,dtLagr);
      pillaz_CheckNaN(&ent);


      //***Check for evaporated, burned or collapsed entity***//
      if (ent.diam<ip.errTol) {
        ent.flag = DFLAG_DISABLED;
      }


      //***Re-calculate material data from database***//
      if(ip.energyCoupl){
        mdd->update(ent.temp,flow.pressure);
      }

      if(ent.flag!=DFLAG_DISABLED){


        //***Perform neighbour element search routine***//
        pillaz_SearchSuccessive(&ent);


        // treat entities that left the domain after unsuccessful element search
        if (ent.flag==DFLAG_LEFT) {
          if (pillaz_FindExitFace(&ent,&ibnd,&ifac) && getWallBndFlag(ibnd)) {

            // perform wall bounce
            pillaz_WallBounce(fp.numDim,ip.wallElasticity,&ent,ibnd,ifac);
            sd.bounce++;
            ent.flag = DFLAG_ENABLED;

          }
          else {

            // entity left computational domain through outlet
            ent.flag = DFLAG_DISABLED;
            sd.out++;

          }
        }


        //***Stochastic collision model***//
        if(fp.numDim==3 && ent.flag==DFLAG_ENABLED && ip.collisionModel && pd[ent.node].volFrac>1e-3){
          pillaz_CollisionModel(&ent,pd[ent.node].numDens,dtLagr);
          pillaz_CheckNaN(&ent);
        }
      }

      if(ent.flag==DFLAG_ENABLED){


        //***Interpolate continuous phase velocity and pressure to new entity position***//
        pillaz_SetElementGeometry(fp.numDim,&ent);
        pillaz_Interpolate(&ent,&flow,(fp.dtEul-dtRemaining)/fp.dtEul);
        pillaz_CalcEntityCoefficients(&ent,&flow);


        //***Compute back-coupling terms to flow solver***//
        const pillaz_flowtype_t flowType(mdd->flowtype);
        if(flowType==FLOW_PARTIC || flowType==FLOW_DROPLET){
          pillaz_CalcCouplingForcesParticle(&ent,&flow,dtLagr/fp.dtEul);
        } else if (flowType==FLOW_BUBBLY){
          pillaz_CalcCouplingForcesBubble(&ent,&flow,dtLagr/fp.dtEul);
        }
      }

    } while(dtRemaining>0.0); /***End subiteration loop***/


    //***Re-calculate and update statistics***//
    if(ent.flag==DFLAG_ENABLED){
      pillaz_SetElementGeometry(fp.numDim,&ent);
      pillaz_Interpolate(&ent,&flow,1.0);
      pillaz_CalcEntityCoefficients(&ent,&flow);
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
      ed[ient].eaddress = ent.eaddress;
      ed[ient].node = ent.node;
    } else{
      ed[ient].eaddress.reset();
      ed[ient].node = -1;
    }

  } //***End trajectory loop***//



  screenOutput("post-processing...");


  //***Update statistics***//
  sd.enabled  -= sd.out;
  sd.enabled  -= sd.lost;
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
  pillaz_CalcCellwiseData();


  //***Write Pillaz output to files***//
  if (ip.writeStatsFilename.length())
    pillaz_WriteStatsFile(ip.writeStatsFilename,fp.iter,fp.time);
  if (ip.writeTecplotFilename.length())
    pillaz_WriteTecplotFile(ip.writeTecplotFilename,fp.iter,fp.time);


  screenOutput("perform Pillaz iteration.");
}


pillaz::~pillaz()
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
    delete[] massResid;
  }
}


void pillaz::pillaz_CalcBackCoupling(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double *force, double tFactor)
{
  double limitFrac = 0.95;
  int idim,jdim;
  double contVol,iFrac,jFrac,rhsTerm;

  std::vector< double > impFac(pillaz_getElmNNodes(getElmType(ent->eaddress.izone,ent->eaddress.ielem)),0.);

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

    contVol = volumeNod[ent->node];
    for(idim=0; idim<fp.numDim; idim++){
      pd[ent->node].dispForce[idim+1] -= tFactor*force[idim]*jFrac/((1.0-iFrac)*contVol * mdc->rho);
    }

  }
  else if(ip.momentumCoupl==FORCE_PROJ){

    pillaz_CalcNodeImpactFactors(ent,&impFac[0]);
    for(size_t idim=0; idim<ent->edata.elmNodes.size(); idim++){
      const int inod = ent->edata.elmNodes[idim];
      contVol = volumeNod[inod];
      for(jdim=0; jdim<fp.numDim; jdim++){
        pd[inod].dispForce[jdim+1] -= tFactor*impFac[idim]*force[jdim]*jFrac/((1.0-iFrac)*contVol * mdc->rho);
      }
    }
  }

  //***Compute back-coupling term for fluid continuity equation***//
  if(ip.volfracCoupl){

    const int inod = ent->node;
    rhsTerm = pd[inod].volFracDt;
    for(idim=0; idim<fp.numDim; idim++)
      rhsTerm += flow->vel[idim]*pd[inod].volFracDx[idim];

    pd[inod].dispForce[0] += tFactor*jFrac*rhsTerm/(1.0-iFrac);
    for(idim=1; idim<fp.numDim; idim++)
      pd[inod].dispForce[idim] += tFactor*jFrac*rhsTerm*flow->vel[idim]/(1.0-iFrac);

  }
}


void pillaz::pillaz_CalcBoundaryUnitNormal(int numDim, int ibnd, int ifac, double *unitVec)
{
  // get normal vector of a boundary face
  const pillaz_elementaddress &in = m_boundtoinner[ibnd][ifac];
  for (int idim=0; idim<numDim; ++idim)
    unitVec[idim] = m_innerelem_normals[in.izone][in.ielem][in.iface][idim];

  // divide normal vector components by vector length
  const double length = pillaz_CalcVectorLength(numDim,unitVec);
  for (int idim=0; idim<numDim; ++idim)
    unitVec[idim] /= length;
}


void pillaz::pillaz_CalcCellwiseData()
{
  //***Initialize cellwise secondary phase data***//
  double volFracOld[fp.numNod];
  for (int inod=0; inod<fp.numNod; ++inod){
    pd[inod].numDens = 0;
    volFracOld[inod] = pd[inod].volFrac;
    pd[inod].volFrac = 0.0;
    pd[inod].volFracDt = 0.0;
    pd[inod].avgDiam = 0.0;
    pd[inod].stdDiam = 0.0;
    pd[inod].avgRespTime = 0.0;
    for (int idim=0; idim<fp.numDim; ++idim) {
      pd[inod].volFracDx[idim] = 0.;
      pd[inod].avgVel[idim]    = 0.;
      pd[inod].stdVel[idim]    = 0.;
    }
  }


  //***Assemble cellwise data by looping over all active entities***//
  for (int ient=0; ient<ip.numMaxEnt; ++ient) {
    if(ed[ient].flag==DFLAG_ENABLED){
      const int inod = ed[ient].node;
      pd[inod].numDens++;

      const double ivol = (fp.numDim==2? PI*ed[ient].diameter*ed[ient].diameter/4. :
                          (fp.numDim==3? PI*ed[ient].diameter*ed[ient].diameter*ed[ient].diameter/6. :
                                         0. ));
      pd[inod].volFrac += ivol/volumeNod[inod];
      pd[inod].avgDiam += ed[ient].diameter;
      pd[inod].avgRespTime += (2.0 * mdd->rho + mdc->rho)
                                   *ed[ient].diameter*ed[ient].diameter/(24.*mdc->mu);
      for (int idim=0; idim<fp.numDim; ++idim) {
        pd[inod].avgVel[idim] += ed[ient].velocity[idim];
      }
    }
  }


  //***Divide averaged quantities by number density***//
  for (int inod=0; inod<fp.numNod; inod++){
    if(pd[inod].numDens>0){
      pd[inod].avgDiam /= pd[inod].numDens;
      pd[inod].avgRespTime /= pd[inod].numDens;
      for (int idim=0; idim<fp.numDim; ++idim) {
        pd[inod].avgVel[idim] /= pd[inod].numDens;
      }
    }
  }


  //***Compute standard deviations of averaged quantities***//
  if(ip.collisionModel){

    for (int ient=0; ient<ip.numMaxEnt; ++ient) {
      if(ed[ient].flag==DFLAG_ENABLED){
        const int inod = ed[ient].node;

        pd[inod].stdDiam += pow(ed[ient].diameter-pd[inod].avgDiam,2.0);
        for (int idim=0; idim<fp.numDim; ++idim) {
          pd[inod].stdVel[idim] += pow(ed[ient].velocity[idim]-pd[inod].avgVel[idim],2.0);
        }
      }
    }

    for (int inod=0; inod<fp.numNod; inod++){
      if(pd[inod].numDens>0){
        pd[inod].stdDiam = sqrt(pd[inod].stdDiam/pd[inod].numDens);
        for (int idim=0; idim<fp.numDim; ++idim) {
          pd[inod].stdVel[idim] = sqrt(pd[inod].stdVel[idim]/pd[inod].numDens);
        }
      }
    }
  }


  //***Calculate volume fraction derivatives in time and space***//
  if (ip.volfracCoupl) {

    for (int n=0; n<fp.numNod; ++n)
      pd[n].volFracDt = (pd[n].volFrac - volFracOld[n])/fp.dtEul;

    for (size_t iz=0; iz<fp.nInnerElements.size(); ++iz) {
      for (size_t ie=0; ie<fp.nInnerElements[iz]; ++ie) {

        // get element nodes and coordinates
        const int     NElmNodes = pillaz_getElmNNodes(getElmType(iz,ie));
        const double dNElmNodes = (double) NElmNodes;
        std::vector< int                   > en(NElmNodes,-1);
        std::vector< std::vector< double > > ec(NElmNodes,std::vector< double >(fp.numDim,0.));

        getElmNodes(iz,ie,&en[0]);
        for (int in=0; in<NElmNodes; ++in)
          getQuantityVector(COORD,en[in],fp.numDim,&(ec[in])[0]);

        for (int idim=0; idim<fp.numDim; ++idim) {
          double
            xi = 0.,
            yi = 0.,
            xixi = 0.,
            yixi = 0.;
          for (int in=0; in<NElmNodes; ++in) {
            xi   += ec[in][idim];
            xixi += ec[in][idim]*ec[in][idim];
            yi   += pd[ en[in] ].volFrac;
            yixi += pd[ en[in] ].volFrac * ec[in][idim];
            pd[ en[in] ].volFracDx[idim] = (dNElmNodes * yixi - xi * yi)
                                         / (dNElmNodes * xixi - xi * xi);
          }
        }

      }
    }

  }
}


double pillaz::pillaz_CalcConcInterf(double pressBubble)
{
  return (pressBubble*mdd->He);
}


void pillaz::pillaz_CalcCouplingForcesBubble(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor)
{
  int idim,jdim;
  double densityRatio,bubbleVol,bubbleMass,vmCoeff;
  double uxw[fp.numDim],dudt[fp.numDim],dvdt[fp.numDim],entForce[fp.numDim];

  //***Compute momentum coupling forces***//

  if(ip.momentumCoupl){

    bubbleVol = PI*ent->diam*ent->diam*ent->diam/6.0;
    bubbleMass = mdd->rho * bubbleVol;
    densityRatio = (mdc->rho / mdd->rho);

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

  pillaz_CalcBackCoupling(ent,flow,entForce,tFactor);

}


void pillaz::pillaz_CalcCouplingForcesParticle(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor)
{
  int idim;
  double densityRatio,particleVol,particleMass;
  double uxw[fp.numDim],entForce[fp.numDim];

  //***Compute momentum coupling forces***//

  if(ip.momentumCoupl){

    particleVol = PI*ent->diam*ent->diam*ent->diam/6.0;
    particleMass = mdd->rho * particleVol;
    densityRatio = (mdc->rho / mdd->rho);

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

  pillaz_CalcBackCoupling(ent,flow,entForce,tFactor);

}


void pillaz::pillaz_CalcCrossProduct_3D(double *value, double *a, double *b)
{
  value[0] = a[1]*b[2] - a[2]*b[1];
  value[1] = a[2]*b[0] - a[0]*b[2];
  value[2] = a[0]*b[1] - a[1]*b[0];
}


double pillaz::pillaz_CalcDispReynolds(double viscosity, double diameter, double normVel)
{
  double reynolds = diameter*normVel/viscosity;

  if(reynolds<1e-4){reynolds = 1e-4;}

  return reynolds;
}


double pillaz::pillaz_CalcDragCoeff(int flowType, double reynolds)
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


void pillaz::pillaz_CalcEntityCoefficients(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow)
{
  //***Compute flow coefficients***//

  ent->reynolds      = pillaz_CalcDispReynolds(mdc->mu / mdc->rho,ent->diam,ent->normVel);
  ent->dragCoeff     = pillaz_CalcDragCoeff(mdd->flowtype,ent->reynolds);
  ent->liftCoeff     = pillaz_CalcLiftCoeff(mdd->flowtype);
  ent->kinRespTime   = pillaz_CalcKinematicResponseTime(ent);
  ent->spalding      = pillaz_CalcSpaldingNumber(flow->pressure);
  ent->prandtl       = pillaz_CalcPrandtlNumber();
  ent->nusselt       = pillaz_CalcNusseltNumber(ip.evapModel,ent->reynolds,ent->spalding,ent->prandtl);
  ent->thermRespTime = pillaz_CalcThermalResponseTime(ent->diam);
  ent->schmidt       = pillaz_CalcSchmidtNumber();
  ent->sherwood      = pillaz_CalcSherwoodNumber(ip.evapModel,ent->reynolds,ent->schmidt,ent->spalding);
  ent->massTrCoeff   = pillaz_CalcMassTransferCoeff(ent->sherwood,ent->spalding);
  ent->pressBubble   = pillaz_CalcPressBubble(ent->diam,flow->pressure);
  ent->rhoBubble     = pillaz_CalcRhoBubble(ent->temp, ent->pressBubble);
  ent->concInterf    = pillaz_CalcConcInterf( ent->pressBubble);
}


double pillaz::pillaz_CalcKinematicResponseTime(LOCAL_ENTITY_VARIABLES *ent)
{
  double tau = 0.;

  const pillaz_flowtype_t flowType(mdd->flowtype);
  if(flowType==FLOW_PARTIC || flowType==FLOW_DROPLET){
    tau = 4.0 * mdd->rho * ent->diam*ent->diam/(3.*mdc->mu*ent->reynolds*ent->dragCoeff);
  } else if(flowType==FLOW_BUBBLY){
    tau = 2.0*ent->diam*ent->diam/(3. * mdc->mu / mdc->rho * ent->reynolds * ent->dragCoeff);
  }

  return tau;
}


double pillaz::pillaz_CalcLiftCoeff(int flowType)
{
  if (flowType==FLOW_PARTIC || flowType==FLOW_DROPLET)
    return 0.;
  else if (flowType==FLOW_BUBBLY)
    return 0.53;
  return 0.;
}


double pillaz::pillaz_CalcMassTransferCoeff(double sherwood, double spalding)
{
  double omega = log(1.0+spalding);

  return (2.0*(mdc->rho / mdd->rho) * mdd->binaryDiffCoeff * sherwood*omega);
}


void pillaz::pillaz_CalcMatVectScalarProduct_2D(double *value, double **m, double *a)
{
  value[0] = m[0][0]*a[0]+m[0][1]*a[1];
  value[1] = m[1][1]*a[0]+m[1][1]*a[1];
}


void pillaz::pillaz_CalcMatVectScalarProduct_3D(double *value, double **m, double *a)
{
  value[0] = m[0][0]*a[0]+m[0][1]*a[1]+m[0][2]*a[2];
  value[1] = m[1][0]*a[0]+m[1][1]*a[1]+m[1][2]*a[2];
  value[2] = m[2][0]*a[0]+m[2][1]*a[1]+m[2][2]*a[2];
}


void pillaz::pillaz_CalcNodeImpactFactors(const LOCAL_ENTITY_VARIABLES *ent, double *imp)
{
  double sum = 0.;
  double distance[fp.numDim];

  for (unsigned in=0; in<ent->edata.elmNodes.size(); ++in) {
    getQuantityVector(COORD,ent->edata.elmNodes[in],fp.numDim,distance);
    imp[in] = pillaz_CalcVectorLength(fp.numDim,distance);
    if (imp[in]<ip.errTol)
      imp[in] = ip.errTol;
    sum += 1.0/imp[in];
  }

  for (unsigned in=0; in<ent->edata.elmNodes.size(); in++)
    imp[in] = 1./(imp[in]*sum);
}


double pillaz::pillaz_CalcNusseltNumber(int evapModel, double reynolds, double spalding, double prandtl)
{
 double Nu = 2.0 + 0.6*sqrt(reynolds)*pow(prandtl,(1.0/3.0));

 if(evapModel){
   Nu *= (log(1.0+spalding)/spalding);
 }

 return Nu;
}


double pillaz::pillaz_CalcPrandtlNumber()
{
  return (mdc->mu * mdc->cp / mdc->k);
}


double pillaz::pillaz_CalcPressBubble(double diameter, double pressure)
{
  return (pressure - mdd->satPres + 4. * mdd->sig / diameter);
}


double pillaz::pillaz_CalcRhoBubble(double temperature, double pressBubble)
{
  return (mdd->molarMass / (Ru*temperature)*pressBubble);
}


void pillaz::pillaz_CalcRotationMatrix_2D(double phi, double **m)
{
  m[0][0] = cos(phi);
  m[0][1] = -sin(phi);
  m[1][0] = sin(phi);
  m[1][1] = cos(phi);
}


void pillaz::pillaz_CalcRotationMatrix_3D(double phi, double **m, int axis)
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


double pillaz::pillaz_CalcSchmidtNumber()
{
  return (( mdd->mu * mdc->k)/( mdd->rho * mdc->rho * mdc->cp));
}


double pillaz::pillaz_CalcSherwoodNumber(int evapModel, double reynolds, double schmidt, double spalding)
{
  double Sh = 2.0 + 0.6*pow(reynolds,(1.0/2.0))*pow(schmidt,(1.0/3.0));

  if(evapModel){
    double F_M = pow(1.0+spalding,0.7)*log(1.0+spalding)/spalding;
    Sh = 2.0+(Sh-2.0)/F_M;
  }

  return Sh;
}


double pillaz::pillaz_CalcSpaldingNumber(double pressure)
{
  if ( mdd->satPres <1.e-20 ||  mdd->molarMass <1.e-20)
    return 0.;

  double Y_s = 1.0/(1.0+(pressure/ mdd->satPres -1.0)*(mdd->molarMassVap / mdd->molarMass ));
  return ((Y_s)/(1.0-Y_s));
}


double pillaz::pillaz_CalcThermalResponseTime(double diameter)
{
  return (1.0/(12. * mdc->k))*( mdd->rho *diameter*diameter* mdd->cp );
}


void pillaz::pillaz_CalcTrajectory(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double dtLagr)
{
  const pillaz_flowtype_t flowType(mdd->flowtype);

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
  double flow_concentration = flow->pressure * mdd->He * alpha;

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

  if(ip.liftForce && flowType==FLOW_BUBBLY && fp.numDim==3){
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

  if(flowType==FLOW_BUBBLY){
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
    if(flowType==FLOW_PARTIC || flowType==FLOW_DROPLET){
      rhsvec[idim] += ip.gravVec[idim];
    } else if(flowType==FLOW_BUBBLY){
      rhsvec[idim] -= 2.0*ip.gravVec[idim];
    }
  }

  //***Temperature equation contribution to matrix and RHS***//

  if(ip.energyCoupl && (flowType==FLOW_PARTIC || flowType==FLOW_DROPLET)){
    convTerm = ent->nusselt/(2.0*ent->thermRespTime);
    radTerm = (6.0*boltzmann* mdd->eps )/( mdd->rho * mdd->cp *ent->diam);
    if(ip.evapModel && flowType==FLOW_DROPLET){
      massTerm = (3.0)*(ent->massTrCoeff/pow(ent->diam,2.0))*(mdd->latHeat / mdd->cp );
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

  if(ip.evapModel && flowType==FLOW_DROPLET ){
    mat[fp.numDim+1][fp.numDim+1] -= theta*(ent->massTrCoeff/pow(ent->diam,2.0))*(1.5-1.0/ent->sherwood);
    rhsvec[fp.numDim+1] += -ent->massTrCoeff/ent->diam;
  }

  //***Concentration Boundary Layer Model Payvar for bubble growth and Epstein & Plesset for Bubble shrink

  if(ip.saturModel && flowType==FLOW_BUBBLY ){
   concTerm  = 4.0 * mdd->massDiffCoeff * mdc->rho * (flow_concentration-ent->concInterf)/(ent->rhoBubble*(mdc->rho - ent->concInterf));

    if(flow_concentration > ent->concInterf){
      rhsvec[fp.numDim+1] += concTerm/ent->diam*(1.0+1.0/(pow(1.0+( flow_concentration - mdc->rho )/(ent->concInterf-flow_concentration)*( ent->rhoBubble / mdc->rho )*(1.0-( mdd->rho /ent->rhoBubble)*pow(2.0e-5/ent->diam,3.0)),0.5)-1.0));
      mat[fp.numDim+1][fp.numDim+1] -= concTerm*((1.0/(ent->diam+eps)+1.0/((ent->diam+eps)*pow(1.0+(flow_concentration - mdc->rho)/(ent->concInterf-flow_concentration)*(ent->rhoBubble / mdc->rho)*(1.0-( mdd->rho /ent->rhoBubble)*pow(2.0e-5/(ent->diam+eps),3.0)),0.5)-1.0))-(1.0/ent->diam+1.0/(ent->diam*pow(1.0+(flow_concentration - mdc->rho)/(ent->concInterf-flow_concentration)*(ent->rhoBubble / mdc->rho)*(1.0-( mdd->rho /ent->rhoBubble)*pow(2.0e-5/ent->diam,3.0)),0.5)-1.0)))/eps;

      }else{
      rhsvec[fp.numDim+1] += concTerm/ent->diam;
      mat[fp.numDim+1][fp.numDim+1] -= ((concTerm/(ent->diam+eps))-(concTerm/ent->diam))/eps;
        }
  }

  //***Solve linear system for velocity, temperature & diameter***//

  pillaz_SolveGaussSeidel(fp.numDim+2,mat,solvec,rhsvec);

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


double pillaz::pillaz_CalcVectorAngle(int numDim, double *a, double *b)
{
  return acos(pillaz_CalcVectorScalarProduct(numDim,a,b)/(pillaz_CalcVectorLength(numDim,a)*pillaz_CalcVectorLength(numDim,b)));
}


double pillaz::pillaz_CalcVectorLength(int numDim, double *a)
{
  return sqrt(pillaz_CalcVectorScalarProduct(numDim,a,a));
}


void pillaz::pillaz_CalcVectorRotation_3D(double phi, double **m, double *a)
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


double pillaz::pillaz_CalcVectorScalarProduct(int numDim, double *a, double *b)
{
  if (numDim==2)
    return a[0]*b[0]+a[1]*b[1];
  else if (numDim==3)
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return 0.;
}


void pillaz::pillaz_CalcVorticity(int numDim, LOCAL_FLOW_VARIABLES *flow)
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


double pillaz::pillaz_CalcWallFaceDistance(int numDim, double *pos, int ibnd, int ifac)
{
  double posVec[numDim],unitVec[numDim];

  // calculate unit normal of boundary segment
  pillaz_CalcBoundaryUnitNormal(numDim,ibnd,ifac,unitVec);

  // calculate and return wall face distance
  std::vector< int > fn(pillaz_getElmNNodes(getElmType(ibnd,ifac)),-1);
  getBndElmNodes(ibnd,ifac,&fn[0]);

  // use first element node as reference (could be done better)
  getQuantityVector(COORD,fn[0],numDim,posVec);

  for (int d=0; d<numDim; ++d)
    posVec[d] = pos[d] - posVec[d];

  return pillaz_CalcVectorScalarProduct(numDim,posVec,unitVec);
}


void pillaz::pillaz_CheckNaN(LOCAL_ENTITY_VARIABLES *ent)
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
    screenWarning("not-a-number detected");
  }
}


int pillaz::pillaz_Coalescence(double dj, double di, double *uijRelPrPr, double *x, double *y, double *z,
  double Mi, double Mj, double *uiPrPr, double *uiPrPrNew, double *uiNew, double *ui)
{
  double h0 = 1e-4;
  double hf = 1e-8;
  double Cc = 0.5;
  double sigma = 0.06;
  double Rij = 2.0*1/(2/dj+2/di);
  double T = sqrt(Rij*Rij*Rij * mdc->rho/(16*sigma))*log(h0/hf);
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

      uiNew[0] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,x);
      uiNew[1] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,y);

    } else if(fp.numDim==3){

      uiPrPrNew[0] = uiPrPr[0]+Mi/(Mi+Mj);
      uiPrPrNew[1] = uiPrPr[1];
      uiPrPrNew[2] = uiPrPr[2];

      uiNew[0] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,x);
      uiNew[1] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,y);
      uiNew[2] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,z);
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


void pillaz::pillaz_CollisionModel(LOCAL_ENTITY_VARIABLES *ent, int numDens, double dtLagr)
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
    uj[idim] = pillaz_RandomGaussian(pd[ent->node].avgVel[idim],pd[ent->node].stdVel[idim]);
  }
  dj = pillaz_RandomGaussian(pd[ent->node].avgDiam,pd[ent->node].stdDiam);

  //***Correlation between the fluctuations of the real and fictitious particle***//

  if(ip.collisionModel==2){

    for(idim=0; idim<fp.numDim; idim++){
      uiFluct[idim] = ui[idim]-pd[ent->node].avgVel[idim];
    }

    tau_f = c_T*getEulerianTimeScale(ent->node);
    stokes = ent->kinRespTime/tau_f;
    epsilon = pillaz_RandomGaussian(0.0,1.0);
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

  conc = numDens/volumeNod[ent->node];
  FreqColl = PI/4.0*(di+dj)*(di+dj)*pillaz_CalcVectorLength(fp.numDim,uijRel)*conc;
  Pcoll = FreqColl*dtLagr;

  //***Perform collision***//

  if(pillaz_RandomDouble()<Pcoll){

    //***Transform into 2nd coordinate system***//

    x[0] = y[1] = z[2] = 1.0;
    x[1] = x[2] = y[0] = y[2] = z[0] = z[1] = 0.0;

    for(idim=0; idim<fp.numDim; idim++){
      xPr[idim] = uijRel[idim]/pillaz_CalcVectorLength(fp.numDim,uijRel);
    }

    pillaz_CalcCrossProduct_3D(yPr,xPr,x);
    pillaz_CalcCrossProduct_3D(zPr,xPr,yPr);

    pos1Pr[0] = pillaz_CalcVectorScalarProduct(fp.numDim,&ent->pos[0],xPr);
    pos1Pr[1] = pillaz_CalcVectorScalarProduct(fp.numDim,&ent->pos[0],yPr);
    pos1Pr[2] = pillaz_CalcVectorScalarProduct(fp.numDim,&ent->pos[0],zPr);

    while(L>1){
      Y = pillaz_RandomDouble();
      Z = pillaz_RandomDouble();
      L = sqrt(Y*Y+Z*Z);
    }

    Phi = asin(L);
    Psi = atan(Y/Z);

    pillaz_CalcRotationMatrix_3D(Phi,mat,2);
    pillaz_CalcMatVectScalarProduct_3D(pos2Pr,mat,xPr);
    pillaz_CalcVectorRotation_3D(Psi,mat,xPr);
    pillaz_CalcMatVectScalarProduct_3D(pos2Pr,mat,pos2Pr);

    for(idim=0; idim<fp.numDim; idim++){
      pos2Pr[idim] *= (dj+di)/2.0;
    }

    pos2[0] = pillaz_CalcVectorScalarProduct(fp.numDim,pos2Pr,x);
    pos2[1] = pillaz_CalcVectorScalarProduct(fp.numDim,pos2Pr,y);
    pos2[2] = pillaz_CalcVectorScalarProduct(fp.numDim,pos2Pr,z);

    //***Transform into 3rd coordinate system***//

    for(idim=0; idim<fp.numDim; idim++){
      posRel[idim] = pos1[idim]-pos2[idim];
    }

    for(idim=0; idim<fp.numDim; idim++){
      xPrPr[idim] = (posRel[idim])/pillaz_CalcVectorLength(fp.numDim,posRel);
    }

    pillaz_CalcCrossProduct_3D(zPrPr,xPrPr,uijRel);

    for(idim=0; idim<fp.numDim; idim++){
      zPrPr[idim] /= pillaz_CalcVectorLength(fp.numDim,zPrPr);
    }

    pillaz_CalcCrossProduct_3D(yPrPr,zPrPr,xPrPr);

    uiPrPr[0] = pillaz_CalcVectorScalarProduct(fp.numDim,ui,xPrPr);
    uiPrPr[1] = pillaz_CalcVectorScalarProduct(fp.numDim,ui,yPrPr);
    uiPrPr[2] = pillaz_CalcVectorScalarProduct(fp.numDim,ui,zPrPr);
    ujPrPr[0] = pillaz_CalcVectorScalarProduct(fp.numDim,uj,xPrPr);
    ujPrPr[1] = pillaz_CalcVectorScalarProduct(fp.numDim,uj,yPrPr);
    ujPrPr[2] = pillaz_CalcVectorScalarProduct(fp.numDim,uj,zPrPr);
    uijRelPrPr[0] = pillaz_CalcVectorScalarProduct(fp.numDim,uijRel,xPrPr);
    uijRelPrPr[1] = pillaz_CalcVectorScalarProduct(fp.numDim,uijRel,yPrPr);
    uijRelPrPr[2] = pillaz_CalcVectorScalarProduct(fp.numDim,uijRel,zPrPr);

    //***Sliding or non-sliding collision***//

    Mi =  mdd->rho *PI*(dj/2.0)*(dj/2.0);
    Mj =  mdd->rho *PI*(dj/2.0)*(dj/2.0);
    Jx = -(1.0-e)*uiPrPr[1]*(Mi*Mj)/(Mi+Mj);
    length = pillaz_CalcVectorLength(fp.numDim,uijRel);

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

    uiNew[0] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,x);
    uiNew[1] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,y);
    uiNew[2] = pillaz_CalcVectorScalarProduct(fp.numDim,uiPrPrNew,z);

    sd.coll++;

    for(idim=0; idim<fp.numDim; idim++){
      ent->vel[idim] = uiNew[idim];
    }
  }

  for (int r=0; r<fp.numDim; ++r)
    delete[] mat[r];
  delete[] mat;
}


void pillaz::pillaz_CreateStatsFile(const std::string &outpString)
{
  std::ofstream f(outpString.c_str());
  f << '\t' << "iter"
    << '\t' << "time"
    << '\t' << "ent"
    << '\t' << "in"
    << '\t' << "out"
    << '\t' << "bounce"
    << '\t' << "coll"
    << '\t' << "lost"
    << '\t' << "reynolds"
    << '\t' << "nusselt"
    << '\t' << "dtLagr"
    << '\t' << "subiter"
    << std::endl;
}


void pillaz::pillaz_CreateTecplotFile(const std::string &outpString)
{
  std::ofstream f(outpString.c_str());
  f << "VARIABLES ="
    << (fp.numDim>2? " \"X\" \"Y\" \"Z\"":" \"X\" \"Y\"")
    << (fp.numDim>2? " \"U\" \"V\" \"W\"":" \"U\" \"V\"")
    << " \"d\" \"T\" \"theta\""
    << std::endl;
}


bool pillaz::pillaz_FindExitFace(LOCAL_ENTITY_VARIABLES *ent, int *i, int *j)
{
  // loop over boundaries and faces to find the exit face
  for (size_t ibnd=0; ibnd<fp.nBoundElements.size(); ++ibnd) {
    for (size_t ifac=0; ifac<fp.nBoundElements[ibnd]; ++ifac) {
      const pillaz_elementaddress& beaddress = m_boundtoinner[ibnd][ifac];
      if (ent->eaddress == beaddress &&
          pillaz_CalcWallFaceDistance(fp.numDim,&ent->pos[0],ibnd,ifac)<0.) {
        *i = ibnd;
        *j = ifac;
        return true;
      }
    }
  }
  return false;
}


void pillaz::pillaz_FindMinimumElementFaceDistance(int numDim, LOCAL_ENTITY_VARIABLES *ent, int *idx, double *dmin)
{
  std::vector< double >
    posvec(numDim,0.),
    normvec(numDim,0.);
  *idx  = -1;
  *dmin = 1.e99;
  for (size_t eface=0; eface<ent->edata.elmFaceVectors.size(); ++eface) {
    for (int jdim=0; jdim<numDim; ++jdim) {
      posvec[jdim]  = ent->pos[jdim] - ent->edata.elmFaceVectors[eface][jdim];
      normvec[jdim] = m_innerelem_normals[ent->eaddress.izone][ent->eaddress.ielem][eface][jdim];
    }
    pillaz_NormalizeVector(numDim,&normvec[0]);
    const double idist = pillaz_CalcVectorScalarProduct(numDim,&posvec[0],&normvec[0]);
    if (idist<*dmin) {
      *idx  = eface;
      *dmin = idist;
    }
  }
}


int pillaz::pillaz_FindNearestElementNode(LOCAL_ENTITY_VARIABLES *ent)
{
  std::vector< double > impFac(pillaz_getElmNNodes(pillaz_getElmType(ent->eaddress)),0.);
  pillaz_CalcNodeImpactFactors(ent,&impFac[0]);

  int    fmax_nod = ent->edata.elmNodes[0];
  double fmax_val = impFac[0];
  for (size_t in=1; in<impFac.size(); ++in) {
    if (impFac[in]>fmax_val) {
      fmax_nod = ent->edata.elmNodes[in];
      fmax_val = impFac[in];
    }
  }

  return fmax_nod;
}


void pillaz::pillaz_ImposeExternal()
{
  LOCAL_ENTITY_VARIABLES ent(fp.numDim);
  LOCAL_FLOW_VARIABLES flow(fp.numDim);
  int ient,idim;
  double *newPos,*newVel,*newDiam,*newTemp;


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

    pillaz_SearchSuccessive(&ent);

    //***Initialize entity***//

    if(ent.flag==DFLAG_ENABLED && sd.enabled<ip.numMaxEnt){

      pillaz_SetElementGeometry(fp.numDim,&ent);
      pillaz_Interpolate(&ent,&flow,0.0);

      ed[sd.enabled].flag = DFLAG_CREATED;
      ed[sd.enabled].eaddress = ent.eaddress;
      ed[sd.enabled].node = pillaz_FindNearestElementNode(&ent);
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

  delete[] newDiam;
  delete[] newTemp;
  delete[] newPos;
  delete[] newVel;
}


void pillaz::pillaz_ImposeProductionDomains()
{
  int ient,ipd,idim,bCtr;
  double p,s,p1[3],p2[3];
  char msg[100];

  std::vector< double >
    newDiam(ip.numMaxEnt,0.),
    newPos(ip.numMaxEnt*fp.numDim);


  //***Creation of bubbles***//
  bCtr = 0;

  //***Loop over production domains***//

  for(ipd=0; ipd<ip.numProdDom; ipd++){

    if(ip.prodDom[ipd]==0){continue;}

    //***Accumulate produced secondary phase mass***//
    massResid[ipd] += ip.massFluxes[ipd]*fp.dtEul;
    if (massResid[ipd]<=0.)
      continue;

    //***Copy production parameters to local data structure***//
    for(idim=0; idim<3; idim++){
      p1[idim] = ip.prodParam[ipd][idim];
      p2[idim] = ip.prodParam[ipd][idim+3];
    }

    //***Iteratively generate entities according to mass flux***//
    for (; massResid[ipd]>0. && bCtr<ip.numMaxEnt; ++bCtr) {

      //***Set diameter***//
      newDiam[bCtr] = pillaz_SetDiameter();

      //***Calculate mass of generated entity***//
      const double mass =
        (fp.numDim==2?  mdd->rho *PI*newDiam[bCtr]*newDiam[bCtr]              /4. :
        (fp.numDim==3?  mdd->rho *PI*newDiam[bCtr]*newDiam[bCtr]*newDiam[bCtr]/6. :
                       0. ));

      massResid[ipd] -= mass;

      //***Compute position***//
      if(ip.prodDom[ipd]==1){

        //***Line production domain***//
        s = pillaz_RandomDouble();
        for(idim=0; idim<fp.numDim; idim++){
          newPos[fp.numDim*bCtr+idim] = p1[idim]+s*(p2[idim]-p1[idim]);
        }

      }
      else if(ip.prodDom[ipd]==2){

        //***Rectangle production domain***//
        for(idim=0; idim<fp.numDim; idim++){
          s = pillaz_RandomDouble();
          newPos[fp.numDim*bCtr+idim] = p1[idim]+s*(p2[idim]-p1[idim]);
        }

      }
      else if(ip.prodDom[ipd]==3){

        //***Ellipse production domain***//
        do{
          p = 0.0;
          for(idim=0; idim<fp.numDim; idim++){
            s = 2.0*pillaz_RandomDouble()-1.0;
            newPos[fp.numDim*bCtr+idim] = p1[idim]+s*p2[idim];
            if (p2[idim]>ip.errTol)
              p += (newPos[fp.numDim*bCtr+idim]-p1[idim])*(newPos[fp.numDim*bCtr+idim]-p1[idim])/(p2[idim]*p2[idim]);
          }
        } while(p>1.0);
      }
    }

  }


  //***Generate entities***//
  LOCAL_FLOW_VARIABLES flow(fp.numDim);
  LOCAL_ENTITY_VARIABLES ent(fp.numDim);
  for(ient=0; ient<bCtr; ient++){

    // find a starting (inner) element
    ent.eaddress.reset();
    for (size_t iz=0; iz<fp.nInnerElements.size() && !ent.eaddress.valid(); ++iz)
      if (fp.nInnerElements[iz])
        ent.eaddress = pillaz_elementaddress(iz,0);

    ent.diam = newDiam[ient];
    ent.temp = ip.iniTempDisp;
    for(idim=0; idim<fp.numDim; idim++){
      ent.pos[idim] = newPos[fp.numDim*ient+idim];
      ent.vel[idim] = ip.iniVel[idim];
    }

    //***Element search***//
    pillaz_SearchSuccessive(&ent);

    //***Initialize entity***//
    if(ent.flag==DFLAG_ENABLED && sd.enabled<ip.numMaxEnt){

      pillaz_SetElementGeometry(fp.numDim,&ent);
      pillaz_Interpolate(&ent,&flow,0.0);

      ed[sd.enabled].flag        = DFLAG_CREATED;
      ed[sd.enabled].eaddress    = ent.eaddress;
      ed[sd.enabled].node        = pillaz_FindNearestElementNode(&ent);
      ed[sd.enabled].diameter    = ent.diam;
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
}


void pillaz::pillaz_Interpolate(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  // node impact factors for interpolation
  const std::vector< int >& en = ent->edata.elmNodes;
  std::vector< double > impFac(en.size(),0.);
  pillaz_CalcNodeImpactFactors(ent,&impFac[0]);

  // interpolation of pressure and temperature
  flow->pressure = 0.;
  flow->temp     = 0.;
  for (size_t i=0; i<en.size(); ++i) {
    flow->pressure += impFac[i] * ( (1.-step)*getQuantityScalarOld(PRESSURE,   en[i]) + step*getQuantityScalar(PRESSURE,   en[i]) );
    flow->temp     += impFac[i] * ( (1.-step)*getQuantityScalarOld(TEMPERATURE,en[i]) + step*getQuantityScalar(TEMPERATURE,en[i]) );
  }

  // interpolation of velocity components, time and space derivatives
  flow->vel  .assign(fp.numDim,0.);
  flow->velDt.assign(fp.numDim,0.);
  std::vector< double >
    vnew(fp.numDim,0.),
    vold(fp.numDim,0.);
  for (size_t i=0; i<en.size(); ++i) {

    // velocity components and time derivatives
    getQuantityVector   (VELOCITY,en[i],fp.numDim,&vnew[0]);
    getQuantityVectorOld(VELOCITY,en[i],fp.numDim,&vold[0]);
    for (int d=0; d<fp.numDim; ++d) {
      flow->vel[d] += impFac[i] * ((1.-step)*vold[d] + step*vnew[d]);
      flow->velDt[d] += (fp.dtEul<1.e-20? 0. : impFac[i]/fp.dtEul * (vnew[d] - vold[d]));
    }

    // velocity space derivatives
    for (int d=0; d<fp.numDim; ++d) {
      pillaz_quantityvector_t VXDX(d==0? VELOCITY_X_D : (d==1? VELOCITY_Y_D : VELOCITY_Z_D));
      getQuantityVector   (VXDX,en[i],fp.numDim,&vnew[0]);
      getQuantityVectorOld(VXDX,en[i],fp.numDim,&vold[0]);
      for (int d2=0; d2<fp.numDim; ++d2)
        flow->velDx[d][d2] += impFac[i] * ((1.-step)*vold[d2] + step*vnew[d2]);
    }

  }

  // calculate vorticity
  pillaz_CalcVorticity(fp.numDim,flow);

  // compute relative velocity vector and its norm
  ent->normVel = 0.;
  for (int d=0; d<fp.numDim; ++d) {
    ent->relVel[d] = flow->vel[d] - ent->vel[d];
    ent->normVel += ent->relVel[d] * ent->relVel[d];
  }
  ent->normVel = sqrt(ent->normVel);

  // compute relative temperature
  ent->relTemp = flow->temp - ent->temp;
}


void pillaz::pillaz_LoadInitialDistribution(const std::string &inpString)
{
#if 0
  LOCAL_ENTITY_VARIABLES ent(fp.numDim);
  int numIni;
  char buffer[200],*search;


  //***Read data file***//
  FILE *inpFile = fopen(inpString.c_str(),"r");


  /* locate zone starting position ("^ZONE.+(dispersed)")) */
  /* set numIni (I) and fp.time (SOLUTIONTIME) */
  /* back to correct position */
  long fpos = 0L;
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


  //***Loop over entities to generate***//
  for (int ient=0, jent=0; ient<numIni && fgets(buffer,200,inpFile)!=NULL; ++ient) {


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
    ent.elm = 0;
    pillaz_SearchSuccessive(&ent);


    //***Initialize entity***//
    if (ent.flag==DFLAG_ENABLED) {
      ed[jent].flag = DFLAG_ENABLED;
      ed[jent].element = ent.elm;
      ed[jent].node = pillaz_FindNearestElementNode(&ent);
      ed[jent].diameter = ent.diam;
      ed[jent].temperature = ent.temp;
      for (int idim=0; idim<fp.numDim; ++idim){
        ed[jent].position[idim] = ent.pos[idim];
        ed[jent].velocity[idim] = ent.vel[idim];
      }
      jent++;
    }
  }


  fclose(inpFile);
#endif
}


void pillaz::pillaz_NormalizeVector(int numDim, double *a)
{
  int idim;
  double length = pillaz_CalcVectorLength(numDim,a);

  for(idim=0; idim<numDim; idim++){
    a[idim] /= length;
  }
}


double pillaz::pillaz_RandomGaussian(float m, float s)
{
  double xx1, xx2, w, yy1;
  static double yy2;
  static int use_last = 0;

  if (use_last){
    yy1 = yy2;
    use_last = 0;
  } else{
    do {
      xx1 = 2.0 * pillaz_RandomDouble() - 1.0;
      xx2 = 2.0 * pillaz_RandomDouble() - 1.0;
      w = xx1 * xx1 + xx2 * xx2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    yy1 = xx1 * w;
    yy2 = xx2 * w;
    use_last = 1;
  }

  return( m + yy1 * s );
}


void pillaz::pillaz_RandomInitialDistribution()
{
  LOCAL_ENTITY_VARIABLES ent (fp.numDim);
  LOCAL_FLOW_VARIABLES   flow(fp.numDim);


  // initializations
  ip.numIniEnt = (ip.numIniEnt>ip.numMaxEnt? ip.numMaxEnt : ip.numIniEnt);


  // calculate total volume and find maximum cell volume
  const double Vtotal = std::accumulate(volumeNod.begin(),volumeNod.end(),0.);
  double maxCellSize = 0.;
  for (size_t iz=0; iz<volumeElm.size(); ++iz)
    if (volumeElm[iz].size())
      maxCellSize = std::max(maxCellSize,*std::max_element(volumeElm[iz].begin(),volumeElm[iz].end()));


  //***Loop over entities to generate***//

  for (int ient=0; ient<ip.numIniEnt; ++ient) {

    // generate entity in a random element
    // 1. calculate total volume (Vtotal)
    // 2. pick a value within ]0.;Vtotal[ (Vpick)
    // 3. sequentially accumulate element volumes (Vacc) until overcoming Vpick
    ent.eaddress.reset();
    double
      Vpick  = pillaz_RandomDouble(0.,Vtotal),
      Vacc   = 0.;
    for (size_t iz=0; iz<fp.nInnerElements.size() && !ent.eaddress.valid(); ++iz)
      for (size_t ie=0; ie<fp.nInnerElements[iz] && !ent.eaddress.valid(); Vacc += volumeElm[iz][ie], ++ie)
        if (Vacc>Vpick)
          ent.eaddress = pillaz_elementaddress(iz,ie);
    pillaz_SetElementGeometry(fp.numDim,&ent);

    // set diameter and random position inside element
    pillaz_RandomElmPosition(ent.eaddress.izone,ent.eaddress.ielem,&ent.pos[0]);

    // initialize entity
    pillaz_Interpolate(&ent,&flow,0.0);
    ed[ient].flag        = DFLAG_CREATED;
    ed[ient].eaddress    = ent.eaddress;
    ed[ient].node        = pillaz_FindNearestElementNode(&ent);
    ed[ient].diameter    = pillaz_SetDiameter();
    ed[ient].temperature = ip.iniTempDisp;
    for (int d=0; d<fp.numDim; ++d) {
      ed[ient].position[d] = ent.pos[d];
      ed[ient].velocity[d] = flow.vel[d]+ip.iniVel[d];
    }
    ++sd.in;
  }
}


double pillaz::pillaz_RandomDouble(double min, double max)
{
  return min + (max-min) * double(rand())/double(RAND_MAX);
}


int pillaz::pillaz_RandomInteger(int min, int max)
{
  return min + (max-min)*int( double(rand()) / double(RAND_MAX) );
}


void pillaz::pillaz_ReadParameters(const XMLNode& x)
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

  ip.writeStatsFilename   = x.getAttribute< std::string >("output.statistics","pillaz.txt");
  ip.writeTecplotFilename = x.getAttribute< std::string >("output.results","pillaz.plt");
  mdd = m::Create< PILLAZ_MATERIAL_DATA >(x.getAttribute< std::string >("entities.material", "air"  ));
  mdc = m::Create< PILLAZ_MATERIAL_DATA >(x.getAttribute< std::string >("continuum.material","water"));

  ip.lagrTimeFactor = x.getAttribute< double >("lagrangeantimefactor",0.3);
  ip.errTol         = x.getAttribute< double >("errortolerance",      1.e-6);
  ip.wallElasticity = x.getAttribute< double >("wallelasticity",      1.);


  // apply corrections and set additional paramters
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


void pillaz::pillaz_SearchSuccessive(LOCAL_ENTITY_VARIABLES *ent)
{
  // successive neighbour search
  int lastNode = ent->node;
  pillaz_elementaddress lastElm;
  bool
    elmFound(false),
    leftDomain(false);
  do {

    // set geometry of current element
    pillaz_SetElementGeometry(fp.numDim,ent);

    // find the face with the minimum distance to the entity
    int eface = -1;
    double dist = 1.e99;
    pillaz_FindMinimumElementFaceDistance(fp.numDim,ent,&eface,&dist);

    // in case the minimum distance is positive, the element is found
    elmFound = (dist > 0. - ip.errTol);

    if (!elmFound) {
      // search the neighbour element in direction of the minimum face distance
      const pillaz_elementaddress neighborElm = m_elemtoelem[ent->eaddress.izone][ent->eaddress.ielem][eface];
      leftDomain = !neighborElm.valid();
      if (!leftDomain) {
        ent->eaddress = neighborElm;
        if (lastElm.valid()) {
          lastElm = neighborElm;
        }
      }
    }

  }
  while (!elmFound && !leftDomain);


  // disable entity in case of not found
  if (elmFound) {
    ent->flag = DFLAG_ENABLED;
    ent->node = pillaz_FindNearestElementNode(ent);
  }
  else {
    ent->flag     = DFLAG_LEFT;
    ent->eaddress = lastElm;
    ent->node     = lastNode;
  }
}


double pillaz::pillaz_SetDiameter()
{
  const double
    muN  = ip.iniDiam,
    sigN = ip.iniDiamStd;

  // compute diameter
  double diameter = 0.;
  do {
    if (ip.iniDiamType==0) {

      // constant diameter
      diameter = muN;

    }
    else if (ip.iniDiamType==1) {

      // gaussian normal distribution
      diameter = pillaz_RandomGaussian(muN,sigN);

    }
    else if (ip.iniDiamType==2) {

      // log-normal distribution
      const double
        sigLn = sqrt(log((sigN*sigN/(muN*muN))+1.0)),
        muLn  = log(muN)-0.5*sigLn*sigLn;
      diameter = exp(pillaz_RandomGaussian(muLn,sigLn));

    }
  }
  while (diameter<ip.errTol);
  return diameter;
}


void pillaz::pillaz_SetElementGeometry(int numDim, LOCAL_ENTITY_VARIABLES *ent)
{
  // set element number of nodes/faces and element nodes
  ent->edata.elmNodes      .resize( pillaz_getElmNNodes(getElmType(ent->eaddress.izone,ent->eaddress.ielem)) );
  ent->edata.elmNorms      .assign( pillaz_getElmNFaces(getElmType(ent->eaddress.izone,ent->eaddress.ielem)), std::vector< double >(fp.numDim,0.) );
  ent->edata.elmFaceVectors.assign( pillaz_getElmNFaces(getElmType(ent->eaddress.izone,ent->eaddress.ielem)), std::vector< double >(fp.numDim,0.) );
  getElmNodes(ent->eaddress.izone,ent->eaddress.ielem,&ent->edata.elmNodes[0]);

  // set faces and normal vectors of the element assigned to a dispersed entity
  for (size_t f=0; f<ent->edata.elmNorms.size(); ++f) {
    pillaz_getElmFaceMiddlePoint(ent->eaddress.izone,ent->eaddress.ielem,f,&(ent->edata.elmFaceVectors[f])[0]);
    for (int d=0; d<numDim; ++d)
      ent->edata.elmNorms[f][d] = m_innerelem_normals[ent->eaddress.izone][ent->eaddress.ielem][f][d];
  }
}


void pillaz::pillaz_SolveGaussSeidel(int numDim, double **mat, double *s, double *rhs)
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


void pillaz::pillaz_TerminateOnError(const std::string& errMessage)
{
  screenWarning(errMessage);
  screenWarning("fatal error!");
  throw 42;
}


void pillaz::pillaz_WallBounce(int numDim, double elasticity, LOCAL_ENTITY_VARIABLES *ent, int ibnd, int ifac)
{
  int idim;
  double unitVec[numDim];
  double wallDistancePositionVector,wallDistanceVelocityVector;


  //***Calculate unit normal of boundary segment***//
  pillaz_CalcBoundaryUnitNormal(numDim,ibnd,ifac,unitVec);


  //***Calculate normal wall distance of entity position***//
  wallDistancePositionVector = pillaz_CalcWallFaceDistance(numDim,&ent->pos[0],ibnd,ifac);


  //***Adapt velocity due to elasticity factor***//
  for(idim=0; idim<numDim; idim++){
    ent->vel[idim] *= elasticity;
  }


  //***Add position vector to velocity vector***//
  for(idim=0; idim<numDim; idim++){
    ent->vel[idim] += ent->pos[idim];
  }


  //***Calculate normal wall distance of entity velocity vector***//
  wallDistanceVelocityVector = pillaz_CalcWallFaceDistance(numDim,&ent->vel[0],ibnd,ifac);


  //***Mirror position and velocity vectors on wall segment***//
  for(idim=0; idim<numDim; idim++){
    ent->pos[idim] = ent->pos[idim]-2.0*wallDistancePositionVector*unitVec[idim];
    ent->vel[idim] = ent->vel[idim]-2.0*wallDistanceVelocityVector*unitVec[idim]-ent->pos[idim];
  }
}


void pillaz::pillaz_WriteStatsFile(const std::string &outpString, int iter, double time)
{
  std::ofstream f(outpString.c_str(),std::ios::app);
  f << '\t' << iter
    << '\t' << time
    << '\t' << sd.enabled
    << '\t' << sd.in
    << '\t' << sd.out
    << '\t' << sd.bounce
    << '\t' << sd.coll
    << '\t' << sd.lost
    << '\t' << sd.reynoldsAvg
    << '\t' << sd.nusseltAvg
    << '\t' << sd.dtLagrAvg
    << '\t' << sd.subIterAvg
    << std::endl;
}


void pillaz::pillaz_WriteTecplotFile(const std::string &outpString, int iter, double time)
{
  // write header (if no entities are active, write dummy)
  std::ofstream f(outpString.c_str(),std::ios::app);
  f << "ZONE T=\"Entities (dispersed)\" I=" << (sd.enabled? sd.enabled:1) << " AUXDATA ITER=\"" << iter << "\" SOLUTIONTIME=" << time << " DATAPACKING=POINT" << std::endl;
  if (!sd.enabled) {
    for (int idim=0; idim<fp.numDim*2+3; ++idim)
      f << " 0";
    f << std::endl;
  }

  for (int ient=0; ient<ip.numMaxEnt && sd.enabled; ++ient) {
    if (ed[ient].flag==DFLAG_ENABLED || ed[ient].flag==DFLAG_CREATED) {
      const PILLAZ_ENTITY_DATA &e = ed[ient];

      for (int d=0; d<fp.numDim; ++d)  f << ' ' << e.position[d];
      for (int d=0; d<fp.numDim; ++d)  f << ' ' << e.velocity[d];
      f << ' ' << e.diameter
        << ' ' << e.temperature
        << ' ' << 0.
        << std::endl;

    }
  }
}

