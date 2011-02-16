
#include "common.h"


/*
 * This file contains routines for the generation of entities,
 * be it to form an initial set of entities or to create
 * entites due to secondary phase boundary conditions.
 *
 * This function generates entites in specified spatial
 * domains of the computational space.
 */

void plas_ImposeProductionDomains(PLAS_DATA *data)
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  int ient,ipd,idim,iterate,bCtr;
  double p,s,mass=0.,p1[3],p2[3];
  double *newPos,*newDiam;
  char msg[100];

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(data->fp.numDim,&ent);
  plas_AllocateLocalFlowVar(data->fp.numDim,&flow);

  newDiam = new double[data->ip.numMaxEnt];
  newPos  = new double[data->ip.numMaxEnt*data->fp.numDim];

  //***Creation of bubbles is performed only by the master process***//

  if(plas_MpiGetRank()==0){

    bCtr = 0;

    //***Loop over production domains***//

    for(ipd=0; ipd<data->ip.numProdDom; ipd++){

      if(data->ip.prodDom[ipd]==0){continue;}

      //***Accumulate produced secondary phase mass***//

      data->rp.massResid[ipd] += data->ip.massFluxes[ipd]*data->fp.dtEul;
      if(data->rp.massResid[ipd]<=0.0){continue;}

      //***Copy production parameters to local data structure***//

      for(idim=0; idim<3; idim++){
        p1[idim] = data->ip.prodParam[ipd][idim];
        p2[idim] = data->ip.prodParam[ipd][idim+3];
      }

      //***Iteratively generate entities according to mass flux***//

      iterate = 1;
      do{

        //***Check if maximum number of entities is not exceedes***//

        if(bCtr==data->ip.numMaxEnt){break;}

        //***Set diameter***//

        newDiam[bCtr] = plas_SetDiameter(data);

        //***Calculate mass of generated entity***//

        if(data->fp.numDim==2){
          mass = data->md.rhoDisp*PI*newDiam[bCtr]*newDiam[bCtr]/4.0;
        } else if(data->fp.numDim==3){
          mass = data->md.rhoDisp*PI*newDiam[bCtr]*newDiam[bCtr]*newDiam[bCtr]/6.0;
        }

        data->rp.massResid[ipd] -= mass;
        if(data->rp.massResid[ipd]<=0.0){iterate = 0;}

        //***Compute position***//

        if(data->ip.prodDom[ipd]==1){

          //***Line production domain***//

          s = plas_RandomDouble();
          for(idim=0; idim<data->fp.numDim; idim++){
            newPos[data->fp.numDim*bCtr+idim] = p1[idim]+s*(p2[idim]-p1[idim]);
          }

        } else if(data->ip.prodDom[ipd]==2){

          //***Rectangle production domain***//

          for(idim=0; idim<data->fp.numDim; idim++){
            s = plas_RandomDouble();
            newPos[data->fp.numDim*bCtr+idim] = p1[idim]+s*(p2[idim]-p1[idim]);
          }

        } else if(data->ip.prodDom[ipd]==3){

          //***Ellipse production domain***//

          do{
            p = 0.0;
            for(idim=0; idim<data->fp.numDim; idim++){
              s = 2.0*plas_RandomDouble()-1.0;
              newPos[data->fp.numDim*bCtr+idim] = p1[idim]+s*p2[idim];
              if(p2[idim]>data->rp.errTol){
                p += (newPos[data->fp.numDim*bCtr+idim]-p1[idim])*(newPos[data->fp.numDim*bCtr+idim]-p1[idim])/(p2[idim]*p2[idim]);
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
  plas_MpiBroadcastDouble(newDiam,data->ip.numMaxEnt,0);
  plas_MpiBroadcastDouble(newPos,data->fp.numDim*data->ip.numMaxEnt,0);
  plas_MpiBroadcastDouble(data->rp.massResid,data->ip.numProdDom,0);

  //***Generate entities//

  for(ient=0; ient<bCtr; ient++){

    ent.diam = newDiam[ient];
    ent.temp = data->ip.iniTempDisp;
    for(idim=0; idim<data->fp.numDim; idim++){
      ent.pos[idim] = newPos[data->fp.numDim*ient+idim];
      ent.vel[idim] = data->ip.iniVel[idim];
    }

    //***Element search***//

    plas_SearchDomainParallel(data,&ent);

    //***Initialize entity***//

    if(ent.flag==DFLAG_ENABLED && data->sd.enabled<data->ip.numMaxEnt){

      plas_SetElementGeometry(data->fp.numDim,&ent);
      plas_Interpolate(data,&ent,&flow,0.0);

      data->ed[data->sd.enabled].flag = DFLAG_CREATED;
      data->ed[data->sd.enabled].element = ent.elm;
      data->ed[data->sd.enabled].node = plas_FindNearestElementNode(data,&ent);
      data->ed[data->sd.enabled].diameter = ent.diam;
      data->ed[data->sd.enabled].temperature = ent.temp;
      for(idim=0; idim<data->fp.numDim; idim++){
        data->ed[data->sd.enabled].position[idim] = ent.pos[idim];
        data->ed[data->sd.enabled].velocity[idim] = flow.vel[idim]+ent.vel[idim];
      }
      data->sd.enabled++;
      data->sd.in++;
    }
  }
  sprintf(msg,"imposed entities: %d\n",bCtr);
  plasinterface_screenOutput(msg);

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(data->fp.numDim,&flow);

  delete[] newDiam;
  delete[] newPos;
}

