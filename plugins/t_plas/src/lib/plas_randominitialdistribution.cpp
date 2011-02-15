
#include "common.h"


/*
 * This file contains routines for the generation of entities,
 * be it to form an initial set of entities or to create
 * entites due to secondary phase boundary conditions.
 *
 * This function performs a random initial distribution of
 * dispersed entites.
 */

void plas_RandomInitialDistribution(PLAS_DATA *data)
{
  LOCAL_ENTITY_VARIABLES ent;
  LOCAL_FLOW_VARIABLES flow;
  int ient,idim,numIni;
  double p,q,rand1,rand2,rand3,rand4,setDiam,partVolume,totalVolume,elmVolume;

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(data->fp.numDim,&ent);
  plas_AllocateLocalFlowVar(data->fp.numDim,&flow);

  //***Initializations***//

  partVolume = data->fp.domainVolume;
  totalVolume = plas_MpiAllSumDouble(partVolume);
  numIni = (int)(data->ip.numIniEnt*partVolume/totalVolume);
  if(numIni>data->ip.numMaxEnt){numIni = data->ip.numMaxEnt;}

  //***Loop over entities to generate***//

  for(ient=0; ient<numIni; ient++){

    //***Generate entity in a random element***//

    do{
      ent.elm = plas_RandomInteger(0,data->fp.numElm-1);
      plas_SetElementGeometry(data->fp.numDim,&ent);
      elmVolume = plasinterface_getElmVol(ent.elm);
      p = plas_RandomDouble();
      q = elmVolume/data->fp.maxElmVolume;
    } while(p>q);

    //***Set diameter***//

    setDiam = plas_SetDiameter(data);

    //***Random position according to element type***//

    if(plasinterface_getElementType(ent.elm)==ELM_SIMPLEX){

      rand1 = plas_RandomDouble();
      rand2 = (1.0-rand1)*plas_RandomDouble();
      rand3 = (1.0-rand1-rand2)*plas_RandomDouble();
      for(idim=0; idim<data->fp.numDim; idim++){
        ent.pos[idim] = plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)
          + rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[1],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[0],idim))
          + rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[0],idim));
        if(data->fp.numDim==3){
          ent.pos[idim] +=
            rand3*(plasinterface_getNodCoord(ent.edata.elmNodes[3],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[0],idim));
        }
      }

    } else if(plasinterface_getElementType(ent.elm)==ELM_PRISM){

      rand1 = plas_RandomDouble();
      rand2 = (1.0-rand1)*plas_RandomDouble();
      rand3 = plas_RandomDouble();

      for(idim=0; idim<data->fp.numDim; idim++){
        ent.pos[idim] = rand3*plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)
          + rand3*rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[1],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[0],idim))
          + rand3*rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[0],idim))
          + (1.0-rand3)*plasinterface_getNodCoord(ent.edata.elmNodes[3],idim)
          + (1.0-rand3)*rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[4],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[3],idim))
          + (1.0-rand3)*rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[5],idim)-plasinterface_getNodCoord(ent.edata.elmNodes[3],idim));
      }

    } else if(plasinterface_getElementType(ent.elm)==ELM_QUAD){

      rand1 = plas_RandomDouble();
      rand2 = plas_RandomDouble();
      rand3 = plas_RandomDouble();

      for(idim=0; idim<data->fp.numDim; idim++){
        ent.pos[idim] = rand3*(plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)
                               + rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[1],idim)
                                        - plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)))
          + (1.0-rand3)*(plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)
                         + rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[3],idim)
                                  - plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)));
      }

    } else if(plasinterface_getElementType(ent.elm)==ELM_HEX){

      rand1 = plas_RandomDouble();
      rand2 = plas_RandomDouble();
      rand3 = plas_RandomDouble();
      rand4 = plas_RandomDouble();

      for(idim=0; idim<data->fp.numDim; idim++){
        ent.pos[idim] = rand4*(rand3*(plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)
                                      + rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[4],idim)
                                               - plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)))
                               + (1.0-rand3)*(plasinterface_getNodCoord(ent.edata.elmNodes[1],idim)
                                              + rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[5],idim)
                                                       - plasinterface_getNodCoord(ent.edata.elmNodes[1],idim))))
          +(1.0-rand4)*(rand3*(plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)
                              + rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[6],idim)
                                       - plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)))
                       + (1.0-rand3)*(plasinterface_getNodCoord(ent.edata.elmNodes[3],idim)
                         + rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[7],idim)
                                  - plasinterface_getNodCoord(ent.edata.elmNodes[3],idim))));
      }

    } else if(plasinterface_getElementType(ent.elm)==ELM_PYRAMID){

      rand1 = plas_RandomDouble();
      rand2 = plas_RandomDouble();
      rand3 = plas_RandomDouble();
      rand4 = plas_RandomDouble();

      for(idim=0; idim<data->fp.numDim; idim++){
        ent.pos[idim] = rand4*plasinterface_getNodCoord(ent.edata.elmNodes[4],idim)
          + (1.0-rand4)*(rand3*(plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)
                                + rand1*(plasinterface_getNodCoord(ent.edata.elmNodes[1],idim)
                                         - plasinterface_getNodCoord(ent.edata.elmNodes[0],idim)))
                         + (1.0-rand3)*(plasinterface_getNodCoord(ent.edata.elmNodes[2],idim)
                                        + rand2*(plasinterface_getNodCoord(ent.edata.elmNodes[3],idim)
                                                 - plasinterface_getNodCoord(ent.edata.elmNodes[2],idim))));
      }
    }

    //***Initialize entity***//

    plas_Interpolate(data,&ent,&flow,0.0);

    data->ed[ient].flag = DFLAG_CREATED;
    data->ed[ient].element = ent.elm;
    data->ed[ient].node = plas_FindNearestElementNode(data,&ent);
    data->ed[ient].diameter = setDiam;
    data->ed[ient].temperature = data->ip.iniTempDisp;
    for(idim=0; idim<data->fp.numDim; idim++){
      data->ed[ient].position[idim] = ent.pos[idim];
      data->ed[ient].velocity[idim] = flow.vel[idim]+data->ip.iniVel[idim];
    }
    data->sd.in++;
  }

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
  plas_DeallocateLocalFlowVar(data->fp.numDim,&flow);
}

