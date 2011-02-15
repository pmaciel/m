
#include "common.h"


/*
 * This file includes a function to manage cellwise data of
 * the secondary phase.
 *
 * This function computes the cellwise data (e.g. the volume
 * fraction) of the secondary phase and corrects the data
 * on the boundaries of a multiprocessor job.
 */

void plas_CalcCellwiseData(PLAS_DATA *data)
{
  int irank = plas_MpiGetRank();
  int ient,inod,ielm,idim,jdim,itype,numElmNodes=0,correctNumDens[data->fp.part.numNodPairs];
  double ivol=0,xi,yi,xixi,xiyi,volFracOld[data->fp.numNod];
  double correctDiam[data->fp.part.numNodPairs],correctRespTime[data->fp.part.numNodPairs];
  double correctVel[data->fp.numDim*data->fp.part.numNodPairs];
  double correctStdDiam[data->fp.part.numNodPairs],correctStdVel[data->fp.numDim*data->fp.part.numNodPairs];
  double correctVolFracDx[data->fp.numDim*data->fp.part.numNodPairs];

  //***Initialize cellwise secondary phase data***//

  for(inod=0; inod<data->fp.numNod; inod++){
    data->pd[inod].numDens = 0;
    volFracOld[inod] = data->pd[inod].volFrac;
    data->pd[inod].volFrac = 0.0;
    data->pd[inod].volFracDt = 0.0;
    data->pd[inod].avgDiam = 0.0;
    data->pd[inod].stdDiam = 0.0;
    data->pd[inod].avgRespTime = 0.0;
    for(idim=0; idim<data->fp.numDim; idim++){
      data->pd[inod].volFracDx[idim] = 0.0;
      data->pd[inod].avgVel[idim] = 0.0;
      data->pd[inod].stdVel[idim] = 0.0;
    }
  }

  //***Assemble cellwise data by looping over all active entities***//

  for(ient=0; ient<data->ip.numMaxEnt; ient++){
    if(data->ed[ient].flag==DFLAG_ENABLED){
      inod = data->ed[ient].node;
      if(data->fp.flowSolver==FLOWSOLVER_SFELES){inod %= data->fp.numNod;}
      data->pd[inod].numDens++;
      if(data->fp.numDim==2){
        ivol = PI*data->ed[ient].diameter*data->ed[ient].diameter/4.0;
      } else if(data->fp.numDim==3){
        ivol = PI*data->ed[ient].diameter*data->ed[ient].diameter*data->ed[ient].diameter/6.0;
      }
      data->pd[inod].volFrac += ivol/plasinterface_getNodVol(inod);
      data->pd[inod].avgDiam += data->ed[ient].diameter;
      data->pd[inod].avgRespTime += (2.0*data->md.rhoDisp+data->fp.rhoCont)
                                   *data->ed[ient].diameter*data->ed[ient].diameter/(24.0*data->fp.muCont);
      for(idim=0; idim<data->fp.numDim; idim++){
        data->pd[inod].avgVel[idim] += data->ed[ient].velocity[idim];
      }
    }
  }

  //***Divide averaged quantities by number density***//

  for(inod=0; inod<data->fp.numNod; inod++){
    if(data->pd[inod].numDens>0){
      data->pd[inod].avgDiam /= data->pd[inod].numDens;
      data->pd[inod].avgRespTime /= data->pd[inod].numDens;
      for(idim=0; idim<data->fp.numDim; idim++){
        data->pd[inod].avgVel[idim] /= data->pd[inod].numDens;
      }
    }
  }

  //***Compute standard deviations of averaged quantities***//

  if(data->ip.collisionModel){

    for(ient=0; ient<data->ip.numMaxEnt; ient++){
      if(data->ed[ient].flag==DFLAG_ENABLED){
        inod = data->ed[ient].node;
        if(data->fp.flowSolver==FLOWSOLVER_SFELES){inod %= data->fp.numNod;}

        data->pd[inod].stdDiam += pow(data->ed[ient].diameter-data->pd[inod].avgDiam,2.0);
        for(idim=0; idim<data->fp.numDim; idim++){
          data->pd[inod].stdVel[idim] += pow(data->ed[ient].velocity[idim]-data->pd[inod].avgVel[idim],2.0);
        }
      }
    }

    for(inod=0; inod<data->fp.numNod; inod++){
      if(data->pd[inod].numDens>0){
        data->pd[inod].stdDiam = sqrt(data->pd[inod].stdDiam/data->pd[inod].numDens);
        for(idim=0; idim<data->fp.numDim; idim++){
          data->pd[inod].stdVel[idim] = sqrt(data->pd[inod].stdVel[idim]/data->pd[inod].numDens);
        }
      }
    }
  }

  //***Calculate volume fraction derivatrives in time and space***//

  if(data->ip.volfracCoupl){

    for(inod=0; inod<data->fp.numNod; inod++){
      data->pd[inod].volFracDt = (data->pd[inod].volFrac-volFracOld[inod])/data->fp.dtEul;
    }

    for(ielm=0; ielm<data->fp.numElm; ielm++){
      for(idim=0; idim<data->fp.numDim; idim++){

        itype = plasinterface_getElementType(ielm);
        if(itype==ELM_SIMPLEX){
          numElmNodes = data->fp.numDim+1;
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
          inod = plasinterface_getElmNode(ielm,jdim);
          xi += plasinterface_getNodCoord(inod,idim);
          yi += data->pd[inod].volFrac;
          xixi += plasinterface_getNodCoord(inod,idim)*plasinterface_getNodCoord(inod,idim);
          xiyi += plasinterface_getNodCoord(inod,idim)*data->pd[inod].volFrac;
        }
        data->pd[inod].volFracDx[idim] = (numElmNodes*xiyi-xi*yi)/(numElmNodes*xixi-xi*xi);
      }
    }
  }

  //***Gather phase data on partition boundaries***//

  for(inod=0; inod<data->fp.part.numNodPairs; inod++){

    correctNumDens[inod] = 0;
    correctDiam[inod] = 0.0;
    correctStdDiam[inod] = 0.0;
    correctRespTime[inod] = 0.0;
    for(idim=0; idim<data->fp.numDim; idim++){
      correctVel[data->fp.numDim*inod+idim] = 0.0;
      correctStdVel[data->fp.numDim*inod+idim] = 0.0;
      correctVolFracDx[data->fp.numDim*inod+idim] = 0.0;
    }

    if(data->fp.part.sendRank[inod]==irank){
      correctNumDens[inod] += data->pd[data->fp.part.sendNodeIdx[inod]].numDens;
      correctDiam[inod] += data->pd[data->fp.part.sendNodeIdx[inod]].avgDiam*correctNumDens[inod];
      correctStdDiam[inod] += data->pd[data->fp.part.sendNodeIdx[inod]].stdDiam*correctNumDens[inod];
      correctRespTime[inod] += data->pd[data->fp.part.sendNodeIdx[inod]].avgRespTime*correctNumDens[inod];
      for(idim=0; idim<data->fp.numDim; idim++){
        correctVel[data->fp.numDim*inod+idim] += data->pd[data->fp.part.sendNodeIdx[inod]].avgVel[idim]*correctNumDens[inod];
        correctStdVel[data->fp.numDim*inod+idim] += data->pd[data->fp.part.sendNodeIdx[inod]].stdVel[idim]*correctNumDens[inod];
        correctVolFracDx[data->fp.numDim*inod+idim] += data->pd[data->fp.part.sendNodeIdx[inod]].volFracDx[idim]*correctNumDens[inod];
      }
    }

    if(data->fp.part.recvRank[inod]==irank){
      correctNumDens[inod] += data->pd[data->fp.part.recvNodeIdx[inod]].numDens;
      correctDiam[inod] += data->pd[data->fp.part.recvNodeIdx[inod]].avgDiam*correctNumDens[inod];
      correctStdDiam[inod] += data->pd[data->fp.part.recvNodeIdx[inod]].stdDiam*correctNumDens[inod];
      correctRespTime[inod] += data->pd[data->fp.part.recvNodeIdx[inod]].avgRespTime*correctNumDens[inod];
      for(idim=0; idim<data->fp.numDim; idim++){
        correctVel[data->fp.numDim*inod+idim] += data->pd[data->fp.part.recvNodeIdx[inod]].avgVel[idim]*correctNumDens[inod];
        correctStdVel[data->fp.numDim*inod+idim] += data->pd[data->fp.part.recvNodeIdx[inod]].stdVel[idim]*correctNumDens[inod];
        correctVolFracDx[data->fp.numDim*inod+idim] += data->pd[data->fp.part.recvNodeIdx[inod]].volFracDx[idim]*correctNumDens[inod];
      }
    }
  }

  //***Sum phase data on partition boundaries***//

  if(data->fp.part.numNodPairs){
    plas_MpiAllSumIntArray(correctNumDens,data->fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctDiam,data->fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctStdDiam,data->fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctRespTime,data->fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctVel,data->fp.numDim*data->fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctStdVel,data->fp.numDim*data->fp.part.numNodPairs);
    plas_MpiAllSumDoubleArray(correctVolFracDx,data->fp.numDim*data->fp.part.numNodPairs);
  }

  //***Correct phase data on partition boundaries***//

  for(inod=0; inod<data->fp.part.numNodPairs; inod++){

    if(data->fp.part.sendRank[inod]==irank){
      data->pd[data->fp.part.sendNodeIdx[inod]].numDens = correctNumDens[inod];
      if(correctNumDens[inod]>0){
        data->pd[data->fp.part.sendNodeIdx[inod]].avgDiam = correctDiam[inod]/correctNumDens[inod];
        data->pd[data->fp.part.sendNodeIdx[inod]].stdDiam = correctStdDiam[inod]/correctNumDens[inod];
        data->pd[data->fp.part.sendNodeIdx[inod]].avgRespTime = correctRespTime[inod]/correctNumDens[inod];
        for(idim=0; idim<data->fp.numDim; idim++){
          data->pd[data->fp.part.sendNodeIdx[inod]].avgVel[idim] = correctVel[data->fp.numDim*inod+idim]/correctNumDens[inod];
          data->pd[data->fp.part.sendNodeIdx[inod]].stdVel[idim] = correctStdVel[data->fp.numDim*inod+idim]/correctNumDens[inod];
          data->pd[data->fp.part.sendNodeIdx[inod]].volFracDx[idim] = correctVolFracDx[data->fp.numDim*inod+idim]/correctNumDens[inod];
        }
      }
    }

    if(data->fp.part.recvRank[inod]==irank){
      data->pd[data->fp.part.recvNodeIdx[inod]].numDens = correctNumDens[inod];
      if(correctNumDens[inod]>0){
        data->pd[data->fp.part.recvNodeIdx[inod]].avgDiam = correctDiam[inod]/correctNumDens[inod];
        data->pd[data->fp.part.recvNodeIdx[inod]].stdDiam = correctStdDiam[inod]/correctNumDens[inod];
        data->pd[data->fp.part.recvNodeIdx[inod]].avgRespTime = correctRespTime[inod]/correctNumDens[inod];
        for(idim=0; idim<data->fp.numDim; idim++){
          data->pd[data->fp.part.recvNodeIdx[inod]].avgVel[idim] = correctVel[data->fp.numDim*inod+idim]/correctNumDens[inod];
          data->pd[data->fp.part.recvNodeIdx[inod]].stdVel[idim] = correctStdVel[data->fp.numDim*inod+idim]/correctNumDens[inod];
          data->pd[data->fp.part.recvNodeIdx[inod]].volFracDx[idim] = correctVolFracDx[data->fp.numDim*inod+idim]/correctNumDens[inod];
        }
      }
    }
  }
}

