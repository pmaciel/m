
#include "common.h"


/*
 * This routine terminates PLaS. It has to be called after
 * the last run of PLaS from the driving flow solver.
 */

void terminatePLaS(PLAS_DATA *data)
{
  int ient,inod,idim;

  ///***Free dynamic memory of PLaS***//

  for(ient=0; ient<data->ip.numMaxEnt; ient++){
    delete[] data->ed[ient].position;
    delete[] data->ed[ient].velocity;
    delete[] data->ed[ient].velocityOld;
  }
  delete[] data->ed;

  for(inod=0; inod<data->fp.numNod; inod++){
    delete[] data->pd[inod].volFracDx;
    delete[] data->pd[inod].avgVel;
    delete[] data->pd[inod].stdVel;
    delete[] data->pd[inod].dispForce;
  }
  delete[] data->pd;

  if(data->ip.numProdDom>0){
    for (idim=0; idim<data->ip.numProdDom; ++idim)
      delete[] data->ip.prodParam[idim];
    delete[] data->ip.prodParam;
    delete[] data->ip.prodDom;
    delete[] data->ip.massFluxes;
    delete[] data->rp.massResid;
  }

  if(data->fp.part.numNodPairs){
    delete[] data->fp.part.sendRank;
    delete[] data->fp.part.sendNodeIdx;
    delete[] data->fp.part.recvRank;
    delete[] data->fp.part.recvNodeIdx;
  }
}

