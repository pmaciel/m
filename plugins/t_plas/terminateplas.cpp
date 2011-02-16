
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
    free(data->ed[ient].position);
    free(data->ed[ient].velocity);
    free(data->ed[ient].velocityOld);
  }
  free(data->ed);

  for(inod=0; inod<data->fp.numNod; inod++){
    free(data->pd[inod].volFracDx);
    free(data->pd[inod].avgVel);
    free(data->pd[inod].stdVel);
    free(data->pd[inod].dispForce);
  }
  free(data->pd);

  if(data->ip.numProdDom>0){
    for(idim=0; idim<data->ip.numProdDom; idim++){
      free(data->ip.prodParam[idim]);
    }
    free(data->ip.prodParam);
    free(data->ip.prodDom);
    free(data->ip.massFluxes);
    free(data->rp.massResid);
  }

  if(data->fp.part.numNodPairs){
    free(data->fp.part.sendRank);
    free(data->fp.part.sendNodeIdx);
    free(data->fp.part.recvRank);
    free(data->fp.part.recvNodeIdx);
  }
}

