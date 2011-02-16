
#include "common.h"


/*
 * This file contains all functionality to read in data from
 * the PLaS.conf data file.
 *
 * This function broadcasts input data for a multi=processor
 * environment by the MPI Broadcast command.
 */

void plas_BroadcastParameters(PLAS_DATA *data)
{
  int irank = plas_MpiGetRank();
  int i;

  //***Broadcast static data***//

  plas_MpiBroadcastInt(&data->ip.numMaxEnt,1,0);
  plas_MpiBroadcastInt(&data->ip.numIniEnt,1,0);
  plas_MpiBroadcastInt(&data->ip.numProdDom,1,0);
  plas_MpiBroadcastInt(&data->ip.iniDiamType,1,0);
  plas_MpiBroadcastDouble(&data->ip.iniDiam,1,0);
  plas_MpiBroadcastDouble(&data->ip.iniDiamStd,1,0);
  plas_MpiBroadcastDouble(data->ip.iniVel,3,0);
  plas_MpiBroadcastDouble(&data->ip.iniTempDisp,1,0);
  plas_MpiBroadcastInt(&data->ip.material,1,0);
  plas_MpiBroadcastInt(&data->rp.flowType,1,0);
  plas_MpiBroadcastInt(&data->ip.momentumCoupl,1,0);
  plas_MpiBroadcastInt(&data->ip.volfracCoupl,1,0);
  plas_MpiBroadcastInt(&data->ip.energyCoupl,1,0);
  plas_MpiBroadcastInt(&data->ip.collisionModel,1,0);
  plas_MpiBroadcastInt(&data->ip.liftForce,1,0);
  plas_MpiBroadcastInt(&data->ip.evapModel,1,0);
  plas_MpiBroadcastInt(&data->ip.perBnd,1,0);
  plas_MpiBroadcastDouble(data->ip.gravVec,3,0);
  plas_MpiBroadcastInt(&data->ip.restart,1,0);

  /*
   * ip.writeTecplotFilename and ip.writeStatsFilename don't need
   * broadcasting because only rank 0 writes them, and ip.confFilename
   * doesn't need reading (again)
   */

  //***Broadcast allocatable data***//

  if (data->ip.numProdDom>0) {
    if (irank!=0)
      data->ip.prodDom = new int[data->ip.numProdDom];
    plas_MpiBroadcastInt(data->ip.prodDom,data->ip.numProdDom,0);
    if (irank!=0) {
      data->ip.prodParam = new double*[data->ip.numProdDom];
      for (i=0; i<data->ip.numProdDom; ++i)
        data->ip.prodParam[i] = new double[6];
    }
    for (i=0; i<data->ip.numProdDom; ++i)
      plas_MpiBroadcastDouble(data->ip.prodParam[i],6,0);
    if (irank!=0)
      data->ip.massFluxes = new double[data->ip.numProdDom];
    plas_MpiBroadcastDouble(data->ip.massFluxes,data->ip.numProdDom,0);
  }
}

