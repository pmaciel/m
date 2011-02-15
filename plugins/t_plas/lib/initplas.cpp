
#include "common.h"


/*
 * This function initializes the PLaS solver. It has to be
 * called before doing any run of PLaS from the driving flow
 * solver.
 */

void initPLaS(PLAS_DATA *data)
{
  int ient,inod;
  FILE *inpFile;

  //***Flow solver parameters that have to be set only once***//

  plasinterface_setFlowSolverParamOnInit(&(data->fp));

  //***Set partitioning data structure in case of parallel computation***//

  plasinterface_setPartitioningData(&(data->fp.part));

  //***Read parameters from input file and broadcast them***//

  if (!plas_MpiGetRank()) {
    plas_ReadParameters(data);
  }
  plas_MpiBarrier();
  plas_BroadcastParameters(data);
  plas_MpiBarrier();

  //***Allocate memory for dispersed phase data***//

  data->ed = (PLAS_ENTITY_DATA*)calloc(data->ip.numMaxEnt,sizeof(PLAS_ENTITY_DATA));
  for(ient=0; ient<data->ip.numMaxEnt; ient++){
    data->ed[ient].flag = DFLAG_DISABLED;
    data->ed[ient].position = (double*)calloc(data->fp.numDim,sizeof(double));
    data->ed[ient].velocity = (double*)calloc(data->fp.numDim,sizeof(double));
    data->ed[ient].velocityOld = (double*)calloc(data->fp.numDim,sizeof(double));
    data->ed[ient].node = -1;
    data->ed[ient].element = -1;
  }

  data->pd = (PLAS_PHASE_DATA*)calloc(data->fp.numNod,sizeof(PLAS_PHASE_DATA));
  for(inod=0; inod<data->fp.numNod; inod++){
    data->pd[inod].volFracDx = (double*)calloc(data->fp.numDim,sizeof(double));
    data->pd[inod].avgVel = (double*)calloc(data->fp.numDim,sizeof(double));
    data->pd[inod].stdVel = (double*)calloc(data->fp.numDim,sizeof(double));
    data->pd[inod].dispForce = (double*)calloc(data->fp.numUnk,sizeof(double));
  }

  //***Initialize fixed runtime parameters***//

  data->rp.lagrTimeFactor = 0.3;
  data->rp.errTol = 1e-6;
  data->rp.wallElasticity = 1.0;

  //***Initialize variable runtime parameters***//

  if(data->ip.numProdDom>0){
    data->rp.massResid = (double*)calloc(data->ip.numProdDom,sizeof(double));
  }

  //***Initialize material data from database***//

  plas_CalcMaterialData(data,data->ip.iniTempDisp,plasinterface_getPressure(0));

  //***Initial entity distribution***//

  /* set output files */
  inpFile = fopen(data->ip.writeTecplotFilename,"r");
  if (inpFile==NULL || !data->ip.restart) {
    data->ip.restart = 0;
    plas_CreateTecplotFile(data,data->ip.writeTecplotFilename);
  }
  else
    fclose(inpFile);
  inpFile = fopen(data->ip.writeStatsFilename,"r");
  if (inpFile==NULL) {
    plas_CreateStatsFile(data,data->ip.writeStatsFilename);
  }
  else
    fclose(inpFile);

  /* loaded or random distribution */
  if (data->ip.restart)
    plas_LoadInitialDistribution(data,data->ip.writeTecplotFilename);
  else if (data->ip.numIniEnt>0)
    plas_RandomInitialDistribution(data);

  plas_MpiBarrier();
}

