
#include "common.h"


/*
 * This file contains all functionality to read in data from
 * the PLaS.conf data file.
 *
 * This function reads the PLaS.conf data file.
 */

void plas_ReadParameters(PLAS_DATA *data)
{
  int i,j,ignore_i;
  char inpText[200],fileString[200],errMessage[200],*ignore_cp;
  FILE *inpFile;

  //***Open input file***//

  inpFile = NULL;
  if (data->ip.confFilename!=NULL) {
    if (strlen(data->ip.confFilename)) {
      inpFile = fopen(data->ip.confFilename,"r");
      if (inpFile==NULL) {
        sprintf(errMessage,"Configuration file: \"%s\" was not found.\n",data->ip.confFilename);
        plasinterface_screenOutput(errMessage);
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

  data->ip.numMaxEnt = plas_ReadIntParam(inpFile);
  if(data->ip.numMaxEnt<0){
    sprintf(errMessage,"Bad value for maximum number of dispersed entities.");
    plas_TerminateOnError(errMessage);
  }

  //***Read number of initially distributed entities***//

  data->ip.numIniEnt = plas_ReadIntParam(inpFile);
  if(data->ip.numIniEnt<0){
    sprintf(errMessage,"Bad value for number of initially distributed dispersed entities.");
    plas_TerminateOnError(errMessage);
  }

  //***Read production domains and mass fluxes***//

  data->ip.numProdDom = plas_ReadIntParam(inpFile);
  if(data->ip.numProdDom<0){
    sprintf(errMessage,"Invalid number of production domains.");
    plas_TerminateOnError(errMessage);
  }

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);

  if(data->ip.numProdDom>0){

    data->ip.prodDom = (int*)calloc(data->ip.numProdDom,sizeof(int));
    if(data->ip.numProdDom>0){
      data->ip.prodParam = (double**)calloc(data->ip.numProdDom,sizeof(double*));
      data->ip.massFluxes = (double*)calloc(data->ip.numProdDom,sizeof(double));
    }

    for(i=0; i<data->ip.numProdDom; i++){
      ignore_i = fscanf(inpFile,"%d",&data->ip.prodDom[i]);
      if(data->ip.prodDom[i]<1 || data->ip.prodDom[i]>3){
        sprintf(errMessage,"Bad value for production domain type #%d.",i);
        plas_TerminateOnError(errMessage);
      }
      data->ip.prodParam[i] = (double*)calloc(6,sizeof(double));
      for(j=0; j<6; j++){
        ignore_i = fscanf(inpFile,"%lf",&data->ip.prodParam[i][j]);
      }
      ignore_i = fscanf(inpFile,"%lf",&data->ip.massFluxes[i]);
      if(data->ip.massFluxes[i]<0.0){
        sprintf(errMessage,"Bad value for mass flux #%d.",i);
        plas_TerminateOnError(errMessage);
      }
      ignore_cp = fgets(inpText,200,inpFile);
    }
  }

  //***Read diameter spectrum of dispersed entities***//

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&data->ip.iniDiamType);
  if(data->ip.iniDiamType<0 || data->ip.iniDiamType>2){
    sprintf(errMessage,"Bad value for initial diameter distribution type.");
    plas_TerminateOnError(errMessage);
  }
  ignore_i = fscanf(inpFile,"%lf",&data->ip.iniDiam);
  if(data->ip.iniDiam<1e-6){
    sprintf(errMessage,"Bad value for entity diameter.");
    plas_TerminateOnError(errMessage);
  }
  ignore_i = fscanf(inpFile,"%lf",&data->ip.iniDiamStd);
  if(data->ip.iniDiamStd<0.0){
    sprintf(errMessage,"Negative value for entity diameter standard deviation.");
    plas_TerminateOnError(errMessage);
  }
  ignore_cp = fgets(inpText,200,inpFile);

  //***Read initial velocity of dispersed entities***//

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);
  for(i=0; i<3; i++){
    ignore_i = fscanf(inpFile,"%lf",&data->ip.iniVel[i]);
  }
  ignore_cp = fgets(inpText,200,inpFile);

  //***Read temperature of dispersed entities***//

  data->ip.iniTempDisp = plas_ReadDoubleParam(inpFile);
  if(data->ip.iniTempDisp<0.0){
    sprintf(errMessage,"Bad value for dispersed phase temperature.");
    plas_TerminateOnError(errMessage);
  }

  //***Read entity material***//

  data->ip.material = plas_ReadIntParam(inpFile);

  if(data->ip.material==MAT_COPPER
     || data->ip.material==MAT_POLY){data->rp.flowType = FLOW_PARTIC;}
  else if(data->ip.material==MAT_WATER
     || data->ip.material==MAT_NHEPTANE){data->rp.flowType = FLOW_DROPLET;}
  else if(data->ip.material==MAT_HYDROGEN
     || data->ip.material==MAT_OXYGEN
     || data->ip.material==MAT_AIR){data->rp.flowType = FLOW_BUBBLY;}
  else{
    sprintf(errMessage,"No match in the entity material database.");
    plas_TerminateOnError(errMessage);
  }

  //***Read momentum back-coupling option***//

  data->ip.momentumCoupl = plas_ReadIntParam(inpFile);
  if(data->ip.momentumCoupl<0 && data->ip.momentumCoupl>2){
    sprintf(errMessage,"Bad value for momentum coupling parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read volume fraction back-coupling option***//

  data->ip.volfracCoupl = plas_ReadIntParam(inpFile);
  if(data->ip.volfracCoupl!=0 && data->ip.volfracCoupl!=1){
    sprintf(errMessage,"Bad value for mass coupling parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read energy coupling option***//

  data->ip.energyCoupl = plas_ReadIntParam(inpFile);
  if(data->ip.energyCoupl!=0 && data->ip.energyCoupl!=1){
    sprintf(errMessage,"Bad value for energy coupling parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read collision model option***//

  data->ip.collisionModel = plas_ReadIntParam(inpFile);
  if(data->ip.collisionModel<0 || data->ip.collisionModel>2){
    sprintf(errMessage,"Bad value for collision model parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read slip-shear lift force option***//

  data->ip.liftForce = plas_ReadIntParam(inpFile);
  if(data->ip.liftForce!=0 && data->ip.liftForce!=1){
    sprintf(errMessage,"Bad value for lift force parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read evaporation model option***//

  data->ip.evapModel = plas_ReadIntParam(inpFile);
  if(data->ip.evapModel!=0 && data->ip.evapModel!=1){
    sprintf(errMessage,"Bad value for evaporation model parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read saturation model option***//

  data->ip.saturModel = plas_ReadIntParam(inpFile);
  if(data->ip.saturModel!=0 && data->ip.saturModel!=1){
    sprintf(errMessage,"Bad value for saturation model parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read periodic boundaries option***//

  data->ip.perBnd = plas_ReadIntParam(inpFile);
  if(data->ip.perBnd!=0 && data->ip.perBnd!=1){
    sprintf(errMessage,"Bad value for periodic boundary parameter.");
    plas_TerminateOnError(errMessage);
  }

  //***Read gravity vector***//

  ignore_cp = fgets(inpText,200,inpFile);
  ignore_cp = fgets(inpText,200,inpFile);
  for(i=0; i<3; i++){
    ignore_i = fscanf(inpFile,"%lf",&data->ip.gravVec[i]);
  }
  ignore_cp = fgets(inpText,200,inpFile);

  //***Read restart flag***//

  data->ip.restart = plas_ReadIntParam(inpFile);
  if(data->ip.restart!=0 && data->ip.restart!=1){
    sprintf(errMessage,"Bad value for restart flag.");
    plas_TerminateOnError(errMessage);
  }

  //***Read the output statistics filename***//

  data->ip.writeStatsFilename = (char *)calloc(200,sizeof(char));
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%s",inpText);
  strncpy(data->ip.writeStatsFilename,inpText,200);
  ignore_cp = fgets(inpText,200,inpFile);
  {
    sprintf(errMessage,"Output statistics filename: \"%s\"\n",data->ip.writeStatsFilename);
    plasinterface_screenOutput(errMessage);
  }

  //***Read the output tecplot filename***//

  data->ip.writeTecplotFilename = (char *)calloc(200,sizeof(char));
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%s",inpText);
  strncpy(data->ip.writeTecplotFilename,inpText,200);
  ignore_cp = fgets(inpText,200,inpFile);
  {
    sprintf(errMessage,"Output tecplot filename: \"%s\"\n",data->ip.writeTecplotFilename);
    plasinterface_screenOutput(errMessage);
  }

  //***Read the configuration filename***//

  data->ip.confFilename = (char *)calloc(200,sizeof(char));
  ignore_cp = fgets(inpText,200,inpFile);
  ignore_i  = fscanf(inpFile,"%s",inpText);
  strncpy(data->ip.confFilename,inpText,200);
  ignore_cp = fgets(inpText,200,inpFile);
  {
    sprintf(errMessage,"Configuration filename: \"%s\"\n",data->ip.confFilename);
    plasinterface_screenOutput(errMessage);
  }

  fclose(inpFile);
}

