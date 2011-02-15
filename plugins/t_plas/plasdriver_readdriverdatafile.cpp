

#include "plas_driver.h"


/*
 * This file contains read and write routines of the PLaS
 * driver program.
 *
 * This function reads the driver data from an input file.
 */

void plasdriver_ReadDriverDataFile(DRIVER_PARAMETERS *dparam)
{
  int i,ignore_i;
  char fileString[100],errMessage[100],text[100],*ignore_cp;
  FILE *inpFile;

  //***Check if file exists***//

  sprintf(fileString,"./plas.driver");
  if(fopen(fileString,"r")==0){
    sprintf(errMessage,"File \"%s\" was not found.",fileString);
    plas_TerminateOnError(errMessage);
  }

  //***Openm file***//

  inpFile = fopen(fileString,"r");

  //***Read case name***//

  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%s",dparam->gridString);
  printf("   Case name is \"%s\"\n",dparam->gridString);
  ignore_cp = fgets(text,100,inpFile);

  //***Read boundaries***//

  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&dparam->numBnd);
  printf("   Number of boundaries is %d\n",dparam->numBnd);
  dparam->bnd = (int*)calloc(dparam->numBnd,sizeof(int));
  ignore_cp  = fgets(text,100,inpFile);
  for(i=0; i<dparam->numBnd; i++){
    ignore_i = fscanf(inpFile,"%d",&dparam->bnd[i]);
  }
  ignore_cp  = fgets(text,100,inpFile);

  //***Read number of iterations***//

  dparam->numIter = plas_ReadIntParam(inpFile);
  if(dparam->numIter<0){
    sprintf(errMessage,"Bad value for numer of iterations.");
    plas_TerminateOnError(errMessage);
  }
  printf("   Number of iterations is %d\n",dparam->numIter);

  //***Read Eulerian time step size***//

  dparam->dtEul = plas_ReadDoubleParam(inpFile);
  if(dparam->dtEul<=0.0){
    sprintf(errMessage,"Bad value for Eulerian time step size.");
    plas_TerminateOnError(errMessage);
  }
  printf("   Eulerian time step is %f\n",dparam->dtEul);

  //***Read flow medium***//

  dparam->material = plas_ReadDoubleParam(inpFile);
  if(dparam->material!=AIR && dparam->material!=WATER && dparam->material!=NITROGEN){
    sprintf(errMessage,"Bad value for material identifier.");
    plas_TerminateOnError(errMessage);
  }
  if(dparam->material==AIR){
    printf("   Flow medium is AIR\n");
  } else if(dparam->material==WATER){
    printf("   Flow medium is WATER\n");
  } else if(dparam->material==NITROGEN){
    printf("   Flow medium is NITROGEN\n");
  }

  fclose(inpFile);
}

