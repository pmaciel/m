
#include "common.h"


/*
 * This file contains all functionality to read in data from
 * the PLaS.conf data file.
 *
 * This function reads integer parameters.
 */

int plas_ReadIntParam(FILE *inpFile)
{
  int i,ignore_i;
  char text[100],*ignore_cp;

  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&i);
  ignore_cp = fgets(text,100,inpFile);

  return i;
}

