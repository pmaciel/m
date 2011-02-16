
#include "common.h"


/*
 * This file contains all functionality to read in data from
 * the PLaS.conf data file.
 *
 * This function reads double parameters.
 */

double plas_ReadDoubleParam(FILE *inpFile)
{
  double d;
  int ignore_i;
  char text[100],*ignore_cp;

  ignore_cp = fgets(text,100,inpFile);
  ignore_i = fscanf(inpFile,"%lf",&d);
  ignore_cp = fgets(text,100,inpFile);

  return d;
}

