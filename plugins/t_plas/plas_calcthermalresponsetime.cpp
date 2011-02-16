
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the thermal response time of the
 * dispersed entity.
 */

double plas_CalcThermalResponseTime(PLAS_DATA *data, double diameter)
{
  return (1.0/(12.0*data->fp.kCont))*(data->md.rhoDisp*diameter*diameter*data->md.cpDisp);
}

