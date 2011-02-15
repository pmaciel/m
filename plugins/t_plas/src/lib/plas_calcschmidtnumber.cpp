
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity Schmidt number.
 */

double plas_CalcSchmidtNumber(PLAS_DATA *data)
{
  return ((data->md.muDisp*data->fp.kCont)/(data->md.rhoDisp*data->fp.rhoCont*data->fp.cpCont));
}

