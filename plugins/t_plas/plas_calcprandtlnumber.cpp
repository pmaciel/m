
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity Prandtl number.
 */

double plas_CalcPrandtlNumber(PLAS_DATA *data)
{
  return (data->fp.muCont*data->fp.cpCont/data->fp.kCont);
}

