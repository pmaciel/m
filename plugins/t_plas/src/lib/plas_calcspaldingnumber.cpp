
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity Spalding number.
 */

double plas_CalcSpaldingNumber(PLAS_DATA *data, double pressure)
{
  if (data->md.satPresDisp<1.e-20 || data->md.molarMassDisp<1.e-20)
    return 0.;

  double Y_s = 1.0/(1.0+(pressure/data->md.satPresDisp-1.0)*(data->md.molarMassDispVap/data->md.molarMassDisp));
  return ((Y_s)/(1.0-Y_s));
}

