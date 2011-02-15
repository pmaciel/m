
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity mass transfer
 * coefficient.
 */

double plas_CalcMassTransferCoeff(PLAS_DATA *data, double sherwood, double spalding)
{
  double omega = log(1.0+spalding);

  return (2.0*(data->fp.rhoCont/data->md.rhoDisp)*data->md.binaryDiffCoeff*sherwood*omega);
}

