
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity Reynolds number.
 */

double plas_CalcDispReynolds(double viscosity, double diameter, double normVel)
{
  double reynolds = diameter*normVel/viscosity;

  if(reynolds<1e-4){reynolds = 1e-4;}

  return reynolds;
}

