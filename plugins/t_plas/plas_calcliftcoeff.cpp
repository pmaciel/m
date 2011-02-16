
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity lift coefficient.
 */

double plas_CalcLiftCoeff(int flowType)
{
  if (flowType==FLOW_PARTIC || flowType==FLOW_DROPLET)
    return 0.;
  else if (flowType==FLOW_BUBBLY)
    return 0.53;
  return 0.;
}

