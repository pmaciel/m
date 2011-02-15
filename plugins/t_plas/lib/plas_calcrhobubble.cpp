
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 */

double plas_CalcRhoBubble(PLAS_DATA *data, double temperature, double pressBubble)
{
  return (data->md.molarMassDisp/(Ru*temperature)*pressBubble);
}

