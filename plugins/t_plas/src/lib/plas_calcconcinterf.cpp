
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the concentation in bubble's surface.
 */

double plas_CalcConcInterf(PLAS_DATA *data, double pressBubble)
{
  return (pressBubble*data->md.HeDisp);
}

