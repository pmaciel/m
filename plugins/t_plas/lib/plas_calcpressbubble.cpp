
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 */

double plas_CalcPressBubble(PLAS_DATA *data, double diameter, double pressure)
{
  return (pressure-data->md.satPresDisp+4.0*data->md.sigDisp/diameter);
}

