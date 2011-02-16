
#include "common.h"


/*
 * This file contains routines to generate random numbers.
 *
 * This functions generates random doubles between 0 and 1.
 */

double plas_RandomDouble()
{
  int r = plas_RandomInteger(0,(int)1e5);

  return (double)(r/1e5);
}

