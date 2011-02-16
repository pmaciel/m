
#include "common.h"


/*
 * This file contains routines to generate random numbers.
 *
 * This functions generates random integers.
 */

int plas_RandomInteger(int min, int max)
{
  double r = (double)max-(double)min+1.0;

  return min + (int)(r*rand()/(RAND_MAX+1.0));
}

