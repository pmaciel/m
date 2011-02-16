
#include <cstdlib>
#include "common.h"


/*
 * This file contains routines to generate random numbers.
 *
 * This functions generates random integers.
 */

int plas_RandomInteger(int min, int max)
{
  return min + (int)( ((double)(max - min) + 1.) * rand()/(RAND_MAX + 1.));
}
