
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes the area of a triangle.
 */

double plasdriver_CalcAreaTriangle(double c[3][2])
{
  return c[0][0]*(c[1][1]-c[2][1]) + c[1][0]*(c[2][1]-c[0][1]) + c[2][0]*(c[0][1]-c[1][1]);
}

