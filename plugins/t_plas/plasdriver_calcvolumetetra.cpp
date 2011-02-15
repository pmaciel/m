
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes the volume of a tetrahedron.
 */

double plasdriver_CalcVolumeTetra(double c[4][3])
{
  double v1[3],v2[3],v3[3],v4[3];

  v1[0] = c[1][0]-c[0][0];
  v1[1] = c[1][1]-c[0][1];
  v1[2] = c[1][2]-c[0][2];
  v2[0] = c[2][0]-c[0][0];
  v2[1] = c[2][1]-c[0][1];
  v2[2] = c[2][2]-c[0][2];

  plas_CalcCrossProduct_3D(v3,v1,v2);

  v4[0] = c[3][0]-c[0][0];
  v4[1] = c[3][1]-c[0][1];
  v4[2] = c[3][2]-c[0][2];

  return (plas_CalcVectScalarProduct(3,v3,v4)/6.0);
}

