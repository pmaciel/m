
#include "common.h"


/*
 * This function computes a rotation matrix in 2D.
 */

void plas_CalcRotationMatrix_2D(double phi, double **m)
{
  m[0][0] = cos(phi);
  m[0][1] = -sin(phi);
  m[1][0] = sin(phi);
  m[1][1] = cos(phi);
}

