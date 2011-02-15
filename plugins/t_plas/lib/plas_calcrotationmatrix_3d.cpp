
#include "common.h"


/*
 * This function computes a rotation matrix in 3D.
 */

void plas_CalcRotationMatrix_3D(double phi, double **m, int axis)
{
  if(axis==0){
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = cos(phi);
    m[1][2] = -sin(phi);
    m[2][0] = 0.0;
    m[2][1] = sin(phi);
    m[2][2] = cos(phi);
  }else if(axis==1){
    m[0][0] = cos(phi);
    m[0][1] = 0.0;
    m[0][2] = -sin(phi);
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[2][0] = sin(phi);
    m[2][1] = 0.0;
    m[2][2] = cos(phi);
  }else if(axis==2){
    m[0][0] = cos(phi);
    m[0][1] = -sin(phi);
    m[0][2] = 0.0;
    m[1][0] = sin(phi);
    m[1][1] = cos(phi);
    m[1][2] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
  }
}

