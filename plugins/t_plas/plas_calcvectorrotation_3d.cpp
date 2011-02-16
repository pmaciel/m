
#include "common.h"


/*
 * This function rotates a 3D vector.
 */

void plas_CalcVectorRotation_3D(double phi, double **m, double *a)
{
  m[0][0]=cos(phi)+(1-cos(phi))*a[0]*a[0];
  m[0][1]=(1-cos(phi))*a[0]*a[1]-sin(phi)*a[2];
  m[0][2]=(1-cos(phi))*a[0]*a[2]+sin(phi)*a[1];
  m[1][0]=(1-cos(phi))*a[0]*a[1]+sin(phi)*a[2];
  m[1][1]=cos(phi)+(1-cos(phi))*a[1]*a[1];
  m[1][2]=(1-cos(phi))*a[1]*a[2]-sin(phi)*a[0];
  m[2][0]=(1-cos(phi))*a[0]*a[2]-sin(phi)*a[1];
  m[2][1]=(1-cos(phi))*a[2]*a[1]+sin(phi)*a[0];
  m[2][2]=cos(phi)+(1-cos(phi))*a[2]*a[2];
}

