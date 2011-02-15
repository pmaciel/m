
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity Nusselt number.
 */

double plas_CalcNusseltNumber(int evapModel, double reynolds, double spalding, double prandtl)
{
 double Nu = 2.0 + 0.6*sqrt(reynolds)*pow(prandtl,(1.0/3.0));

 if(evapModel){
   Nu *= (log(1.0+spalding)/spalding);
 }

 return Nu;
}

