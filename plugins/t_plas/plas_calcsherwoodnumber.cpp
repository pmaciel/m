
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity Sherwood number.
 */

double plas_CalcSherwoodNumber(int evapModel, double reynolds, double schmidt, double spalding)
{
  double Sh = 2.0 + 0.6*pow(reynolds,(1.0/2.0))*pow(schmidt,(1.0/3.0));

  if(evapModel){
    double F_M = pow(1.0+spalding,0.7)*log(1.0+spalding)/spalding;
    Sh = 2.0+(Sh-2.0)/F_M;
  }

  return Sh;
}

