
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes the entity drag coefficient.
 */

double plas_CalcDragCoeff(int flowType, double reynolds)
{
  double dragCoeff = 0.;

  if(flowType==FLOW_BUBBLY){

    //***Empirical drag coefficient correlations for a fluid sphere***//

    if(reynolds<1.5){
      dragCoeff = 16.0/reynolds;
    } else if(reynolds>=1.5 && reynolds<80.0){
      dragCoeff = 14.9/pow(reynolds,0.78);
    } else if(reynolds>=80.0 && reynolds<1530.0){
      dragCoeff = (49.9/reynolds)*(1.0-(2.21/pow(reynolds,0.5))) + (1.17e-5)*pow(reynolds,1.675);
    } else{
      dragCoeff = 2.61;
    }

  } else if(flowType==FLOW_PARTIC || flowType==FLOW_DROPLET){

    //Empirical drag coefficient correlations for a rigid sphere***//

    if(reynolds<0.5){
      dragCoeff = 24.0/reynolds;
    } else if(reynolds>=0.5 && reynolds<1000.0){
      dragCoeff = 24.0*(1.0+0.15*pow(reynolds,0.687))/reynolds;
    } else{
      dragCoeff = 0.44;
    }
  }

  return dragCoeff;
}

