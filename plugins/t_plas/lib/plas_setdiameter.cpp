
#include "common.h"


/*
 * This file contains routines for the generation of entities,
 * be it to form an initial set of entities or to create
 * entites due to secondary phase boundary conditions.
 *
 * This function sets an entity diameter according to a
 * distribution function.
 */

double plas_SetDiameter(PLAS_DATA *data)
{
  double muN = data->ip.iniDiam;
  double sigN = data->ip.iniDiamStd;
  double diameter=0.,muLn,sigLn;

  //***Compute diameter***//

  do{

    if(data->ip.iniDiamType==0){

      //***Constant diameter***//

      diameter = muN;

    } else if(data->ip.iniDiamType==1){

      //***Gaussian normal distribution***/

      diameter = plas_RandomGaussian(muN,sigN);

    } else if(data->ip.iniDiamType==2){

      //***Log-normal distribution***//

      sigLn = sqrt(log((sigN*sigN/(muN*muN))+1.0));
      muLn = log(muN)-0.5*sigLn*sigLn;
      diameter = exp(plas_RandomGaussian(muLn,sigLn));
    }

  } while(diameter<data->rp.errTol);

  return diameter;
}

