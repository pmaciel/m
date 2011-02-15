
#include <cmath>

#include "plas_driver.h"


/*
 * This file contains functionality to define or read a
 * steady-state primary flow field.
 *
 * This function initializes a steady-state flow field.
 */

void plasdriver_InitFlowField(DRIVER_GAMBIT_MESH *dmesh, int material, DRIVER_FLOW_FIELD *dflow)
{
  // Define pressure, velocity and temperature
  double p = 101325.0;
  double u = 1.0;
  double v = 0.0;
  double w = 0.0;
  double T = 373.124;

  if (material==AIR) {

    // Material properties of air
    dflow->rho = 1.225;
    dflow->mu  = 1.7894e-5;
    dflow->cp  = 1.006;
    dflow->k   = 0.0242;

  }
  else if (material==WATER) {

    // Material properties of water
    dflow->rho = 998.2;
    dflow->mu  = 1.003e-3;
    dflow->cp  = 4.182;
    dflow->k   = 0.6;

  }
  else if (material==NITROGEN) {

    // Material properties of nitrogen
    dflow->rho = 1.25;
    dflow->mu  = (6.5592e-7*pow(T,0.6081))/(1.0+54.715/T);
    dflow->cp  = (6.50+0.001*T)*4.184/(28.01e-3);
    dflow->k   = 2.5*(6.5592e-7*pow(T,0.6081))/(1.0+54.715/T)*((6.50+0.001*T)*4.184/(28.01e-3)-8.314)/28.01;
  }

  // Impose flow variables
  for (int n=0; n<dmesh->numNod; ++n) {
    dflow->p[n]    = p;
    dflow->u[n][0] = u;
    dflow->u[n][1] = v;
    dflow->u[n][2] = w;
    dflow->T[n]    = T;
  }
}

