
/* subroutines for turbulence models */

#include "common.h"

void turb_init(ITid model)
{
  switch (model) {

    case ITKEHR:    // k-epsilon, high-Re with standard wall functions
    case ITKEHG:    // k-epsilon, high-Re with generalized wall functions
    case ITKE2L:    // k-epsilon, two-layer model
    case ITKELB: {  // k-epsilon, Lam-Bremhorst low-Re
      Cmu = 0.09;
      Ceps1 = 1.44;
      Ceps2 = 1.92;
      Sig1 = 1.0;
      Sig2 = 1.3;
      if (model==ITKELB) {
        Amu = 1.0;
        Bmu = 0.0165;
        Dmu = 20.5;
        Aeps1 = 0.05;
        Aeps2 = 1.0;
        Beps2 = 1.0;
      }
    } break;

    case ITKENA: {  // k-epsilon, Abe-Kondoh-Nagano low-Re
      Cmu = 0.09;
      Ceps1 = 1.5;
      Ceps2 = 1.9;
      Sig1 = 1.4;
      Sig2 = 1.4;
    } break;

    case ITKELS: {  // k-epsilon, Launder-Sharma (Durbin version)
      Cmu = 0.19;
      Ceps1 = 1.55;
      Ceps2 = 1.9;
      Sig1 = 1.0;
      Sig2 = 1.3;
      Ceta = 70.0;
      Ceta2 = Ceta*Ceta;
    } break;

    case ITKWHR:    // k-omega, Wilcox high-Re
    case ITKWLR: {  // k-omega, Wilcox low-Re
      Ck  = 0.09;
      Cw1 = 0.5555;
      Cw2 = 0.075;
      Sig1 = 2.0;
      Sig2 = 2.0;
    } break;

    case ITKWBS:    // k-omega, Menter baseline model
    case ITKWSS: {  // k-omega, Menter shear-stress transport
      Betas = 0.09;

      Beta1 = 0.075;
      Sigk1 = model==ITKWSS? 0.85 : 0.5;
      Sigw1 = 0.5;

      Beta2 = 0.0828;
      Sigk2 = 1.0;
      Sigw2 = 0.856;

      Gamma1 = Beta1/Betas-0.168*Sigw1/sqrt(Betas);
      Gamma2 = Beta2/Betas-0.168*Sigw2/sqrt(Betas);

      Cw2 = Beta1;
      A1 = model==ITKWSS? 0.31 : 0.;
    } break;

    case ITKWPD: {  // k-omega, Peng-Davidson-Holmberg
      Ck  = 0.09;
      Cw1 = 0.42;
      Cw2 = 0.075;
      Cw  = 0.75;
      Sig1 = 0.8;
      Sig2 = 1.35;
    } break;

    case ITKWWF: {  // k-omega, with wall functions
    } break;

    // dummy fallback
    case ITNONE: case ITMGKE: case ITMGKW: default: {
    } break;
  }
}

