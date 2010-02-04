
/* subroutines for turbulence models */

#include "common.h"

void turb_init(int turmod)
{
  if ( turmod==ITKEHR || turmod==ITKEHG ) {       /* high-Re k-e model */
    Cmu = 0.09;
    Ceps1 = 1.44;
    Ceps2 = 1.92;
    Sig1 = 1.0;
    Sig2 = 1.3;
  }
  else if ( turmod==ITKE2L ) {                    /* two-layer model */
    Cmu = 0.09;
    Ceps1 = 1.44;
    Ceps2 = 1.92;
    Sig1 = 1.0;
    Sig2 = 1.3;
  }
  else if ( turmod==ITKELB ) {                    /* Lam-Bremhorst model */
    Cmu = 0.09;
    Ceps1 = 1.44;
    Ceps2 = 1.92;
    Sig1 = 1.0;
    Sig2 = 1.3;
    Amu = 1.0;
    Bmu = 0.0165;
    Dmu = 20.5;
    Aeps1 = 0.05;
    Aeps2 = 1.0;
    Beps2 = 1.0;
  }
  else if ( turmod==ITKENA ) {                    /* Abe-Kondoh-Nagano model */
    Cmu = 0.09;
    Ceps1 = 1.5;
    Ceps2 = 1.9;
    Sig1 = 1.4;
    Sig2 = 1.4;
  }
  else if ( turmod==ITKELS ) {                    /* Launder-Sharma model (Durbin version) */
    Cmu = 0.19;
    Ceps1 = 1.55;
    Ceps2 = 1.9;
    Sig1 = 1.0;
    Sig2 = 1.3;
    Ceta = 70.0;
    Ceta2 = Ceta*Ceta;
  }
  else if ( turmod==ITKWHR || turmod==ITKWLR ) {  /* Wilcox k-w model */
    Ck  = 0.09;
    Cw1 = 0.5555;
    Cw2 = 0.075;
    Sig1 = 2.0;
    Sig2 = 2.0;
  }
  else if ( turmod==ITKWBS ) {                    /* k-w BSL models */
    Betas = 0.09;

    Beta1 = 0.075;
    Sigk1 = 0.5;
    Sigw1 = 0.5;

    Beta2 = 0.0828;
    Sigk2 = 1.0;
    Sigw2 = 0.856;

    Gamma1 = Beta1/Betas-0.168*Sigw1/sqrt(Betas);
    Gamma2 = Beta2/Betas-0.168*Sigw2/sqrt(Betas);

    Cw2 = Beta1;
  }
  else if ( turmod==ITKWSS ) {                    /* k-w SST models */
    Betas = 0.09;

    Beta1 = 0.075;
    Sigk1 = 0.85;
    Sigw1 = 0.5;

    Beta2 = 0.0828;
    Sigk2 = 1.0;
    Sigw2 = 0.856;

    Gamma1 = Beta1/Betas-0.168*Sigw1/sqrt(Betas);
    Gamma2 = Beta2/Betas-0.168*Sigw2/sqrt(Betas);

    Cw2 = Beta1;

    A1=0.31;
  }
  else if ( turmod==ITKWPD ) {                    /* Peng-Davidson-Holmberg k-w model */
    Ck  = 0.09;
    Cw1 = 0.42;
    Cw2 = 0.075;
    Cw  = 0.75;
    Sig1 = 0.8;
    Sig2 = 1.35;
  }
}

