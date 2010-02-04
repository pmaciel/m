
/* subroutines for turbulence models */

#include "common.h"

/* blending function for BSL and SST k-w models */
double F1_function(double k, double w, double nu, double y, double dkdw)
{
  const double CDkw = std::max< double >(2.*Sigw2*dkdw/w,1.e-20);
  const double a = std::min< double >(
    std::max< double >(
      sqrt(k)/(0.09*w*y),
      500.*nu/(y*y*w) ),
    4.*Sigw2*k/(CDkw*y*y) );
  return tanh(pow(a,4));
}

