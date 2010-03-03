
#include <cmath>
#include "common.h"

/* compute L1, L2 and Li error norms */
void rescalc(int iv, int iv_local, double *res, int Nsize)
{
  double resL1new = 0.,
         resL2new = 0.,
         resLinew = 0.;
  int nnint = 0;
  for (int n=0; n<Nnode; ++n) {
    if (!No_group[n]) {
      ++nnint;
      const double r = std::abs(res[n*Nsize+iv_local]);
      resL1new += r;
      resL2new += r*r;
      resLinew = std::max< double >(resLinew,r);
    }
  }
  resL1new /= (double) nnint;
  resL2new  = sqrt(resL2new/(double) nnint);

  resL1new = resL1new<epsilon? epsilon:resL1new;
  resL2new = resL2new<epsilon? epsilon:resL2new;
  resLinew = resLinew<epsilon? epsilon:resLinew;

  if (iter<=1) {
    logL1[iv] = 0.;
    logL2[iv] = 0.;
    logLi[iv] = 0.;
  }
  else {
    logL1[iv] += log10(resL1new)-log10(resL1[iv]);
    logL2[iv] += log10(resL2new)-log10(resL2[iv]);
    logLi[iv] += log10(resLinew)-log10(resLi[iv]);
  }

  resL1[iv] = resL1new;
  resL2[iv] = resL2new;
  resLi[iv] = resLinew;
}
