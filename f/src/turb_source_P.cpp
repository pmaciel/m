
/* subroutines for turbulence models */

#include "common.h"

void turb_source_P(
  int model,
  double k,
  double turb2,
  double nu_t,
  double nu_l,
  double wd,
  double G,
  double gradkw,
  double *source_k,
  double *source_ew )
{
  double Pk;
  double fmu;
  double feps1;
  double fw;
  double Rt;
  double Ry;
  double T;
  double F1;

  Pk = nu_t*G ;

  if ( model==ITKELB ) {
    Rt = k*k/(nu_l*turb2) ;
    Ry = sqrt(k)*wd/nu_l ;

    fmu = (1.-exp(-Bmu*Ry)) ;
    fmu = fmu*fmu*(1.+Dmu/(Rt+epsilon)) ;

    feps1 = Aeps1/(fmu+epsilon) ;
    feps1 = 1. + feps1*feps1*feps1 ;
    if (iter<=turb_iterinit)
      feps1 = 1. ;

    *source_k  = Pk ;
    *source_ew = Ceps1*feps1*turb2*Pk/k ;
  }
  else if ( model==ITKENA ) {
    *source_k  = Pk ;
    *source_ew = Ceps1*turb2*Pk/k ;
  }
  else if ( model==ITKELS ) {
    T = std::max< double >(k/turb2,6.0*sqrt(nu_l/turb2)) ;

    *source_k  = Pk ;
    *source_ew = Ceps1*Pk/T ;
  }
  else if ( model==ITKWHR ) {
    *source_k  = Pk ;
    *source_ew = Cw1*turb2*Pk/k ;
  }
  else if ( model==ITKWLR ) {
    Rt = k/(nu_l*turb2) ;
    fmu = (0.15+Rt)/(6.0+Rt) ;
    fw = (0.27+Rt)/((2.7+Rt)*fmu) ;

    *source_k  = Pk ;
    *source_ew = fw*Cw1*turb2*Pk/k ;
  }
  else if ( model==ITKWBS ) {
    F1 = F1_function(k,turb2,nu_l,wd,gradkw) ;
    if (iter<=turb_iterinit)
      F1 = 1. ;
    const double Gamma = F1*Gamma1 + (1.-F1)*Gamma2 ;

    *source_k  = Pk ;
    *source_ew = (Gamma*Pk/nu_t) + (2.*Sigw2*(1.-F1)*gradkw/turb2) ;
  }
  else if ( model==ITKWSS ) {
    F1 = F1_function(k,turb2,nu_l,wd,gradkw) ;
    const double Gamma = F1*Gamma1 + (1.-F1)*Gamma2 ;

    *source_k  = Pk ;
    *source_ew = (Gamma*Pk/nu_t) + (2.*Sigw2*(1.-F1)*gradkw/turb2) ;
  }

  else if ( model==ITKWPD ) {
    Rt     = k/(nu_l*turb2) ;
    fw = 1. + 4.3*exp(-sqrt(0.6667*Rt)) ;

    *source_k  = Pk ;
    *source_ew = fw*Cw1*turb2*Pk/k + Cw*nu_t*gradkw/k ;
  }
}

