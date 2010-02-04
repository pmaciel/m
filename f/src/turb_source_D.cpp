
/* subroutines for turbulence models */

#include "common.h"

void turb_source_D(
  int model,
  double k,
  double turb2,
  double nu_l,
  double wd,
  double len,
  double gradkw,
  double *source_k,
  double *source_ew,
  double *deriv_k,
  double *deriv_ew,
  double *deriv_kew,
  double *deriv_ewk )
{
  double ystar;
  double Rt;
  double Rt4;
  double feps2;
  double fk;
  double T;
  double F1;
  double Beta;

  if ( model==ITKELB ) {
    Rt = k*k/(nu_l*turb2) ;
    feps2 = 1. - exp(-Rt*Rt) ;
    if ( iter<=iter_init && !restart )
      feps2 = 1. ;

    *source_k  = -turb2 ;
    *source_ew = -Ceps2*feps2*turb2*turb2/k ;

    *deriv_k  = -2.*turb2/k ;
    *deriv_ew = -2.*Ceps2*feps2*turb2/k ;
    *deriv_kew = -1. ;

/* Use k-l model for early iterations */
    if ( iter<=iter_init && !restart ) {
      *source_k = -0.08*pow(k,1.5)/len ;
      *deriv_k  = -0.12*sqrt(k)/len ;
    }
  }
  else if ( model==ITKENA ) {
    ystar = pow(nu_l*turb2,0.25)*wd/nu_l ;
    Rt = k*k/(nu_l*turb2) ;

    feps2 = ( 1. - exp(-ystar/3.1) ) ;
    feps2 = feps2*feps2*( 1. - 0.3*exp(-Rt*Rt/42.25) ) ;
    if ( iter<=iter_init && !restart )
      feps2=1. ;

    *source_k  = -turb2 ;
    *source_ew = -Ceps2*feps2*turb2*turb2/k ;

    *deriv_k  = -2.*turb2/k ;
    *deriv_ew = -2.*Ceps2*feps2*turb2/k ;
    *deriv_kew = -1. ;
    *deriv_ewk =  Ceps2*feps2*turb2*turb2/(k*k) ;

/* Use k-l model for early iterations */
    if ( iter<=iter_init && !restart ) {
      *source_k = -0.08*pow(k,1.5)/len ;
      *deriv_k  = -0.12*sqrt(k)/len ;
    }
  }
  else if ( model==ITKELS ) {
    T = std::max< double >(k/turb2,6.0*sqrt(nu_l/turb2)) ;

    *source_k  = -turb2 ;
    *source_ew = -Ceps2*turb2/T ;

    *deriv_k  = -2.*turb2/k ;
    *deriv_ew = -Ceps2/T ;
    *deriv_kew = -1. ;

/* Use k-l model for early iterations */
    if ( iter<=iter_init && !restart ) {
      *source_k = -0.08*pow(k,1.5)/len ;
      *deriv_k  = -0.12*sqrt(k)/len ;
    }
  }
  else if ( model==ITKWHR ) {
    *source_k  = -Ck*k*turb2 ;
    *source_ew = -Cw2*turb2*turb2 ;

    *deriv_k  = -Ck*turb2 ;
    *deriv_ew = -2.*Cw2*turb2 ;
    *deriv_kew = -Ck*k ;
  }
  else if ( model==ITKWLR ) {
    Rt = k/(nu_l*turb2) ;
    Rt4 = Rt*Rt*Rt*Rt ;
    fk = (1137.8+Rt4)/(4096.0+Rt4) ;

    *source_k  = -fk*Ck*k*turb2 ;
    *source_ew = -Cw2*turb2*turb2 ;

    *deriv_k  = -fk*Ck*turb2 ;
    *deriv_ew = -2.*Cw2*turb2 ;
    *deriv_kew = -fk*Ck*k ;
  }
  else if ( model==ITKWBS ) {
    F1 = F1_function(k,turb2,nu_l,wd,gradkw) ;
    if ( iter<=iter_init && !restart )
      F1 = 1. ;
    Beta = F1*Beta1 + (1.-F1)*Beta2 ;

    *source_k  = -Betas*k*turb2 ;
    *source_ew = -Beta*turb2*turb2 ;

    *deriv_k  = -Betas*turb2 ;
    *deriv_ew = -2.*Beta*turb2 ;
    *deriv_kew = -Betas*k ;
  }
  else if ( model==ITKWSS ) {
    F1 = F1_function(k,turb2,nu_l,wd,gradkw) ;
    Beta = F1*Beta1 + (1.-F1)*Beta2 ;

    *source_k  = -Betas*k*turb2 ;
    *source_ew = -Beta*turb2*turb2 ;

    *deriv_k  = -Betas*turb2 ;
    *deriv_ew = -2.*Beta*turb2 ;
    *deriv_kew = -Betas*k ;
  }

  else if ( model==ITKWPD ) {
    Rt = k/(nu_l*turb2) ;
    fk = 1. - 0.722*exp(-1.e-4*Rt*Rt*Rt*Rt) ;

    *source_k  = -fk*Ck*k*turb2 ;
    *source_ew = -Cw2*turb2*turb2 ;

    *deriv_k  = -fk*Ck*turb2 ;
    *deriv_ew = -2.*Cw2*turb2 ;
    *deriv_kew = -fk*Ck*k ;
  }
}

