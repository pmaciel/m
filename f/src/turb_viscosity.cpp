
/* subroutines for turbulence models */

#include "common.h"

double turb_viscosity(local_node_struct *No_loc, int turmod, int cell_type, double vol)
{
  double nuturb = 0.;
  double eturb  = 0.;
  double wturb  = 0.;
  double ystar;
  double dwallc = 0.;
  double len    = 0.;
  double fmu    = 1.;
  double Ry;
  double Rt;
  double T;
  double F2;
  double arg2;
  double Omega;

/* If wall distance is required calculate its cell-average value */
  if ( walldist ) {
    dwallc = 0.;
    len    = 0.;
    for (int inc=0; inc<Nvtcell; inc++ ) {
      dwallc += No_wd[No_loc[inc].node];
      len    += No_lenturb[No_loc[inc].node];
    }
    dwallc /= dNvtcell;
    len    /= dNvtcell;
  }

  double kturb = 0.;
  double turb2 = 0.;
  for (int inc=0; inc<Nvtcell; inc++ ) {
    kturb += No_loc[inc].W[iv_turb1];
    turb2 += No_loc[inc].W[iv_turb2];
  }
  kturb /= dNvtcell;
  turb2 /= dNvtcell;

  if ( (turmod/10)==ITMGKE )
    eturb=turb2;
  else if ( (turmod/10)==ITMGKW )
    wturb=turb2;

/* Calculate cell-average turbulent viscosity */
  if ( turmod==ITKEHR || turmod==ITKEHG ) {       /* high-Re k-e model */
    if ( cell_type==1 )
      nuturb = 0. ;
    else
      nuturb = Cmu*kturb*kturb/eturb ;
  }
  else if ( turmod==ITKE2L ) {                    /* two-layer model */
    Ry = sqrt(kturb)*dwallc/nulam ;
    nuturb = Cmu*kturb*kturb/eturb ;
    if ( nuturb/nulam<20.0 )
      nuturb = pow(Cmu,0.25)*sqrt(kturb)*0.41*dwallc*(1.-exp(-Ry/70.0)) ;
  }
  else if ( turmod==ITKELB ) {                    /* Lam-Bremhorst */
    Rt = kturb*kturb/(nulam*eturb) ;
    Ry = sqrt(kturb)*dwallc/nulam ;

    fmu = (1.-exp(-Bmu*Ry)) ;
    fmu = fmu*fmu*(1.+Dmu/(Rt+epsilon)) ;

    nuturb = fmu*Cmu*kturb*kturb/eturb ;
  }
  else if ( turmod==ITKENA ) {                    /* Abe-Kondoh-Nagano */
    ystar = pow(nulam*eturb,0.25)*dwallc/nulam ;
    Rt = kturb*kturb/(nulam*eturb) ;

    fmu = 1. - exp(-ystar/14.0) ;
    fmu = fmu*fmu*( 1. + 5.0*pow(Rt,-0.75)*exp(-Rt*Rt/4.0e4) ) ;

    nuturb = Cmu*fmu*kturb*kturb/eturb ;
  }
  else if ( turmod==ITKELS ) {                    /* Launder-Sharma (Durbin version) */
    kturb = std::max< double >(1.e-10,kturb) ;
    eturb = std::max< double >(1.e-10,eturb) ;
    T = std::max< double >(kturb/eturb,6.0*sqrt(nulam/eturb)) ;
    fmu = 1.-exp(-0.01*std::abs(kturb*T/nulam)) ;
    nuturb = 0.09*kturb*T*fmu ; // should be Cmu, but Cmu=0.19 for ITKELS?
  }
  else if ( turmod==ITKWHR || turmod==ITKWWF ) {  /* Standard k-w model */
    nuturb = kturb/wturb ;
  }
  else if ( turmod==ITKWBS ) {                    /* k-w BSL model */
    nuturb = kturb/wturb ;
  }
  else if ( turmod==ITKWSS ) {                    /* k-w SST model */
    double gradv[3][3];
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        gradv[i][j] = 0.;
    for (int i=0; i<Ndim; ++i) {
      for (int j=0; j<Ndim; ++j) {
        for (int inc=0 ; inc<Nvtcell ; ++inc)
          gradv[i][j] += No_loc[inc].W[i+1]*No_loc[inc].norm[j];
        gradv[i][j] = gradv[i][j]/(dNvtfce*vol);
      }
    }
    Omega = (gradv[1][0]-gradv[0][1])*(gradv[1][0]-gradv[0][1]);
    if ( Ndim==3 ) {
      Omega += (gradv[2][1]-gradv[1][2])*(gradv[2][1]-gradv[1][2]);
      Omega += (gradv[0][2]-gradv[2][0])*(gradv[0][2]-gradv[2][0]);
    }
    Omega = sqrt(Omega);

    arg2 = std::max< double >(2.0*sqrt(kturb)/(0.09*wturb*dwallc),500.0*nulam/(dwallc*dwallc*wturb)) ;
    F2 = tanh(arg2*arg2) ;
    nuturb = A1*kturb/std::max< double >(A1*wturb,Omega*F2) ;
  }
  else if ( turmod==ITKWLR ) {                    /* Low-Re k-w model */
    Rt     = kturb/(nulam*wturb) ;
    fmu = (0.15+Rt)/(6.0+Rt) ;
    nuturb = fmu*kturb/wturb ;
  }
  else if ( turmod==ITKWPD ) {                    /* Peng-Davidson-Holmberg k-w model */
    Rt     = kturb/(nulam*wturb) ;
    fmu    = pow(0.1*Rt,0.75) ;
    fmu    = 1. - exp(-fmu) ;
    fmu   *= 0.975+(0.001/Rt)*exp(-Rt*Rt/4.e4) ;
    fmu   += 0.025 ;
    nuturb = fmu*kturb/wturb ;
  }
  else
    nrerror("Turbulence model not defined !!!\n") ;

  if (turmod/10==ITMGKE && iter<=turb_iterinit && !wall_functions)
    nuturb = 1.125*fmu*sqrt(kturb)*len;

  if (turmod/10==ITMGKW && iter<=turb_iterinit)
    nuturb = 1000.0*nulam ;

  return(nuturb) ;
}

