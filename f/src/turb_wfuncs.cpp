
/* subroutines for turbulence models */

#include "common.h"

void turb_wfuncs(double *res)
{
  int inwf;
  int inu;
  int id;
  int i;
  int iv;
  int nyp_near = 0;
  int nyp_far  = 0;
  int in_ypmax = -1;
  int in_ypmin = -1;
  double vdotn;
  double vt;
  double vtan[3];
  double norm[3];
  double modn;
  double utau2;
  double wdist;
  double kturb;
  double lnEystar;
  double yplus;
  double ypmax = 0.;
  double ypmin = 1.e20;

  const double cmu25 = pow(Cmu,0.25);
  const double cmu50 = sqrt(Cmu);
  const double cmu75 = pow(Cmu,0.75);

  if ( turmod==ITKEHR )
  for ( inwf=0 ; inwf<(int) WFnodes.size() ; inwf++ ) {
    inu = WFnodes[inwf].node ;

    wdist = WFnodes[inwf].dist ;
    modn  = WFnodes[inwf].modn ;
    for ( id=0 ; id<Ndim ; id++ )
      norm[id]=WFnodes[inwf].n[id]/modn ;

    for ( id=0, vdotn=0. ; id<Ndim ; id++ )
      vdotn += No_W[id+1][inu]*norm[id] ;

    for ( id=0, vt=0. ; id<Ndim ; id++ ) {
      vtan[id] = No_W[id+1][inu] - vdotn*norm[id] ;
      vt += vtan[id]*vtan[id] ;
    }
    vt = sqrt(vt) ;

    kturb = No_W[iv_turb1][inu] ;

    for ( i=1 ; i<=10 ; i++ ) {
      lnEystar = log(9.535*cmu25*wdist*sqrt(kturb)/nulam) ;
      if ( lnEystar>0. )
        utau2 = vt*0.41*cmu25*sqrt(kturb)/lnEystar ;
      else
        utau2 = vt*nulam/wdist ;
      kturb = std::max< double >(1.e-10,utau2/cmu50) ;
    }

    yplus = wdist*sqrt(utau2)/nulam ;

/* Calculate y+ statistics */

    if ( yplus>ypmax ) {
      ypmax=yplus ;
      in_ypmax=inu ;
    }

    if ( yplus<ypmin ) {
      ypmin=yplus ;
      in_ypmin=inu ;
    }

    if ( yplus>100.0 )
      nyp_far++ ;
    else if ( yplus<10.0 )
      nyp_near++ ;

/* Set values of k and epsilon */
    No_W[iv_turb1][inu]=kturb ;
    No_W[iv_turb2][inu]=cmu75*pow(kturb,1.5)/(0.41*wdist) ;

/* Add shear force to momentum equation */
    for ( iv=1 ; iv<=Ndim ; iv++ )
      res[inu*Ncoupled+iv] -= utau2*modn*(vtan[iv-1]+epsilon)/(vt+epsilon) ;
  }

  const double dNwfnode = (double) WFnodes.size();
  std::cout << "turb_wfuncs: wall-function report..." << std::endl
            << "Maximum y+=" << ypmax << " at n=" << in_ypmax << std::endl
            << "Minimum y+=" << ypmin << " at n=" << in_ypmin << std::endl
            << "N. nodes with y+>100.0: " << nyp_far  << " (" << 100.*(double)nyp_far /dNwfnode << "%)" << std::endl
            << "N. nodes with y+<10.0:  " << nyp_near << " (" << 100.*(double)nyp_near/dNwfnode << "%)" << std::endl
            << "turb_wfuncs: wall-function report." << std::endl;
}

