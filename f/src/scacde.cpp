
/* distribution of scalar convection/diffusion equation */

#include "common.h"

void scacde(
  int ic,
  int iv,
  local_node_struct *No_loc,
  double vol,
  int inc_min,
  int scheme,
  double diffco,
  double source,
  double coeff,
  int coeff_calc )
{
  int i,id,inc,jnc;
  double u[3], k[4], kplus[4], kminus[4], Win, sumkplus, sumkabs, recipsum, beta[4];
  double convres, diffres, diffac, peclet, res[4], q, ninj, grad[3], h, lwfac;

  for ( inc=0; inc<Nvtcell; inc++ )
    No_loc[inc].Res[iv] = 0.;

  /*** Calculate cell-average velocity ***/
  for ( id=0, q=0.; id<Ndim; id++ ) {
    u[id] = 0.;
    for ( inc=0; inc<Nvtcell; inc++ )
      u[id] += No_loc[inc].W[id+1];
    u[id] = u[id]/dNvtcell;
    q += u[id]*u[id];
  }
  q = sqrt(q);

  /*** Calculate k[i]=a.dot.n[i]/Nvtfce and convective Dt ***/
  for ( inc=0, sumkabs=0., sumkplus=0., Win=0.; inc<Nvtcell; inc++ ) {
    k[inc] = 0.;
    for ( id=0; id<Ndim; id++ )
      k[inc] += u[id]*No_loc[inc].norm[id];
    k[inc] = k[inc]/dNvtfce;

    kplus[inc] = std::max< double >(k[inc],epsilon);
    kminus[inc] = std::min< double >(k[inc],epsilon);

    sumkabs += std::abs(k[inc]);
    sumkplus += kplus[inc];

    Win += No_loc[inc].W[iv]*kminus[inc];
  }
  Win /= -sumkplus;

  /* Calculate gradient of variable */
  if (diffusion)
  for ( id=0; id<Ndim; id++ ) {
    grad[id] = 0.;
    for ( inc=0; inc<Nvtcell; inc++ )
      grad[id] += No_loc[inc].W[iv]*No_loc[inc].norm[id];
    grad[id] = grad[id]/((double)Ndim*vol);
  }

  /*** Cell length-scale = Vc/(length or area of smallest edge) ***/
  h = sqrt(No_loc[inc_min].norm2);
  if (q<1.e-10)
    h = 2.*vol/h;
  else
    h = 2.*q*vol/sumkabs;

  peclet = q*h/diffco;

//  if (peclet<2.)
//    scheme = ISSGAL;
//  else if (peclet<20.)
//    scheme = ISSLWS;

  /*** CONVECTION ***/

  for ( inc=0, convres=0.; inc<Nvtcell; inc++ )
    convres -= No_loc[inc].W[iv]*k[inc];

  if ( scheme==ISNONE || scheme==ISSFOU || scheme==ISSNSC )
    for ( inc=0; inc<Nvtcell; inc++ )
      No_loc[inc].Res[iv] += source/dNvtcell;
  else
    convres += source;

/* Compute convective nodal residual and coefficents of schemes */

  switch ( scheme ) {
    case ISNONE :
     /* No convection */

      for ( inc=0; inc<Nvtcell; inc++ )
        for ( jnc=0; jnc<Nvtcell; jnc++ )
          No_loc[inc].C[jnc] = 0.;

      break;

    case ISSFOU :                                 /* FOU finite-volume scheme */

      for ( inc=0; inc<Nvtcell; inc++ )
      for ( i=1; i<Nvtcell; i++ ) {
        jnc = (inc+i)%Nvtcell;
        No_loc[inc].Res[iv] += std::max< double >(k[inc]-k[jnc],0.)*
          (No_loc[jnc].W[iv]-No_loc[inc].W[iv])/(double)(Ndim+1);
      }

      break;

    case ISSNSC :                                 /* N-scheme */

      for ( inc=0; inc<Nvtcell; inc++ ) {
        res[inc] = kplus[inc]*(No_loc[inc].W[iv]-Win);
        No_loc[inc].Res[iv] -= res[inc];
        if ( coeff_calc ) {
          No_loc[inc].C[inc] = -kplus[inc];
          for ( i=1; i<Nvtcell; i++ ) {
            jnc = (inc+i)%Nvtcell;
            No_loc[jnc].C[inc] = -kplus[jnc]*kminus[inc]/sumkplus;
            No_loc[inc].C[jnc] = -kplus[inc]*kminus[jnc]/sumkplus;
          }
        }
      }

      break;

    case ISSGAL :                                 /* Galerkin scheme */

      for ( inc=0; inc<Nvtcell; inc++ ) {
        beta[inc] = 1./dNvtcell;
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if ( coeff_calc )
      for ( inc=0; inc<Nvtcell; inc++ ) {
        No_loc[inc].C[inc] = -beta[inc]*k[inc];
        for ( i=1; i<Nvtcell; i++ ) {
          jnc = (inc+i)%Nvtcell;
          No_loc[jnc].C[inc] = -beta[jnc]*k[inc];
          No_loc[inc].C[jnc] = -beta[inc]*k[jnc];
        }
      }

      break;

    case ISSLDA :                                 /* LDA-scheme */

      recipsum = 0.;
      if ( sumkplus>epsilon ) recipsum=1./sumkplus;

      for ( inc=0; inc<Nvtcell; inc++ ) {
        beta[inc] = kplus[inc]*recipsum;
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if ( coeff_calc )
      for ( inc=0; inc<Nvtcell; inc++ ) {
        No_loc[inc].C[inc] = -beta[inc]*k[inc];
        for ( i=1; i<Nvtcell; i++ ) {
          jnc = (inc+i)%Nvtcell;
          No_loc[jnc].C[inc] = -beta[jnc]*k[inc];
          No_loc[inc].C[jnc] = -beta[inc]*k[jnc];
        }
      }

      break;

    case ISSLWS :                                 /* Lax-Wendroff */

      /*
       * for ( inc=0; inc<Nvtcell; inc++ )
       *   sumkplus += fabs(k[inc]);
       */

      lwfac = 0.;
      if ( sumkplus>epsilon )
        lwfac = 0.5*0.5/sumkplus;
      else
        convres = 0.;

      for ( inc=0; inc<Nvtcell; inc++ ) {
        beta[inc] = 1./dNvtcell + lwfac*k[inc];
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if ( coeff_calc )
      for ( inc=0; inc<Nvtcell; inc++ ) {
        No_loc[inc].C[inc] = -beta[inc]*k[inc];
        for ( i=1; i<Nvtcell; i++ ) {
          jnc = (inc+i)%Nvtcell;
          No_loc[jnc].C[inc] = -beta[jnc]*k[inc];
          No_loc[inc].C[jnc] = -beta[inc]*k[jnc];
        }
      }

      break;

    case ISSPSI :                                 /* PSI-scheme */

      for ( inc=0, recipsum=0.; inc<Nvtcell; inc++ ) {
        res[inc] = No_loc[inc].W[iv] - Win;
        recipsum += kplus[inc]*std::min< double >(epsilon,res[inc]*convres);
      }
      recipsum = (std::abs(recipsum)>epsilon? 1/recipsum : 1.);

      for ( inc=0; inc<Nvtcell; inc++ ) {
        beta[inc]= kplus[inc]*std::min< double >(epsilon/*0.*/,res[inc]*convres)*recipsum;
        No_loc[inc].Res[iv] += beta[inc]*convres;
      }

      if ( coeff_calc )
      for ( inc=0; inc<Nvtcell; inc++ ) {
        /*
         * C[inc][inc] = -beta[inc]*k[inc];
         * for ( i=1; i<Nvtcell; i++ ) {
         *   jnc = (inc+i)%Nvtcell;
         *   C[jnc][inc] = -beta[jnc]*k[inc];
         *   C[inc][jnc] = -beta[inc]*k[jnc];
         * }
         */
        /* N-scheme Jacobian entries */
        No_loc[inc].C[inc] = -kplus[inc];
        for ( i=1; i<Nvtcell; i++ ) {
          jnc = (inc+i)%Nvtcell;
          No_loc[jnc].C[inc] = -kplus[jnc]*kminus[inc]/sumkplus;
          No_loc[inc].C[jnc] = -kplus[inc]*kminus[jnc]/sumkplus;
        }
      }

      recipsum = 0.;
      if ( sumkplus>epsilon ) recipsum=1./sumkplus;

      for ( inc=0; inc<Nvtcell; inc++ )
        beta[inc] = kplus[inc]*recipsum;

      break;

    default   : nrerror("Scalar scheme not defined in scacde.c !!!\n");
    break;
  }

  /*** DIFFUSION ***/

  /* Compute and distribute diffusive residual */
  if ( diffusion )
  for ( inc=0; inc<Nvtcell; inc++ ) {
    for ( id=0, diffres=0.; id<Ndim; id++ )
      diffres += grad[id]*No_loc[inc].norm[id];
    diffres = -diffco*diffres/dNvtfce;

    No_loc[inc].Res[iv] += diffres;
  }

  /* Add source-term and diffusion contributions to Jacobian entries */
  if ( coeff_calc ) {
    /*** SOURCE TERM COEFFICIENTS ***/
    if ( scheme==ISNONE || scheme==ISSFOU || scheme==ISSNSC )
      for ( inc=0; inc<Nvtcell; inc++ )
        for ( jnc=0; jnc<Nvtcell; jnc++ )
          No_loc[inc].C[jnc] += coeff/dNvtcell;
    else
      for ( inc=0; inc<Nvtcell; inc++ )
        for ( jnc=0; jnc<Nvtcell; jnc++ )
          No_loc[inc].C[jnc] += beta[inc]*coeff;

    /* Add diffusion terms */
    if ( diffusion ) {
      diffac = diffco/(dNvtfce*dNvtfce*vol);
      for ( inc=0; inc<Nvtcell; inc++ )
      for ( jnc=0; jnc<Nvtcell; jnc++ ) {
        for ( id=0, ninj=0.; id<Ndim; id++ )
          ninj += No_loc[inc].norm[id]*No_loc[jnc].norm[id];
        No_loc[inc].C[jnc] -= diffac*ninj;
      }
    }
  }
}
