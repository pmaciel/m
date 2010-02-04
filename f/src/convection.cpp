
/*
 * compute and distribute cell residual, lax-wendroff system distribution scheme
 *
 * this subroutine operates on a single cell and is called within a loop over all the cells i.e. 0 to Ncell
 *
 * inputs:  No_loc    local node structure
 *          vol       cell volume
 *          inc_min   node opposite smallest face (local)
 *          celltype  marker for cell types
 *
 * outputs: Beta2     artificial compressibility factor
 *          LWfactor  Lax-Wendroff factor
 *          k[inc]    scalar influence coefficients
 *          sbuoy[iv] coefficients for buoyancy source term
 */

#include "common.h"

void convection(
  local_node_struct *No_loc,
  double vol,
  int inc_min,
  int celltype,
  double *Beta2,
  double *LWfactor,
  std::vector< double >& k,
  std::vector< double >& sbuoy,
  double ***A,
  double ***K,
  double ***B,
  std::vector< double >& Rc )
{
  int iv;
  int inc;
  int id;
  int irow;
  int jcol;
  double KiUi[10];
  double u[3];
  double lwfac;
  double betasq;
  double zeta[4];
  double sumkabs;
  double q;
  double h;
  double Tave;
  double nueff;

  /* initialize local nodal residuals to zero */
  for (inc=0; inc<Nvtcell; inc++)
    for (iv=0; iv<=Ndim; iv++)
      No_loc[inc].Res[iv] = 0.;

  /* calculate cell-average velocity */
  q = 0.;
  for (id=0; id<Ndim; id++) {
    u[id] = 0.;
    for (inc=0; inc<Nvtcell; inc++)
      u[id] += No_loc[inc].W[id+1];
    u[id] = u[id]/dNvtcell;
    q += u[id]*u[id];
  }
  q = sqrt(q);

  betasq=C2;
  *Beta2=betasq;

  /* calculate k[i]=a.dot.n[i]/Nvtfce and convective Dt */
  sumkabs = 0.;
  for (inc=0; inc<Nvtcell; inc++) {
    k[inc] = 0.;
    for (id=0; id<Ndim; id++)
      k[inc] += u[id]*No_loc[inc].norm[id];
    k[inc] = k[inc]/dNvtfce;
    No_loc[inc].Dt = std::max< double >(k[inc],0.);
    sumkabs += std::abs(k[inc]);
  }

  /* cell length-scale = Vc/(length or area of smallest edge) */
  h = sqrt(No_loc[inc_min].norm2);

  /* use Umax/2 for velocity scale */
  lwfac = 1./(Umax*h);
  /* lwfac = dt/(2*Vc) = h/(2*U*Vc) = 1/(2*U*A); */
  /* lwfac = h/(Umax*vol); */

  lwfac *= 0.5;
  *LWfactor = lwfac;

  /*
   * if (q<1.e-10)
   *   h = 2.*vol/h;
   * else
   *   h = 2.*q*vol/sumkabs;
   */

  nueff = nulam;

  /*
   * if (turmod)
   *   nueff += turb_viscosity(No_loc,turmod,celltype,vol);
   */

  /*
   * Reu = Umax*h/nueff;
   * zeta[0] = 1.;
   * zeta[0] = Reu/(1.+Reu);
   * Reu = q*h/nueff;
   * zeta[1] = Reu/(1.+Reu);
   * for (iv=2; iv<=Ndim; iv++)
   *   zeta[iv] = zeta[1];
   */

  /* use nuc=1 for p and nuc=1/2 for vel. */
  zeta[0] = 1.;
  for (iv=1; iv<=Ndim; iv++)
    zeta[iv] = 0.5;

  /* calculate jacobian matrices for cell */
  for (id=0; id<Ndim; id++) {
    for (irow=0; irow<=Ndim; irow++)
      for (jcol=0; jcol<=Ndim; jcol++)
        A[id][irow][jcol] = 0.;
    for (irow=1; irow<=Ndim; irow++)
      A[id][irow][irow] = u[id];

    /*
     * additional terms for conservative form
     *
     * for (irow=1; irow<=Ndim; irow++)
     *   A[id][irow][id+1] += u[irow-1];
     */

    A[id][id+1][0] = 1.;
    A[id][0][id+1] = betasq;
  }

  /* calculate Ki matrices for each node: Ki=ni.A/Nvtfce */
  for (inc=0; inc<Nvtcell; inc++)
    for (irow=0; irow<=Ndim; irow++)
    for (jcol=0; jcol<=Ndim; jcol++) {
      K[inc][irow][jcol] = 0.;
    for (id=0; id<Ndim; id++)
      K[inc][irow][jcol] += No_loc[inc].norm[id]*A[id][irow][jcol]/dNvtfce;
  }

  /* calculate vector residual for cell */
  for (iv=0; iv<=Ndim; iv++)
    Rc[iv] = 0.;
  for (inc=0; inc<Nvtcell; inc++) {
    for (irow=0; irow<=Ndim; irow++) {
      KiUi[irow] = 0.;
      for (jcol=0; jcol<=Ndim; jcol++)
        KiUi[irow] += K[inc][irow][jcol]*No_loc[inc].W[jcol];
      Rc[irow] -= KiUi[irow];
    }
  }

  /*
   * add -2/3*grad(k) term into momentum residuals
   *
   * if (turmod && iter>10) {
   *   for (id=0; id<Ndim; id++) {
   *     gradk[id] = 0.;
   *     for (inc=0; inc<Nvtcell; inc++)
   *       gradk[id] += No_loc[inc].W[iv_turb1]*No_loc[inc].norm[id];
   *     gradk[id] /= dNvtfce;
   *     Rc[id+1] -= 0.66667*gradk[id];
   *   }
   * }
   */

  if (periodic)
    Rc[periodic_dirn+1] -= periodic_pgrad*vol/rho;

  if (buoyancy) {
    Tave = 0.;
    for (inc=0; inc<Nvtcell; inc++)
      Tave += No_loc[inc].W[iv_temp];
    Tave = Tave/dNvtcell;
    for (iv=1; iv<=Ndim; iv++) {
      Rc[iv] -= vol*grav[iv-1]*rhofac*(Tave-To);
      sbuoy[iv]=vol*grav[iv-1]*rhofac/dNvtcell;
    }
    sbuoy[0]=sbuoy[iv_temp]=0.;
  }

  /*
   * calculate old-time nodal residuals
   * convection and pressure contribution: Rci=Bci*Kci*Rc
   */
  for (inc=0; inc<Nvtcell; inc++) {
    for (irow=0; irow<=Ndim; irow++)
      for (jcol=0; jcol<=Ndim; jcol++)
        No_loc[inc].Res[irow] += zeta[irow]*lwfac*K[inc][irow][jcol]*Rc[jcol];
    for (irow=0; irow<=Ndim; irow++)
      No_loc[inc].Res[irow] += Rc[irow]/dNvtcell;
  }

  /* fill distribution matrices (for approximate jacobians) */
  if (Jacobian==0 || Jacobian==1) {
    for (inc=0; inc<Nvtcell; inc++) {
      for (irow=0; irow<=Ndim; irow++)
        for (jcol=0; jcol<=Ndim; jcol++)
          B[inc][irow][jcol] = -zeta[irow]*lwfac*K[inc][irow][jcol];
      for (irow=0; irow<=Ndim; irow++)
        B[inc][irow][irow] += -1./dNvtcell;
    }
  }

}

