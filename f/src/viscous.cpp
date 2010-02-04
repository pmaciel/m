
/* galerkin treatment of viscous terms */

#include "common.h"

void viscous(local_node_struct *No_loc, double vol, int celltype)
{
  double gradv[4][3];
  double tau[3][3];

  /* gradients of velocity components */
  for (int iv=1; iv<=Ndim; ++iv)
    for (int id=0; id<Ndim; ++id) {
    gradv[iv][id] = 0.;
    for (int inc=0; inc<Nvtcell; ++inc)
      gradv[iv][id] += No_loc[inc].W[iv]*No_loc[inc].norm[id];
    gradv[iv][id] = gradv[iv][id]/(dNvtfce*vol);
  }

  /* effective viscosity */
  double nueff = nulam;
  if (turmod) {
    nueff += turb_viscosity(No_loc,turmod,celltype,vol);
    if (wall_functions && celltype==1)
      nueff = 0.;
  }

  for (int id=0; id<Ndim; ++id)
    for (int jd=0; jd<Ndim; ++jd)
      tau[id][jd] = nueff*(gradv[id+1][jd]+gradv[jd+1][id]);

  double viscfac = nueff/dNvtfce;

  /* diffusive Dt */
  for (int inc=0; inc<Nvtcell; ++inc)
    No_loc[inc].Dt += viscfac*No_loc[inc].norm2/(dNvtfce*vol);

  /* assemble and distribute viscous residual to nodes */
  for (int iv=1; iv<=Ndim; ++iv) {
    for (int inc=0; inc<Nvtcell; ++inc) {
      double viscres = 0.;
      for (int id=0; id<Ndim; ++id)
        viscres += tau[iv-1][id]*No_loc[inc].norm[id];
      No_loc[inc].Res[iv] -= viscres/dNvtfce;
    }
  }

  /* add jacobian entries */
  if (Jacobian<2) {

    viscfac /= dNvtfce*vol;
    for (int inc=0; inc<Nvtcell; ++inc) {
      for (int jnc=0; jnc<Nvtcell; ++jnc) {

        double ninj = 0.;
        for (int id=0; id<Ndim; ++id)
          ninj += No_loc[inc].norm[id]*No_loc[jnc].norm[id];
        for (int iv=1; iv<=Ndim; ++iv)
          ls_aztec_coupled->A(No_loc[inc].node,No_loc[jnc].node,iv,iv) += -viscfac*ninj;

        /* instead of the previous one
         * for (iv=1; iv<=Ndim; ++iv)
         *   for (jv=1; jv<=Ndim; ++jv)
         *     ls_aztec_coupled->A(No_loc[inc].node,No_loc[jnc].node,iv,jv) += -viscfac*No_loc[inc].norm[jv-1]*No_loc[jnc].norm[iv-1];
         */

      }
    }
  }
}

