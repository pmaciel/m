
/* main body of implicit scheme */

#include "common.h"

void Update_turbulence_coupled()
{
  /* zero jacobian, residual and vector */
  ls_aztec_turb->reset();

  /* assemble linear system */
  IJacobian_turb_coupled();
  if (dtrelax_turb) {
    for (int n=0; n<Nnode; ++n)
      if (!No_group[n]) {
        ls_aztec_turb->A(n,n,0,0) += No_dt[n]/CFL_turb;
        ls_aztec_turb->A(n,n,1,1) += No_dt[n]/CFL_turb;
      }
  }

  /* solve linear system */
  ls_aztec_turb->solve();

  /* calculate global residuals */
  rescalc(iv_turb1,0,&(ls_aztec_turb->m_B)[0],2);
  rescalc(iv_turb2,1,&(ls_aztec_turb->m_B)[0],2);

  /* update scalar solution values */
  int nkneg = 0;
  int neneg = 0;
  for (int n=0; n<Nnode; ++n) {
    const double k = No_W[iv_turb1][n] - linrlx_turb*ls_aztec_turb->X(n,0);
    const double e = No_W[iv_turb2][n] - linrlx_turb*ls_aztec_turb->X(n,1);
    if (k<0.) ++nkneg; else No_W[iv_turb1][n] = k;
    if (e<0.) ++neneg; else No_W[iv_turb2][n] = e;
    /*
     * No_W[iv_turb1][n] = std::max< double >(1.e-10,No_W[iv_turb1][n]);
     * No_W[iv_turb2][n] = std::max< double >(1.e-10,No_W[iv_turb2][n]);
     */
  }
  if (periodic) {
    for (int g=1; g<=Nbcgroup; ++g)
      if (BCgroup[g].type==IBPERE)
        for (int inb=1; inb<=BCgroup[g].nnode; ++inb) {
          No_W[iv_turb1][Nobg[g][inb].node] = No_W[iv_turb1][Nobg[g][inb].twin];
          No_W[iv_turb2][Nobg[g][inb].node] = No_W[iv_turb2][Nobg[g][inb].twin];
        }
  }
  std::cout << "*** negative k/e nodes: " << nkneg << '/' << neneg << std::endl;
}

