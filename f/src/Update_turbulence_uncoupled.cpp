
/* main body of implicit scheme */

#include "common.h"

void Update_turbulence_uncoupled()
{
  /* zero jacobian, residual and vector */
  ls_aztec_turb1->reset();
  ls_aztec_turb2->reset();

  /* assemble linear systems */
  IJacobian_turb_uncoupled();
  if (dtrelax_turb)
    for (int n=0; n<Nnode; ++n) {
      ls_aztec_turb1->A(n,n) += No_dt[n]/CFL_turb;
      ls_aztec_turb2->A(n,n) += No_dt[n]/CFL_turb;
    }

  /* solve linear systems */
  ls_aztec_turb1->solve();
  ls_aztec_turb2->solve();

  /* calculate global residuals */
  rescalc(iv_turb1,0,&(ls_aztec_turb1->m_B)[0],1);
  rescalc(iv_turb2,0,&(ls_aztec_turb2->m_B)[0],1);

  /* update scalar solution values */
  int nkneg = 0;
  int neneg = 0;
  for (int n=0; n<Nnode; ++n) {
    const double k = No_W[iv_turb1][n] - linrlx_turb*ls_aztec_turb1->X(n);
    const double e = No_W[iv_turb2][n] - linrlx_turb*ls_aztec_turb2->X(n);
    if (k<0.) ++nkneg; else No_W[iv_turb1][n] = k;
    if (e<0.) ++neneg; else No_W[iv_turb2][n] = e;
  }
  std::cout << "*** negative k/e nodes: " << nkneg << '/' << neneg << std::endl;

  if (periodic) {
    for (int g=1; g<=Nbcgroup; ++g)
      if (BCgroup[g].type==IBPERE)
        for (int inb=1; inb<=BCgroup[g].nnode; ++inb) {
          No_W[iv_turb1][Nobg[g][inb].node] = No_W[iv_turb1][Nobg[g][inb].twin];
          No_W[iv_turb2][Nobg[g][inb].node] = No_W[iv_turb2][Nobg[g][inb].twin];
        }
  }
}

