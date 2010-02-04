
/* main body of implicit scheme */

#include "common.h"

void Update_turbulence_sequential()
{
  /* zero jacobian, residual and vector */
  ls_aztec_turb1->reset();
  ls_aztec_turb2->reset();

  /* assemble linear system (1 and 2) */
  IJacobian_turb_uncoupled();
  if (dtrelax_turb)
    for (int n=0; n<Nnode; ++n) {
      ls_aztec_turb1->A(n,n) += No_dt[n]/CFL_turb;
      ls_aztec_turb2->A(n,n) += No_dt[n]/CFL_turb;
    }

  /* solve linear system (1) */
  ls_aztec_turb1->solve();

  /* calculate global residuals (1) */
  rescalc(iv_turb1,0,&(ls_aztec_turb1->m_B)[0],1);

  /* update solution values (1) */
  int nkneg = 0;
  for (int n=0; n<Nnode; ++n) {
    const double k = No_W[iv_turb1][n] - linrlx_turb*ls_aztec_turb1->X(n);
    if (k<0.) ++nkneg; else No_W[iv_turb1][n] = k;
  }
  if (periodic) {
    for (int g=1; g<=Nbcgroup; ++g)
      if (BCgroup[g].type==IBPERE)
        for (int inb=1; inb<=BCgroup[g].nnode; ++inb)
          No_W[iv_turb1][Nobg[g][inb].node] = No_W[iv_turb1][Nobg[g][inb].twin];
  }
  std::cout << "*** negative k nodes: " << nkneg << std::endl;

  /* solve linear system (2) */
  ls_aztec_turb2->solve();

  /* calculate global residuals (2) */
  rescalc(iv_turb2,0,&(ls_aztec_turb2->m_B)[0],1);

  /* update solution values (2) */
  int neneg = 0;
  for (int n=0; n<Nnode; ++n) {
    const double e = No_W[iv_turb2][n] - linrlx_turb*ls_aztec_turb2->X(n);
    if (e<0.) ++neneg; else No_W[iv_turb2][n] = e;
  }
  if (periodic) {
    for (int g=1; g<=Nbcgroup; ++g)
      if (BCgroup[g].type==IBPERE)
        for (int inb=1; inb<=BCgroup[g].nnode; ++inb)
          No_W[iv_turb2][Nobg[g][inb].node] = No_W[iv_turb2][Nobg[g][inb].twin];
  }
  std::cout << "*** negative e nodes: " << neneg << std::endl;

}

