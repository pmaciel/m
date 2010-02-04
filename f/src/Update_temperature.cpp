
/* main body of implicit scheme */

#include "common.h"

void Update_temperature()
{
  /* zero jacobian, residual and vector */
  ls_aztec_scalar->reset();

  /* assemble linear system */
  IJacobian_scalar(iv_temp);
  if (dtrelax_scalar)
    for (int n=0; n<Nnode; ++n)
      ls_aztec_scalar->A(n,n) += No_dt[n]/CFL_scalar;

  /* solve linear system */
  ls_aztec_scalar->solve();

  /* calculate global residuals */
  rescalc(iv_temp,0,&(ls_aztec_scalar->m_B)[0],1);

  /* update scalar solution values */
  for (int n=0; n<Nnode; ++n)
    No_W[iv_temp][n] -= linrlx_scalar*ls_aztec_scalar->X(n);

  if (periodic) {
    for (int g=1; g<=Nbcgroup; ++g)
      if (BCgroup[g].type==IBPERE)
        for (int inb=1; inb<=BCgroup[g].nnode; ++inb)
          No_W[iv_temp][Nobg[g][inb].node] = No_W[iv_temp][Nobg[g][inb].twin];
  }
}

