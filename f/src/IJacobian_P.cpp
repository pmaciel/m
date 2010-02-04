
/* picard jacobian calculation */

#include <fstream>
#include "common.h"

void IJacobian_P()
{
  int inc_min;
  double vol;
  double betasq;
  double lwfac;

  double Cscalar[4][4];
  double Cturb1[4][4];
  double Cturb2[4][4];
  local_node_struct No_local[4];


  std::cout << "IJacobian_P..." << std::endl;


  /* allocate memory for temporary cell arrays */
  std::vector< double > k    (Nvtcell,    0.);
  std::vector< double > sbuoy(Ncoupled+1, 0.);
  std::vector< double > Rc   (1+Ndim,     0.);
  double ***A = d3tensor(Ndim,   Ndim+1,  Ndim+1);
  double ***K = d3tensor(Nvtcell,Ndim+1,  Ndim+1);
  double ***B = d3tensor(Nvtcell,Ndim+1,  Ndim+1);
  double ***J = d3tensor(Nvtcell,Ncoupled,Ncoupled);


  /* loop over cells */
  for (int ic=0; ic<Ncell; ++ic) {

    /* compute cell normals and volume */
    cellgeom(ic,No_local,&vol,&inc_min);

    /* switch to modified connectivity if periodic */
    if (periodic) {
      e2n.swap(e2n_periodic);
      for (int inc=0; inc<Nvtcell; ++inc)
        No_local[inc].node = e2n[ic].n[inc];
    }

    /* copy nodal values from global to local structure */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv]=No_W[iv][e2n[ic].n[inc]];

    /* calculate convective contributions */
    convection(No_local,vol,inc_min,Ce_type[ic],&betasq,&lwfac,k,sbuoy,A,K,B,Rc);

    /* calculate viscous contributions */
    viscous(No_local,vol,Ce_type[ic]);

    /* calculate scalar equation contributions */
    if (temperature && scalar_coupling) {
      scacde(ic,iv_temp,No_local,vol,inc_min,scaconv,nulam/prandtl,0.,0.,1);

      for (int inc=0; inc<Nvtcell; ++inc)
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          Cscalar[inc][jnc]=No_local[inc].C[jnc];
    }

    /*
     * if (turmod && turbulence_coupling==2) {
     *   nuturb = turb_viscosity(No_local,turmod,Ce_type[ic],vol);
     *   turb_source_cell(No_local,vol,turmod,&tsrce1,&tsrce2,&tcoeff1,&tcoeff2);
     *
     *   scacde(ic,iv_turb1,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig1,tsrce1,tcoeff1,1);
     *   for (int inc=0; inc<Nvtcell; ++inc)
     *     for (int jnc=0; jnc<Nvtcell; ++jnc)
     *       Cturb1[inc][jnc]=No_local[inc].C[jnc];
     *
     *   scacde(ic,iv_turb2,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig2,tsrce2,tcoeff2,1);
     *   for (int inc=0; inc<Nvtcell; ++inc)
     *     for (int jnc=0; jnc<Nvtcell; ++jnc)
     *       Cturb2[inc][jnc]=No_local[inc].C[jnc];
     * }
     */

    /* add residuals and dt's into global No structure */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int inu = e2n[ic].n[inc];
      for (int iv=0; iv<Ncoupled; ++iv)
        ls_aztec_coupled->B(inu,iv) += No_local[inc].Res[iv];
      No_dt[inu] += No_local[inc].Dt;
    }


    /* fill jacobian matrix */
    for (int inc=0; inc<Nvtcell; ++inc) {

      /* zero entries */
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        for (int irow=0; irow<Ncoupled; ++irow)
          for (int jcol=0; jcol<Ncoupled; ++jcol)
            J[jnc][irow][jcol] = 0.;

      /* calculate Bi*Kj entries for node (row) i */
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        for (int irow=0; irow<=Ndim; ++irow)
          for (int jcol=0; jcol<=Ndim; ++jcol)
            for (int ij=0; ij<=Ndim; ++ij)
              J[jnc][irow][jcol] += B[inc][irow][ij]*K[jnc][ij][jcol];

      /* add jacobian entries for dependence of T on T (from scacde) */
      if (temperature && scalar_coupling)
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          J[jnc][iv_temp][iv_temp] = Cscalar[inc][jnc];

      /* add jacobian entries for self-dependence of turbulence quantities */
      if (turmod && turbulence_coupling==2)
      for (int jnc=0; jnc<Nvtcell; ++jnc) {
        J[jnc][iv_turb1][iv_turb1] = Cturb1[inc][jnc];
        J[jnc][iv_turb2][iv_turb2] = Cturb2[inc][jnc];
      }

      /* add jacobian entries for dependence of turbulence on V (uses N-scheme coeffs) */
      if (turmod && turbulence_coupling==2) {
        double turb1in   = 0.;
        double turb2in   = 0.;
        double sumkminus = 0.;
        for (int jnc=0; jnc<Nvtcell; ++jnc) {
          turb1in += std::min< double >(0.,k[jnc])*No_local[jnc].W[iv_turb1];
          turb2in += std::min< double >(0.,k[jnc])*No_local[jnc].W[iv_turb2];
          sumkminus += std::min< double >(0.,k[jnc]);
        }
        turb1in /= sumkminus;
        turb2in /= sumkminus;
        if (k[inc]>0.) {
          for (int jnc=0; jnc<Nvtcell; ++jnc)
            for (int iv=1; iv<=Ndim; ++iv) {
              J[jnc][iv_turb1][iv] -= No_local[inc].norm[iv-1]/(dNvtfce*dNvtcell)*(No_local[inc].W[iv_turb1]-turb1in);
              J[jnc][iv_turb2][iv] -= No_local[inc].norm[iv-1]/(dNvtfce*dNvtcell)*(No_local[inc].W[iv_turb2]-turb2in);
            }
        }
      }

      /* add jacobian entries for dependence of T on V (for N-scheme or LWS) */
      if (temperature && scalar_coupling) {
        if (scaconv==ISSNSC || scaconv==ISSPSI) {
          double Tin       = 0.;
          double sumkminus = 0.;
          for (int jnc=0; jnc<Nvtcell; ++jnc) {
            Tin += std::min< double >(0.,k[jnc])*No_local[jnc].W[iv_temp];
            sumkminus += std::min< double >(0.,k[jnc]);
          }
          Tin /= sumkminus;
          if (k[inc]>0.) {
            for (int jnc=0; jnc<Nvtcell; ++jnc)
              for (int iv=1; iv<=Ndim; ++iv)
                J[jnc][iv_temp][iv] -= No_local[inc].norm[iv-1]/(dNvtfce*dNvtcell)*(No_local[inc].W[iv_temp]-Tin);
          }
        }
        else if (scaconv==ISSLWS) {
          for (int jnc=0; jnc<Nvtcell; ++jnc)
            for (int iv=1; iv<=Ndim; ++iv)
              for (int knc=0; knc<Nvtcell; ++knc) {
                J[jnc][iv_temp][iv]-=(1./dNvtcell+lwfac*k[inc])*No_local[knc].norm[iv-1]*No_local[knc].W[iv_temp]/(dNvtcell*dNvtfce);
                J[jnc][iv_temp][iv]-=lwfac*No_local[inc].norm[iv-1]*No_local[knc].W[iv_temp]*k[knc]/(dNvtcell*dNvtfce);
              }
        }
      }

      /* add jacobian entries for dependence of p and V on T (buoyancy) */
      if (buoyancy && scalar_coupling) {
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          for (int iv=1; iv<=Ndim; ++iv) {
            J[jnc][0][iv_temp]  -= (lwfac*betasq*No_local[jnc].norm[iv-1]/dNvtfce)*sbuoy[iv];
            J[jnc][iv][iv_temp] -= (1./dNvtcell+lwfac*k[jnc])*sbuoy[iv];
          }
      }

      /* add entries in row inc to Jacobian */
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        for (int i=0; i<Ncoupled; ++i)
          for (int j=0; j<Ncoupled; ++j)
            ls_aztec_coupled->A(No_local[inc].node,No_local[jnc].node,i,j) += J[jnc][i][j];

  }

    if (periodic)
      e2n.swap(e2n_periodic);
  }


  /* free memory for temporary cell arrays */
  free_d3tensor(A);
  free_d3tensor(K);
  free_d3tensor(B);
  free_d3tensor(J);


  /* boundaries and sources */
  /*
   * if (turmod && turbulence_coupling==2)
   *   turb_source();
   */
  if (wall_functions)
    turb_wfuncs(&ls_aztec_coupled->m_B[0]);
  Iboundary();


  std::cout << "IJacobian_P." << std::endl;
}

