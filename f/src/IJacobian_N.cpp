
/*
 * compute the Jacobian matrix of the residual by finite difference
 * approximations as well as the residual itself. Final update for boundary
 * elements are computed elsewhere
 */

#include <fstream>
#include "common.h"

void IJacobian_N()
{
  int inc_min;
  double vol;
  double betasq;
  double lwfac;

  int scaconv_save;
  local_node_struct No_local[4];
  local_node_struct No_eps[4];


  std::cout << "IJacobian_N..." << std::endl;


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

    for (int inc=0; inc<Nvtcell; ++inc) {
      No_eps[inc].node=No_local[inc].node;
      for (int id=0; id<Ndim; ++id)
        No_eps[inc].norm[id]=No_local[inc].norm[id];
      No_eps[inc].norm2=No_local[inc].norm2;
    }

    /* Switch to modified connectivity if periodic */
    if (periodic)
      e2n.swap(e2n_periodic);

    /* Compute unperturbed residual */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv] = No_W[iv][e2n[ic].n[inc]];

    /* Calculate convective contributions */
    convection(No_local,vol,inc_min,Ce_type[ic],&betasq,&lwfac,k,sbuoy,A,K,B,Rc);

    /* Calculate scalar equation contributions */
    if (temperature && scalar_coupling)
      scacde(ic,iv_temp,No_local,vol,inc_min,scaconv,nulam/prandtl,0.,0.,0);

    /*
     * if (turmod && turbulence_coupling==2) {
     *   nuturb = turb_viscosity(No_local,turmod,Ce_type[ic],vol);
     *   turb_source_cell(No_local,vol,turmod,&tsrce1,&tsrce2,&tcoeff1,&tcoeff2);
     *   scacde(ic,iv_turb1,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig1,tsrce1,0.,0);
     *   scacde(ic,iv_turb2,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig2,tsrce2,0.,0);
     * }
     */

    /* Calculate viscous contributions */
    viscous(No_local,vol,Ce_type[ic]);

    /*
     * for (int inc=0; inc<Nvtcell; ++inc)
     *   for (int iv=0; iv<Neqns; ++iv)
     *     No_eps[inc].Res[iv]=No_local[inc].Res[iv];
     */

    /* Add unperturbed residuals and Dt's into global No structure */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int inu = e2n[ic].n[inc];
      for (int iv=0; iv<Ncoupled; ++iv)
        ls_aztec_coupled->B(inu,iv) += No_local[inc].Res[iv];
      No_dt[inu] += No_local[inc].Dt;
    }

    /*
     * for (int inc=0; inc<Nvtcell; ++inc)
     *   for (int iv=0; iv<Neqns; ++iv)
     *     No_local[inc].Res[iv]=No_eps[inc].Res[iv];
     */

    scaconv_save = scaconv;
    if (temperature && scalar_coupling && scaconv==ISSPSI) {
      scaconv = ISSNSC;
      scacde(ic,iv_temp,No_local,vol,inc_min,scaconv,nulam/prandtl,0.,0.,0);
    }

    /* loop over nodes within cell */
    int node[4];
    int jnc[4];
    double U[6];
    for (int inc=0; inc<Nvtcell; ++inc) {
      /* loop over variables */
      for (int iv=0; iv<Ncoupled; ++iv) {

        /* Find node numbers of other nodes in cell */
        for (int i=0; i<Nvtcell; ++i) {
          jnc[i] = (inc+i)%Nvtcell;
          node[i] = e2n[ic].n[jnc[i]];
        }

        /* Perturb value of variable at inc */
        copy(No_local[inc].W,U,Nsys);
        const double tmp = U[iv];
        const double typ = .001;
        U[iv] += (U[iv]<0.? -1.:1.) * newton_eps * std::max< double >(std::abs(U[iv]),typ);;
        const double perturb = 1./(U[iv]-tmp);
        copy(U,No_eps[inc].W,Nsys);

        /* Set remaining nodal values in No_eps to unperturbed state */
        for (int i=1; i<Nvtcell; ++i)
          copy(No_local[jnc[i]].W,No_eps[jnc[i]].W,Nsys);

        /* Compute perturbed residual */
        convection(No_eps,vol,inc_min,Ce_type[ic],&betasq,&lwfac,k,sbuoy,A,K,B,Rc);

        /* Calculate viscous contributions */
        viscous(No_eps,vol,Ce_type[ic]);

        /* Compute perturbed residual for scalar equation */
        if (temperature && scalar_coupling)
          scacde(ic,iv_temp,No_eps,vol,inc_min,scaconv,nulam/prandtl,0.,0.,0);

        /* Convection and diffusion of turbulent variables */
        /*
         * if (turmod && turbulence_coupling==2) {
         *   nuturb = turb_viscosity(No_eps,turmod,Ce_type[ic],vol);
         *   turb_source_cell(No_eps,vol,turmod,&tsrce1,&tsrce2,&tcoeff1,&tcoeff2);
         *   scacde(ic,iv_turb1,No_eps,vol,inc_min,scaconv,nulam+nuturb/Sig1,tsrce1,0.,0);
         *   scacde(ic,iv_turb2,No_eps,vol,inc_min,scaconv,nulam+nuturb/Sig2,tsrce2,0.,0);
         * }
         */

        /* Calculate own-node (diagonal) Jacobian contributions */
        for (int jv=0; jv<Ncoupled; ++jv)
          J[inc][jv][iv] = (No_eps[inc].Res[jv]-No_local[inc].Res[jv])*perturb;

        /* Calculate neighbour (off-diagonal) Jacobian contributions */
        for (int i=1; i<Nvtcell; ++i)
          for (int jv=0; jv<Ncoupled; ++jv)
            J[jnc[i]][jv][iv]=(No_eps[jnc[i]].Res[jv]-No_local[jnc[i]].Res[jv])*perturb;

        if (turmod && turbulence_coupling==2)
          for (int knc=0; knc<Nvtcell; ++knc) {
            J[knc][iv_turb1][iv_turb2] = 0.;
            J[knc][iv_turb2][iv_turb1] = 0.;
          }

      /* loop over variables */
      }

      /* Add entries in column inc to Jacobian */
      for (int knc=0; knc<Nvtcell; ++knc)
        for (int i=0; i<Ncoupled; ++i)
          for (int j=0; j<Ncoupled; ++j)
            ls_aztec_coupled->A(e2n[ic].n[knc],e2n[ic].n[inc],i,j) += J[knc][i][j];

    /* loop over nodes within cell */
    }

    scaconv = scaconv_save;
    if (periodic)
      e2n.swap(e2n_periodic);

  /* loop over cells */
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


  std::cout << "IJacobian_N." << std::endl;
}

