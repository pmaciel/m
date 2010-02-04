
/*
 * compute the Jacobian matrix of the residual by finite difference
 * approximations as well as the residual itself. Final update for boundary
 * elements are computed elsewhere
 */

#include "common.h"

void IJacobian_A()
{
  int inc_min;
  double vol;
  double betasq;
  double lwfac;

  int scaconv_save;
  double U[6];
  local_node_struct No_local[4];
  local_node_struct No_eps[4];


  std::cout << "IJacobian_A..." << std::endl;


  /* allocate memory for temporary cell arrays */
  std::vector< double > k    (Nvtcell,    0.);
  std::vector< double > sbuoy(Ncoupled+1, 0.);
  std::vector< double > Rc    (Ncoupled,  0.);
  std::vector< double > Rc_org(Ncoupled,  0.);
  double ***A     = d3tensor(Ndim,   Ndim+1,  Ndim+1);
  double ***K     = d3tensor(Nvtcell,Ndim+1,  Ndim+1);
  double ***B     = d3tensor(Nvtcell,Ndim+1,  Ndim+1);
  double ***J     = d3tensor(Nvtcell,Ncoupled,Ncoupled);
  double ***B_org = d3tensor(Nvtcell,Ndim+1,  Ndim+1);
  double ***dPdU  = d3tensor(Nvtcell,Ncoupled,Ncoupled);


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

    /* compute unperturbed residual */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv]=No_W[iv][e2n[ic].n[inc]];

    /* calculate convective contributions */
    convection(No_local,vol,inc_min,Ce_type[ic],&betasq,&lwfac,k,sbuoy,A,K,B,Rc);

    /* calculate viscous contributions */
    viscous(No_local,vol,Ce_type[ic]);

    /* add unperturbed residuals and Dt's into global No structure */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int inu = e2n[ic].n[inc];
      for (int iv=0; iv<Ncoupled; ++iv)
        ls_aztec_coupled->B(inu,iv) += No_local[inc].Res[iv];
      No_dt[inu] += No_local[inc].Dt;
    }

    scaconv_save=scaconv;
    if ( (temperature || turmod) && scaconv==6 ) {
      scaconv=2;
      convection(No_local,vol,inc_min,Ce_type[ic],&betasq,&lwfac,k,sbuoy,A,K,B,Rc);
    }

    Rc_org = Rc;
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Ncoupled; ++iv)
        for (int jv=0; jv<Ncoupled; ++jv)
          B_org[inc][iv][jv] = B[inc][iv][jv];

    /* loop over nodes within cell */
    for (int inc=0; inc<Nvtcell; ++inc) {
      /* loop over variables */
      for (int iv=0; iv<Ncoupled; ++iv) {

        /* perturb value of variable at inc */
        copy(No_local[inc].W,U,Nsys);
        const double tmp = U[iv];
        const double typ = .001;
        U[iv] += (U[iv]<0.? -1.:1.) * newton_eps * std::max< double >(std::abs(U[iv]),typ);
        const double perturb = 1./(U[iv]-tmp);
        copy(U,No_eps[inc].W,Nsys);

        /* set remaining nodal values in No_eps to unperturbed state */
        for (int i=1; i<Nvtcell; ++i) {
          const int jnc=(inc+i)%Nvtcell;
          copy(No_local[jnc].W,No_eps[jnc].W,Nsys);
        }

        /* compute perturbed residual */
        convection(No_eps,vol,inc_min,Ce_type[ic],&betasq,&lwfac,k,sbuoy,A,K,B,Rc);

        /* calculate own-node (diagonal) Jacobian contributions */
        for (int jv=0; jv<Ncoupled; ++jv)
          dPdU[inc][jv][iv] = (Rc[jv]-Rc_org[jv])*perturb;

      /* loop over variables */
      }
    /* loop over nodes within cell */
    }

    for (int inc=0; inc<Nvtcell; ++inc) {
      for (int jnc=0; jnc<Nvtcell; ++jnc) {
        for (int iv=0; iv<Ncoupled; ++iv)
          for (int jv=0; jv<Ncoupled; ++jv)
            J[jnc][iv][jv] = 0.;

        for (int iv=0; iv<Ncoupled; ++iv)
          for (int jv=0; jv<Ncoupled; ++jv)
            for (int ij=0; ij<Ncoupled; ++ij)
              J[jnc][iv][jv] -= B_org[inc][iv][ij]*dPdU[jnc][ij][jv];

        for (int i=0; i<Ncoupled; ++i)
          for (int j=0; j<Ncoupled; ++j)
            ls_aztec_coupled->A(No_local[inc].node,No_local[jnc].node,i,j) += J[jnc][i][j];
      }
    }

    scaconv = scaconv_save;

  /* loop over cells */
  }


  /* free memory for temporary cell arrays */
  free_d3tensor(A);
  free_d3tensor(K);
  free_d3tensor(B);
  free_d3tensor(J);
  free_d3tensor(B_org);
  free_d3tensor(dPdU);


  /* boundaries and sources */
  if (wall_functions)
    turb_wfuncs(&(ls_aztec_coupled->m_B)[0]);
  Iboundary();


  std::cout << "IJacobian_A." << std::endl;
}

