
#include "common.h"

void IJacobian_scalar(int iv_scalar)
{
  int inc_min;
  double vol;
  double Prt;
  local_node_struct No_local[4];


  std::cout << "Scalar Jacobian..." << std::endl;


  /*******************/
  /* LOOP OVER CELLS */
  /*******************/
  for (int ic=0; ic<Ncell; ++ic) {

    /* compute cell normals and volume */
    cellgeom(ic,No_local,&vol,&inc_min);

    /* switch to modified connectivity if periodic */
    if ( periodic ) {
      e2n.swap(e2n_periodic);
      for (int inc=0; inc<Nvtcell; ++inc)
        No_local[inc].node = e2n[ic].n[inc];
    }

    /* copy nodal values from global to local structure */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv] = No_W[iv][e2n[ic].n[inc]];

    /* calculate residuals arising from cell */
    double scadiff = nulam/prandtl;
    if ( turmod ) {
      const double nuturb = turb_viscosity(No_local,turmod,Ce_type[ic],vol);
      Prt = 1. - exp(-5.165*nulam/(nuturb+1.e-10));
      Prt = 1./(0.5882+(nuturb/nulam)*(0.228-0.0441*(nuturb/nulam)*Prt));
      Prt = 0.9;
      scadiff += nuturb/Prt;
    }
    scacde(ic,iv_scalar,No_local,vol,inc_min,scaconv,scadiff,0.,0.,1);

    /* add residuals and Jacobian entries to global arrays */
    for (int inc=0; inc<Nvtcell; ++inc) {
      ls_aztec_scalar->B(e2n[ic].n[inc]) += No_local[inc].Res[iv_scalar];
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        ls_aztec_scalar->A(No_local[inc].node,No_local[jnc].node) += No_local[inc].C[jnc];
    }

    if ( periodic )
      e2n.swap(e2n_periodic);
  }


  /**************************/
  /* BOUNDARIES AND SOURCES */
  /**************************/
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    if ( BCgroup[ig].type==IBFIXV || BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBPERE )
      for (int  inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;
        ls_aztec_scalar->zerorow(inu);
        ls_aztec_scalar->A(inu,inu) = 1.;
      }
    else if ( BCgroup[ig].type==IBWALQ )
      for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;
        ls_aztec_scalar->B(inu) += Nobg[ig][inb].modn*BCgroup[ig].invals[0];
      }
  }

  std::cout << "Scalar Jacobian." << std::endl;
}

