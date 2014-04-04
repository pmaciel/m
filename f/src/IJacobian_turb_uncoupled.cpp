
#include "boost/progress.hpp"
#include "common.h"

void IJacobian_turb_uncoupled()
{
  int inc_min;
  double vol;
  local_node_struct No_local[4];


  std::cout << "IJacobian_turb_uncoupled..." << std::endl;


  /*******************/
  /* LOOP OVER CELLS */
  /*******************/
  boost::progress_display pbar(Ncell);
  for (int ic=0; ic<Ncell; ++ic, ++pbar) {

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

    /* calculate turbulent viscosity */
    const double nuturb = turb_viscosity(No_local,turmod,Ce_type[ic],vol);

    /* convection and diffusion of turbulent variables */

    scacde(ic,iv_turb1,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig1,0.,0.,1);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        ls_turb1->A(No_local[inc].node,No_local[jnc].node) += No_local[inc].C[jnc];

    scacde(ic,iv_turb2,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig2,0.,0.,1);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        ls_turb2->A(No_local[inc].node,No_local[jnc].node) += No_local[inc].C[jnc];

    /* add residuals to global arrays */
    for (int inc=0; inc<Nvtcell; ++inc) {
      ls_turb1->B(e2n[ic].n[inc]) += No_local[inc].Res[iv_turb1];
      ls_turb2->B(e2n[ic].n[inc]) += No_local[inc].Res[iv_turb2];
    }

    if (periodic)
      e2n.swap(e2n_periodic);

  }


  /**************************/
  /* BOUNDARIES AND SOURCES */
  /**************************/

  /* sources for turbulence equations (nodal version) */
  turb_source_node();
  turb_wallbc(0,0,ls_turb1,ls_turb2);

  /* main bc groups */
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    if (BCgroup[ig].type==IBFIXV || BCgroup[ig].type==IBPERE)
      for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;
        ls_turb1->zerorow(inu); ls_turb1->A(inu,inu) = 1.;
        ls_turb2->zerorow(inu); ls_turb2->A(inu,inu) = 1.;
      }
  }


  std::cout << "IJacobian_turb_uncoupled." << std::endl;
}

