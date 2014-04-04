
#include <fstream>
#include "boost/progress.hpp"
#include "common.h"

void IJacobian_turb_coupled()
{
  double gradk[3],
         gradw[3];
  double vol;
  int inc_min;
  local_node_struct No_local[4];

  LS *ls = (turbulence_coupling==1? ls_turb    :
           (turbulence_coupling==2? ls_coupled : NULL ));
  const int iv1 = turbulence_coupling==2? iv_turb1 : 0;
  const int iv2 = turbulence_coupling==2? iv_turb2 : 1;


  std::cout << "IJacobian_turb_coupled..." << std::endl;


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
        No_local[inc].W[iv] = No_W[iv][e2n[ic].n[inc]];

    /* calculate turbulent viscosity */
    const double nuturb = turb_viscosity(No_local,turmod,Ce_type[ic],vol);

    /* convection and diffusion of turbulent variables */
    if (turmod==ITKWBS || turmod==ITKWSS) {
      double wdist = 0.,
             turb1 = 0.,
             turb2 = 0.;
      for (int inc=0; inc<Nvtcell; ++inc) {
        turb1 += No_local[inc].W[iv_turb1];
        turb2 += No_local[inc].W[iv_turb2];
        wdist += No_wd[No_local[inc].node];
      }
      turb1 /= dNvtcell;
      turb2 /= dNvtcell;
      wdist /= dNvtcell;

      for (int id=0; id<Ndim; ++id) {
        gradk[id] = 0.;
        gradw[id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc) {
          gradk[id] += No_W[iv_turb1][No_local[inc].node]*No_local[inc].norm[id];
          gradw[id] += No_W[iv_turb2][No_local[inc].node]*No_local[inc].norm[id];
        }
        gradk[id] /= dNvtfce*vol;
        gradw[id] /= dNvtfce*vol;
      }

      double dkdw = 0.;
      for (int id=0; id<Ndim; ++id)
        dkdw += gradk[id]*gradw[id];

      const double F1 = iter<=turb_iterinit? 1.
                                       : F1_function(turb1,turb2,nulam,wdist,dkdw);

      Sig1 = 1./(F1*Sigk1 + (1.-F1)*Sigk2);
      Sig2 = 1./(F1*Sigw1 + (1.-F1)*Sigw2);
    }

    scacde(ic,iv_turb1,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig1,0.,0.,1);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        (ls->A)(No_local[inc].node,No_local[jnc].node,iv1,iv1) += No_local[inc].C[jnc];

    scacde(ic,iv_turb2,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig2,0.,0.,1);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        (ls->A)(No_local[inc].node,No_local[jnc].node,iv2,iv2) += No_local[inc].C[jnc];

    /* add residuals to global arrays */
    for (int inc=0; inc<Nvtcell; ++inc) {
      (ls->B)(e2n[ic].n[inc],iv1) += No_local[inc].Res[iv_turb1];
      (ls->B)(e2n[ic].n[inc],iv2) += No_local[inc].Res[iv_turb2];
    }

    if (periodic)
      e2n.swap(e2n_periodic);

  }

  /**************************/
  /* BOUNDARIES AND SOURCES */
  /**************************/

  /* sources for turbulence equations (nodal version) */
  turb_source_node();
  turb_wallbc(iv1,iv2,ls,ls);

  /* main bc groups */
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    if (BCgroup[ig].type==IBFIXV || BCgroup[ig].type==IBPERE)
      for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;
        ls->zerorow(inu,iv1);  (ls->A)(inu,inu,iv1,iv1) = 1.;
        ls->zerorow(inu,iv2);  (ls->A)(inu,inu,iv2,iv2) = 1.;
      }
  }


  std::cout << "IJacobian_turb_coupled." << std::endl;
}

