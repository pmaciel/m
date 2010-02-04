
/* implicit boundary treatment */

#include "common.h"

void massflux()
{
  local_node_struct No_local[4];
  double vol;
  int inc_min;

  for (int ig=1; ig<=Nbcgroup; ++ig)
    BCgroup[ig].qflux=0.;

  for (int ifc=0; ifc<Nbface; ++ifc) {
    const int ic=Fab_cell[ifc];

    std::vector< double > vavge(Ndim,0.);
    for (int inf=1; inf<=Nvtfce; ++inf) {
      const int inc=(inf+Fab_inc[ifc])%Nvtcell;
      for (int id=0; id<Ndim; ++id)
        vavge[id] += No_W[id+1][e2n[ic].n[inc]];
    }

    cellgeom(ic,No_local,&vol,&inc_min);

    double vdotn = 0.;
    for (int id=0; id<Ndim; ++id)
      vdotn += (vavge[id]/dNvtfce)*No_local[Fab_inc[ifc]].norm[id];
    BCgroup[Fab_group[ifc]].qflux += vdotn;
  }

  // calculate total inlet and outlet mass fluxes
  Qin  = 0.;
  Qout = 0.;
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    if ( BCgroup[ig].qflux>=0. )
      Qin += BCgroup[ig].qflux;
    else
      Qout += BCgroup[ig].qflux;
  }
}

