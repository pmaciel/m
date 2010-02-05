
/* subroutines for turbulence models */

#include "common.h"

void turb_source_node()
{
  int inc_min;
  double vol;
  local_node_struct No_loc[4];
  double gradk[3];
  double gradw[3];
  double gradv[4][3];

  std::vector< double > No_dkdw;
  if (turmod==ITKWSS || turmod==ITKWBS || turmod==ITKWPD)
    No_dkdw.assign(Nnode,0.);

  No_nuturb.assign(Nnode,0.);

  /*** Calculate nodal G function values ***/
  for (int ic=0; ic<Ncell; ++ic) {
    cellgeom(ic,No_loc,&vol,&inc_min);

    /* use periodic connectivity if appropriate */
    if (periodic)
      for (int inc=0; inc<Nvtcell; ++inc)
        No_loc[inc].node = e2n_periodic[ic].n[inc];

    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Nsys; ++iv)
        No_loc[inc].W[iv]=No_W[iv][No_loc[inc].node];

    /* Calculate turbulent viscosity */
    const double nuturb = turb_viscosity(No_loc,turmod,Ce_type[ic],vol);

    for (int inc=0; inc<Nvtcell; ++inc)
      No_nuturb[No_loc[inc].node] += nuturb*vol/dNvtcell;

    /* Average turbulence values */
    double kturb = 0.;
    double turb2 = 0.;
    for (int inc=0; inc<Nvtcell; ++inc) {
      kturb += No_loc[inc].W[iv_turb1];
      turb2 += No_loc[inc].W[iv_turb2];
    }
    kturb /= dNvtcell;
    turb2 /= dNvtcell;

    /* Calculate velocity gradient and G function */
    for (int iv=1; iv<=Ndim; ++iv) {
      for (int id=0; id<Ndim; ++id) {
        gradv[iv][id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc)
          gradv[iv][id] += No_loc[inc].W[iv]*No_loc[inc].norm[id];
        gradv[iv][id] = gradv[iv][id]/(dNvtfce*vol);
      }
    }

    const double Gfunc = Gfunction(Ndim,gradv);

    /* Average wall distance for cell */
    double dwall = 0.;
    if (walldist) {
      for (int inc=0; inc<Nvtcell; ++inc)
        dwall += No_wd[No_loc[inc].node];
      dwall /= dNvtcell;
    }

    /*** Calculate grad(k).grad(w) for BSL, SST and PDH models ***/
    double dkdw = 0.;
    if (turmod==ITKWSS || turmod==ITKWBS || turmod==ITKWPD) {
      for (int id=0; id<Ndim; ++id) {
        gradk[id] = 0.;
        gradw[id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc) {
          gradk[id] += No_W[iv_turb1][No_loc[inc].node]*No_loc[inc].norm[id];
          gradw[id] += No_W[iv_turb2][No_loc[inc].node]*No_loc[inc].norm[id];
        }
        gradk[id] /= dNvtfce*vol;
        gradw[id] /= dNvtfce*vol;
      }

      dkdw = 0.;
      for (int id=0; id<Ndim; ++id)
        dkdw += gradk[id]*gradw[id];

      for (int inc=0; inc<Nvtcell; ++inc)
        No_dkdw[No_loc[inc].node] += dkdw*vol/dNvtcell;
    }

    double v2turb = 0.;
    double SP1 = 0.;
    double SP2 = 0.;
    turb_source_P(turmod,kturb,turb2,nuturb,nulam,dwall,Gfunc,dkdw,&SP1,&SP2,v2turb);

    /* Add cell-based sources and Jacobian entries */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int inu = No_loc[inc].node;
      ls_turb->B(inu,0) += SP1*vol/dNvtcell;
      ls_turb->B(inu,1) += SP2*vol/dNvtcell;
    }

  }

  //////////////////////////////

  for (int inu=0; inu<Nnode; ++inu)
    No_nuturb[inu] /= No_vol[inu];

  if (turmod==ITKWSS || turmod==ITKWBS || turmod==ITKWPD)
    for (int inu=0; inu<Nnode; ++inu)
      No_dkdw[inu] /= No_vol[inu];

  /* Add nodal sources and Jacobian entries */
  for (int inu=0; inu<Nnode; ++inu) {
    double dwall   = 0.;
    double lenturb = 0.;
    if (walldist) {
      dwall = No_wd[inu];
      if (iter<=turb_iterinit)
        lenturb = No_lenturb[inu];
    }

    const double dkdw = (turmod==ITKWSS || turmod==ITKWBS || turmod==ITKWPD?
      No_dkdw[inu] : 0.);

    double SD1  = 0., SD2  = 0.,
           JD1  = 0., JD12 = 0.,
           JD21 = 0., JD2  = 0.;
    turb_source_D(turmod,No_W[iv_turb1][inu],No_W[iv_turb2][inu],nulam,dwall,lenturb,dkdw,&SD1,&SD2,&JD1,&JD2,&JD12,&JD21);

    No_dissipation[inu] = -SD1;

    ls_turb->B(inu,0) += No_vol[inu]*SD1;
    ls_turb->B(inu,1) += No_vol[inu]*SD2;

    ls_turb->A(inu,inu,0,0) += No_vol[inu]*JD1;
    ls_turb->A(inu,inu,1,1) += No_vol[inu]*JD2;
    ls_turb->A(inu,inu,0,1) += No_vol[inu]*JD12;

  }
}

