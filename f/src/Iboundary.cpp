
/* implicit boundary treatment */

#include <cstdio>
#include "common.h"

void Iboundary()
{
  /* set BC's for pressure-velocity system */
  /* (loop over all boundary groups and all nodes in each group) */
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    switch (BCgroup[ig].type) {

      /* FIXED PRESSURE */
      case IBFIXP :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls_aztec_coupled->zerorow(inu,0);
          ls_aztec_coupled->A(inu,inu,0,0) = -No_vol[inu];
          ls_aztec_coupled->B(inu,0)       = -No_vol[inu]*( No_W[0][inu]
            - 0.66667*(turmod? No_W[iv_turb1][inu]:0.) );

          /* if.option is active then freeze all velocity components except one */
          if (BCgroup[ig].option) {
            for (int iv=1; iv<=Ndim; ++iv)
              if (iv!=BCgroup[ig].option) {
                ls_aztec_coupled->zerorow(inu,iv);
                ls_aztec_coupled->A(inu,inu,iv,iv) = 1.;
              }
          }
        }
        break;

      // FIXED FLUX and WALL (FIXED TEMPERATURE, SPECIFIED HEAT FLUX)
      case IBFIXV : case IBWALL : case IBWALQ :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          for (int iv=1; iv<=Ndim; ++iv) {
            ls_aztec_coupled->zerorow(inu,iv);
            ls_aztec_coupled->A(inu,inu,iv,iv) = 1.;
          }
        }
        break;

      // SYMMETRY PLANE (X, Y or Z DIRN.)
      case IBSYMM :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int iv = BCgroup[ig].option;
          const int inu = Nobg[ig][inb].node;
          ls_aztec_coupled->zerorow(inu,iv);
          ls_aztec_coupled->A(inu,inu,iv,iv) = 1.;
        }
        break;

      // PERIODIC EXIT
      case IBPERE :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          for (int  iv=0; iv<Ncoupled; ++iv) {
            No_W[iv][inu]=No_W[iv][Nobg[ig][inb].twin];
            ls_aztec_coupled->zerorow(inu,iv);
            ls_aztec_coupled->A(inu,inu,iv,iv) = 1.;
          }
        }
        break;

      // GROUP NOT ACTIVE
      case IBNONE : default :
        break;

    }
  }


  /* set temperature BC for inlets and fixed-temperature walls */
  if (temperature && scalar_coupling) {
    double resT_plus  = 0.;
    double resT_minus = 0.;
    for (int ig=1; ig<=Nbcgroup; ++ig) {
      if (BCgroup[ig].type==IBFIXV || BCgroup[ig].type==IBWALL)
      for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;

        if (BCgroup[ig].type==IBWALL) {
          const double t = ls_aztec_coupled->B(inu,iv_temp);
          resT_plus  += (t>0.? t:0.);
          resT_minus += (t>0.? 0.:t);
        }

        ls_aztec_coupled->zerorow(inu,iv_temp);  ls_aztec_coupled->A(inu,inu,iv_temp,iv_temp) = 1.;
      }
    }
    std::cout << "Total fixed-temperature-wall heat fluxes > 0: " << resT_plus  << std::endl
              << "Total fixed-temperature-wall heat fluxes < 0: " << resT_minus << std::endl;
  }


  /* set BC's for turbulence system */
  if (turmod && turbulence_coupling==2) {
    for (int ig=1; ig<=Nbcgroup; ++ig)
      if (BCgroup[ig].type==IBFIXV)
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls_aztec_coupled->zerorow(inu,iv_turb1); ls_aztec_coupled->A(inu,inu,iv_turb1,iv_turb1) = 1.;
          ls_aztec_coupled->zerorow(inu,iv_turb2); ls_aztec_coupled->A(inu,inu,iv_turb2,iv_turb2) = 1.;
        }
      else if (BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ)
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls_aztec_coupled->zerorow(inu,iv_turb1); ls_aztec_coupled->A(inu,inu,iv_turb1,iv_turb1) = 1.;
          ls_aztec_coupled->zerorow(inu,iv_turb2); ls_aztec_coupled->A(inu,inu,iv_turb2,iv_turb2) = 1.;
        }

    /* fix k-e values at wall-function nodes */
    if (wall_functions)
      for (int inwf=0; inwf<(int) WFnodes.size(); ++inwf) {
        const int inu = WFnodes[inwf].node;
        ls_aztec_coupled->zerorow(inu,iv_turb1); ls_aztec_coupled->A(inu,inu,iv_turb1,iv_turb1) = 1.;
        ls_aztec_coupled->zerorow(inu,iv_turb2); ls_aztec_coupled->A(inu,inu,iv_turb2,iv_turb2) = 1.;
      }
  }

}

