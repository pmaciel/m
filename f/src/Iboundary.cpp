
/* implicit boundary treatment */

#include <cstdio>
#include "common.h"

void Iboundary()
{
  std::cout << "Iboundary: set BC's for pressure-velocity system..." << std::endl;

  /* (loop over all boundary groups and all nodes in each group) */
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    switch (BCgroup[ig].type) {

      /* FIXED PRESSURE */
      case IBFIXP :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls_coupled->zerorow(inu,0);
          ls_coupled->A(inu,inu,0,0) = -No_vol[inu];
          ls_coupled->B(inu,0)       = -No_vol[inu]*( No_W[0][inu]
            - 0.66667*(turmod? No_W[iv_turb1][inu]:0.) );

          /* if.option is active then freeze all velocity components except one */
          if (BCgroup[ig].option) {
            for (int iv=1; iv<=Ndim; ++iv)
              if (iv!=BCgroup[ig].option) {
                ls_coupled->zerorow(inu,iv);
                ls_coupled->A(inu,inu,iv,iv) = 1.;
              }
          }
        }
        break;

      // FIXED FLUX and WALL (FIXED TEMPERATURE, SPECIFIED HEAT FLUX)
      case IBFIXV : case IBWALL : case IBWALQ :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          for (int iv=1; iv<=Ndim; ++iv) {
            ls_coupled->zerorow(inu,iv);
            ls_coupled->A(inu,inu,iv,iv) = 1.;
          }
        }
        break;

      // SYMMETRY PLANE (X, Y or Z DIRN.)
      case IBSYMM :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int iv = BCgroup[ig].option;
          const int inu = Nobg[ig][inb].node;
          ls_coupled->zerorow(inu,iv);
          ls_coupled->A(inu,inu,iv,iv) = 1.;
        }
        break;

      // PERIODIC EXIT
      case IBPERE :
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          for (int  iv=0; iv<Ncoupled; ++iv) {
            No_W[iv][inu]=No_W[iv][Nobg[ig][inb].twin];
            ls_coupled->zerorow(inu,iv);
            ls_coupled->A(inu,inu,iv,iv) = 1.;
          }
        }
        break;

      // GROUP NOT ACTIVE
      case IBNONE : default :
        break;

    }
  }
  std::cout << "Iboundary: set BC's for pressure-velocity system." << std::endl;


  if (temperature && scalar_coupling) {
    std::cout << "Iboundary: set temperature BC for inlets and fixed-temperature walls..." << std::endl;
    double resT_plus  = 0.;
    double resT_minus = 0.;
    for (int ig=1; ig<=Nbcgroup; ++ig) {
      if (BCgroup[ig].type==IBFIXV || BCgroup[ig].type==IBWALL)
      for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;

        if (BCgroup[ig].type==IBWALL) {
          const double t = ls_coupled->B(inu,iv_temp);
          resT_plus  += (t>0.? t:0.);
          resT_minus += (t>0.? 0.:t);
        }

        ls_coupled->zerorow(inu,iv_temp);  ls_coupled->A(inu,inu,iv_temp,iv_temp) = 1.;
      }
    }
    std::cout << "total fixed-temperature-wall heat fluxes > 0: " << resT_plus  << std::endl
              << "total fixed-temperature-wall heat fluxes < 0: " << resT_minus << std::endl;
    std::cout << "Iboundary: set temperature BC for inlets and fixed-temperature walls." << std::endl;
  }


  if (turmod && turbulence_coupling==2) {
    std::cout << "Iboundary: set BC's for turbulence system..." << std::endl;

    for (int ig=1; ig<=Nbcgroup; ++ig)
      if (BCgroup[ig].type==IBFIXV)
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls_coupled->zerorow(inu,iv_turb1); ls_coupled->A(inu,inu,iv_turb1,iv_turb1) = 1.;
          ls_coupled->zerorow(inu,iv_turb2); ls_coupled->A(inu,inu,iv_turb2,iv_turb2) = 1.;
        }
      else if (BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ)
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls_coupled->zerorow(inu,iv_turb1); ls_coupled->A(inu,inu,iv_turb1,iv_turb1) = 1.;
          ls_coupled->zerorow(inu,iv_turb2); ls_coupled->A(inu,inu,iv_turb2,iv_turb2) = 1.;
        }

    if (wall_functions) {
      std::cout << "Iboundary: fix k-e values at wall-function nodes..." << std::endl;
      for (int inwf=0; inwf<(int) WFnodes.size(); ++inwf) {
        const int inu = WFnodes[inwf].node;
        ls_coupled->zerorow(inu,iv_turb1); ls_coupled->A(inu,inu,iv_turb1,iv_turb1) = 1.;
        ls_coupled->zerorow(inu,iv_turb2); ls_coupled->A(inu,inu,iv_turb2,iv_turb2) = 1.;
      }
      std::cout << "Iboundary: fix k-e values at wall-function nodes." << std::endl;
    }

    std::cout << "Iboundary: set BC's for turbulence system." << std::endl;
  }

}

