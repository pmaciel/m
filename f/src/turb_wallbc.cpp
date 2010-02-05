
/* subroutines for turbulence models */

#include "common.h"

void turb_wallbc(int iv1, int iv2, LS *ls1, LS *ls2)
{
  /* Set k residuals on wall nodes to zero for all models */
  for (int ig=1; ig<=Nbcgroup; ++ig)
    if ( BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ )
      for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
        const int inu = Nobg[ig][inb].node;
        ls1->zerorow(inu,iv1); (ls1->A)(inu,inu,iv1,iv1) = 1.;
      }

  /* Main BC groups */
  for (int ig=1; ig<=Nbcgroup; ++ig) {
    if ( BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ ) {

      if ( turmod==ITKEHR || turmod==ITKEHG ) {
        /* high-Re k-e models (with wall functions) */
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          ls1->zerorow(inu,iv1); (ls1->A)(inu,inu,iv1,iv1) = 1.;
          ls2->zerorow(inu,iv2); (ls2->A)(inu,inu,iv2,iv2) = 1.;
        }
      }

      else if ( turmod==ITKENA || turmod==ITKELB || turmod==ITKELS ) {
        /* low-Re k-e models */
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu  = Nobg[ig][inb].node;
          const int intw = Nobg[ig][inb].twin;
          const double epswall = 2.*nulam*No_W[iv_turb1][Nobg[ig][inb].twin]/(Nobg[ig][inb].dist*Nobg[ig][inb].dist);
          ls1->zerorow(inu,iv1); (ls1->A)(inu,inu,iv1,iv1) = 1.;
          ls2->zerorow(inu,iv2); (ls2->B)(inu,iv2)         = No_vol[inu]*(epswall-No_W[iv_turb2][inu]);
          (ls2->A)(inu,inu, iv2,iv2)   = -No_vol[inu];
          if (iv2<(int) ls1->Nb)  // if there is room for the coupling term
            (ls1->A)(inu,intw,iv2,iv1) = No_vol[inu]*2.*nulam/(Nobg[ig][inb].dist*Nobg[ig][inb].dist);
        }
      }

      else if ( turmod==ITKWHR || turmod==ITKWBS || turmod==ITKWSS || turmod==ITKWLR || turmod==ITKWPD ) {
        /* high and low-Re k-omega models */
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].node;
          No_W[iv_turb2][inu] = 60.0*nulam/(Cw2*Nobg[ig][inb].dist*Nobg[ig][inb].dist);
          ls1->zerorow(inu,iv1); (ls1->A)(inu,inu,iv1,iv1) = 1.;
          ls2->zerorow(inu,iv2); (ls2->A)(inu,inu,iv2,iv2) = 1.;
        }
      }

    }
  }

  if (!wall_functions) {

    double vtan[3];
    double norm[3];
    double ypmax = 0.,    ypmin = 1.e20;
    double utmax,         utmin;
    int    in_ypmax = -1, in_ypmin = -1;
    int Nwnode   = 0;
    int nyp_near = 0;
    int nyp_far  = 0;
    for (int ig=1; ig<=Nbcgroup; ++ig) {
      if ( BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ )
        for (int inb=1; inb<=BCgroup[ig].nnode; ++inb) {
          const int inu = Nobg[ig][inb].twin;

          const double wdist = Nobg[ig][inb].dist;
          const double modn  = Nobg[ig][inb].modn;
          for (int id=0; id<Ndim; ++id)
            norm[id] = Nobg[ig][inb].n[id]/modn;

          double vdotn = 0.;
          for (int id=0; id<Ndim; ++id)
            vdotn += No_W[id+1][inu]*norm[id];

          double vt = 0.;
          for (int id=0; id<Ndim; ++id) {
            vtan[id] = No_W[id+1][inu] - vdotn*norm[id] - No_W[id+1][Nobg[ig][inb].node];
            vt += vtan[id]*vtan[id];
          }
          vt = sqrt(vt);

          const double utau2 = nulam*vt/wdist;
          const double yplus = wdist*sqrt(utau2)/nulam;

          /* y+ statistics */
          Nwnode++;
          if (yplus>ypmax) {
            ypmax    = yplus;
            utmax    = sqrt(utau2);
            in_ypmax = inu;
          }
          else if (yplus<ypmin && yplus>1.e-8) {
            ypmin    = yplus;
            utmin    = sqrt(utau2);
            in_ypmin = inu;
          }
          nyp_far  += yplus>2.? 1:0;
          nyp_near += yplus<1.? 1:0;
        }
    }

    const double dNwnode = (double) Nwnode;
    std::cout << "turb_wallbc: near-wall nodes report..." << std::endl
              << "Maximum y+=" << ypmax << " at n=" << in_ypmax << " with utau=" << utmax << std::endl
              << "Minimum y+=" << ypmin << " at n=" << in_ypmin << " with utau=" << utmin << std::endl
              << "N. nodes with y+>2.0: " << nyp_far  << " (" << 100.*(double)nyp_far /dNwnode << "%)" << std::endl
              << "N. nodes with y+<1.0: " << nyp_near << " (" << 100.*(double)nyp_near/dNwnode << "%)" << std::endl
              << "turb_wallbc: near-wall nodes report." << std::endl;

  }
  else {

    /*
     * Fix k and epsilon values at wall-function nodes
     * N.B. This is just to zero the residuals and links,
     * the actual values of k and e are set in turb_wfuncs()
     * which is called from the main Jacobian routines as
     * this must be carried out at the same time as for the momentum
     */
    for (int inwf=0; inwf<(int) WFnodes.size(); ++inwf) {
      const int inu = WFnodes[inwf].node;
      ls1->zerorow(inu,iv1); (ls1->A)(inu,inu,iv1,iv1) = 1.;
      ls2->zerorow(inu,iv2); (ls2->A)(inu,inu,iv2,iv2) = 1.;
    }
  }

}

