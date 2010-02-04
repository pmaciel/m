
#include "common.h"

int main(int argc, char **argv)
{
  using namespace std;

  GetPot o(argc,argv);
  if (o.search(2,"--help", "-h")) {
    cout << endl
         << "Usage (" << o[0] << "):" << endl
         << "  --help, -h: this help" << endl
         << "  --case, -c: testcase file (default: \"testcases.c.xml\")" << endl
         << endl;
    return 0;
  }


  // setup from testcase file
  readstart(o.follow("testcases.c.xml",2,"--case","-c"));

  // - read grid/solution and set connectivity/grid-geometry data
  // - write out read data
  readsoltp(file_input,restart);
  writesoltp(file_output);


  // basic procedure on each timestep:
  // - compute residual in cell
  // - distribute residual in cell to nodes
  // - update nodal values using distributed residuals

  cout << "starting main update loop..." << endl;

  /* start of main loop */
  do {
    iter++;

    /* find maximum absolute velocity in field */
    Umax = 0.;
    C2   = 1.;
    for (int inu=0; inu<Nnode; inu++) {
      double Uabs=0.;
      for (int id=0; id<Ndim; id++)
        Uabs += No_W[id+1][inu]*No_W[id+1][inu];
      Umax = std::max< double >(Uabs,Umax);
    }
    Umax = sqrt(Umax);
    if (Umax<epsilon)
      Umax=1.;
    C2 = 0.25*Umax*Umax;


    /* non-linear update loop */

    /* 1) compute jacobian and residual */

    if (Jacobian==0 && logL2[iverr]<newton_thresh && iter>2) {
      cout << "*** Switching to Newton ***" << endl;
      Jacobian = 2;
      CFL *= 100.;
      linrlx = 1.;
/*
      ls_aztec_coupled->set_params(0)  /= 100.;  // AZ_tol
      ls_aztec_coupled->set_options(6) *= 2;     // AZ_max_iter
*/
      if (temperature && scalar_coupling)
        linrlx_scalar = 1.;
      if (turmod && turbulence_coupling)
        linrlx_turb = 1.;
    }

    /* zero nodal turbulent viscosity store */
    if (turmod)
      No_nuturb.assign(Nnode,0.);

    /* calculate bulk temperature for boussinesq buoyancy case */
    if (temperature && buoyancy) {
      To  = 0.;
      for (int inu=0; inu<Nnode; inu++)
        To += No_W[iv_temp][inu]*No_vol[inu];
      To /= Vtot;
      double sum = 0.;
      for (int inu=0; inu<Nnode; inu++)
        sum += No_vol[inu]/(287.0*No_W[iv_temp][inu]);
      if (iter<=1)
        masstot = 1.e5*sum;
      Po = masstot/sum;
      cout << "Total mass=" << masstot << " To=" << To << " Po=" << Po << endl;
    }


    cout << "Computing nodal residuals and steady Jacobian..." << endl;

    /* zero jacobian, residual and vector */
    ls_aztec_coupled->reset();

    No_dt.assign(Nnode,0.);

    /* build jacobian and residual vector */
    if      (Jacobian==0)  IJacobian_P();
    else if (Jacobian==1)  IJacobian_A();
    else if (Jacobian==2)  IJacobian_N();

    cout << "Computing nodal residuals and steady Jacobian." << endl;


    /* calculate residuals for coupled system */
    for (int iv=0; iv<Ncoupled; iv++)
      rescalc(iv,iv,&(ls_aztec_coupled->m_B)[0],Ncoupled);


    for (int n=0; n<Nnode; n++)
      for (int iv=0; iv<Ncoupled; iv++)
        ls_aztec_coupled->B(n,iv) *= linrlx;


    if (dtrelax) {
      cout << "adding time-step (backward-euler) to coupled Jacobian..." << endl;

      if (dtrelax==1) {
        /* global time-step */
        double tmax     = 0.;
        double kplusmax = 0.;
        for (int inu=0; inu<Nnode; inu++) {
          kplusmax = std::max< double >(0.,No_dt[inu]);
          tmax     = std::max< double >(tmax,kplusmax/No_vol[inu]);
        }
        tmax = 1./tmax;

        for (int inu=0; inu<Nnode; inu++)
          No_dt[inu] = -No_vol[inu]/tmax;
      }
      else if (dtrelax==2) {
        /* local time step */
        for (int inu=0; inu<Nnode; inu++)
          No_dt[inu] = -No_dt[inu];
      }

      for (int inu=0; inu<Nnode; ++inu)
        for (int iv=1; iv<Ncoupled; ++iv)
          ls_aztec_coupled->A(inu,inu,iv,iv) += No_dt[inu]/CFL;

      if (logL2[iverr]<newton_thresh)
        CFL *= 1.5;

      cout << "adding time-step (backward-euler) to coupled Jacobian." << endl;
    }

    cout << "solve coupled system..." << endl;
    ls_aztec_coupled->solve();
    cout << "solve coupled system." << endl;


    /*
     * Carry out decoupled turbulence and temperature updates using
     * old velocities (before update). This is slower but more stable
     * and handy for fully-decoupled temperature solutions
     *
     * if (temperature && !scalar_coupling)
     *   Update_temperature();
     */


    cout << "update solution..." << endl;

    /* pressure-velocity system */
    for (int n=0; n<Nnode; ++n)
      for (int iv=0; iv<=Ndim; ++iv)
        No_W[iv][n] -= /*linrlx**/ls_aztec_coupled->X(n,iv);

    /* temperature */
    if (temperature && scalar_coupling) {
      for (int n=0; n<Nnode; ++n)
        No_W[iv_temp][n] -= linrlx_scalar*ls_aztec_coupled->X(n,iv_temp);
    }

    /* turbulence */
    if (turmod && turbulence_coupling==2) {
      int nkneg = 0;
      int neneg = 0;
      for (int n=0; n<Nnode; ++n) {
        const double k = No_W[iv_turb1][n] - linrlx_turb*ls_aztec_coupled->X(n,iv_turb1);
        const double e = No_W[iv_turb2][n] - linrlx_turb*ls_aztec_coupled->X(n,iv_turb2);
        if (k<0.) ++nkneg; else No_W[iv_turb1][n] = k;
        if (e<0.) ++neneg; else No_W[iv_turb2][n] = e;
      }
      cout << "*** negative k/e nodes: " << nkneg << '/' << neneg << endl;
    }

    if (temperature && !scalar_coupling)
      Update_temperature();

    if (periodic) {
      for (int ig=1; ig<=Nbcgroup; ig++)
        if (BCgroup[ig].type==IBPERE)
          for (int inb=1; inb<=BCgroup[ig].nnode; inb++)
            for (int iv=0; iv<Ncoupled; iv++)
              No_W[iv][Nobg[ig][inb].node] = No_W[iv][Nobg[ig][inb].twin];
    }

    if (turmod) {
      if (turbulence_coupling==0)
        Update_turbulence_sequential();  /* Update_turbulence_uncoupled(); */
      else if (turbulence_coupling==1) {
        Update_turbulence_coupled();
        Update_turbulence_coupled();
      }
    }

    cout << "update solution." << endl;


    massflux();

    /*
     * // zero nodal turbulent viscosity store
     * if (turmod)
     *   for (int inu=0; inu<Nnode; inu++)
     *     No_nuturb[inu] /= No_vol[inu];
     */

    /* write residuals and solution */
    writeres();
    writesoltp(file_output);
  }
  while
    ((logL2[iverr]>conv_thresh) && (iter<Niter));
  /* end of main loop */


  /* linear systems deallocation */
  if (temperature && !scalar_coupling)  delete ls_aztec_scalar;
  if (turbulence_coupling==1)           delete ls_aztec_turb;
  if (turmod && turbulence_coupling==0) delete ls_aztec_turb1;
  if (turmod && turbulence_coupling==0) delete ls_aztec_turb2;
  delete ls_aztec_coupled;


  return 0;
}

