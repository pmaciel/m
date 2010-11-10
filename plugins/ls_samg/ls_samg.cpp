
#include "mfactory.h"
#include "ls_samg.h"
#include "samg.h"


m::Register< m::mlinearsystem< double >,m::ls_samg > ls_samg("ls_samg","linear system solver, using SAMG (double p.)");


namespace m {


void ls_samg::initialize(const std::vector< std::vector< unsigned > >& nz)
{
  using std::cout;
  using std::endl;

  cout << "info: initialize matrix sparsity structure..." << endl;
  m_A.initialize(nz);

  // correct matrix diagonal entries to first position in the column indices
  int    *ia = m_A.ia;
  int    *ja = m_A.ja;
  double *a  = m_A.a;
  for (int r=0; r<m_A.nnu; ++r)
    if (r != ja[ia[r]-1]-1) {
      // diagonal is not first entry, so locate it...
      const int k1 = ia[r]-1;
      for (int k2=ia[r]; k2<ia[r+1]-1; ++k2)
        if (r == ja[k2]-1) {
          // ... and put it to the first position
          std::swap( a[k1], a[k2]);
          std::swap(ja[k1],ja[k2]);
          break;
        }
    }
  cout << "info: initialize matrix sparsity structure." << endl;


  cout << "info: initialize maps and parameters..." << endl;
  init_maps();
  init_params();
  cout << "info: initialize maps and parameters." << endl;


  char* nth = getenv("OMP_NUM_THREADS");
  char* lic = getenv("LM_LICENSE_FILE");
  cout << "info: number of threads: " << (nth? nth:"")            << " (OMP_NUM_THREADS: " << (nth? "":"not ") << "set)" << endl;
  cout << "info: license file: "      << (lic? lic:"license.dat") << " (LM_LICENSE_FILE: " << (nth? "":"not ") << "set)" << endl;
}


void ls_samg::solve()
{
  //FIXME: necessary before a second SAMG run if secondary parameters have to be reset (see manual)
  //SAMG_RESET_HIDDEN();

  p.ierr = 0;
  SAMG(
    &m_A.nnu, &m_A.nnz, &p.nsys, m_A.ia, m_A.ja, m_A.a, &m_B[0], &m_X[0],
    &(p.iu)[0], &p.ndiu, &(p.ip)[0], &p.ndip, &p.matrix, &(p.iscale)[0],
    &p.res_in, &p.res_out, &p.ncyc_done, &p.ierr, &p.nsolve, &p.ifirst, &p.eps,
    &p.ncyc, &p.iswtch, &p.a_cmplx, &p.g_cmplx, &p.p_cmplx, &p.w_avrge,
    &p.chktol, &p.idump, &p.iout );
  if (p.ierr) {
    const int e = std::max(p.ierr,-p.ierr);
    std::cerr << (p.ierr>0? "error: ":"warning: ") << p.ierr << ": "
      << (m_ierr.find(e)!=m_ierr.end()? m_ierr[e] : "(unknown)") << std::endl;
    if (p.ierr>0)
      throw 42;
  }
}


int ls_samg::digit(const int &number, const unsigned &place, unsigned radix) const
{
  // build vector of the number digits, reversed (!)
  std::vector< int > n;
  unsigned numberp = std::max(number,-number);
  while (numberp) {
    n.push_back(numberp % radix);
    numberp /= radix;
  }
  return (place>n.size()? -1 :
         (place<1?         0 : n[n.size()-place] ));
}


void ls_samg::init_block_iu()
{
  p.nsys = (int) Nb;
  p.ndiu = m_A.nnu;
  p.iu.assign(p.ndiu,0);
  for (int i=0, k=0; i<m_A.nnu/p.nsys; ++i)
    for (int j=0; j<p.nsys; ++j, ++k)
      p.iu[k] = j+1;
}


void ls_samg::init_block_ip()
{
  p.nsys = (int) Nb;
  p.ndip = m_A.nnu;
  p.ip.assign(p.ndip,0);
  for (int i=0, k=0; i<m_A.nnu/p.nsys; ++i)
    for (int j=0; j<p.nsys; ++j, ++k)
      p.ip[k] = i+1;
}


void ls_samg::init_maps()
{
  using std::make_pair;
  using aux::s_parm;

  // list primary main and sub-parameters pointers
  m_prim.push_back(s_parm( "matrix",        &p.matrix,   22   ));
  m_prim.push_back(s_parm(   "isym",        &p.matrix,  1,1   ));
  m_prim.push_back(s_parm(   "irow0",       &p.matrix,  2,2   ));
  m_prim.push_back(s_parm( "ifirst",        &p.ifirst,    1   ));
  m_prim.push_back(s_parm(   "itypu",       &p.ifirst,  1,1   ));
  m_prim.push_back(s_parm(   "ifirst_",     &p.ifirst,  2,0   ));
  m_prim.push_back(s_parm( "nsolve",        &p.nsolve,    2   ));
  m_prim.push_back(s_parm(   "napproach",   &p.nsolve,  1,1   ));
  m_prim.push_back(s_parm(   "nxtyp",       &p.nsolve,  2,2   ));
  m_prim.push_back(s_parm(   "internal",    &p.nsolve,  3,3   ));
  m_prim.push_back(s_parm(   "nprim",       &p.nsolve,  4,5   ));
  m_prim.push_back(s_parm(   "npr_is_dummy",&p.nsolve,  6,6   ));
  m_prim.push_back(s_parm(   "nint_weights",&p.nsolve,  7,7   ));
  m_prim.push_back(s_parm(   "nint_pat",    &p.nsolve,  8,8   ));
  m_prim.push_back(s_parm( "ncyc",          &p.ncyc,  11030   ));
  m_prim.push_back(s_parm(   "igam",        &p.ncyc,    1,1   ));
  m_prim.push_back(s_parm(   "ngrad",       &p.ncyc,    2,2   ));
  m_prim.push_back(s_parm(   "nkdim",       &p.ncyc,    3,3   ));
  m_prim.push_back(s_parm(   "ncycle",      &p.ncyc,    4,0   ));
  m_prim.push_back(s_parm( "iswtch",        &p.iswtch, 5140   ));
  m_prim.push_back(s_parm(   "iswit",       &p.iswtch,  1,1   ));
  m_prim.push_back(s_parm(   "iextent",     &p.iswtch,  2,2   ));
  m_prim.push_back(s_parm(   "n_default",   &p.iswtch,  3,4   ));
  m_prim.push_back(s_parm(   "norm_typ",    &p.iswtch,  5,5   ));
  m_prim.push_back(s_parm(   "ioscratch",   &p.iswtch,  6,7   ));
  m_prim.push_back(s_parm( "iout",          &p.iout,      2   ));
  m_prim.push_back(s_parm(   "iout1",       &p.iout,    1,1   ));
  m_prim.push_back(s_parm(   "iout2",       &p.iout,    2,2   ));
  m_prim.push_back(s_parm( "idump",         &p.idump,     0   ));
  m_prim.push_back(s_parm(   "idmp",        &p.idump,   1,1   ));
  m_prim.push_back(s_parm(   "igdp",        &p.idump,   2,2   ));
  m_prim.push_back(s_parm(   "iadp",        &p.idump,   3,3   ));
  m_prim.push_back(s_parm(   "iwdp",        &p.idump,   4,4   ));
  m_prim.push_back(s_parm(   "icdp",        &p.idump,   5,5   ));
  m_prim.push_back(s_parm( "eps",           &p.eps,     1.e-8 ));
  m_prim.push_back(s_parm( "a_cmplx",       &p.a_cmplx, 0.    ));
  m_prim.push_back(s_parm( "g_cmplx",       &p.g_cmplx, 0.    ));
  m_prim.push_back(s_parm( "p_cmplx",       &p.p_cmplx, 0.    ));
  m_prim.push_back(s_parm( "w_avrge",       &p.w_avrge, 0.    ));
  m_prim.push_back(s_parm( "chktol",        &p.chktol, -1.e0  ));

  // list secondary/hidden parameters function pointers
  m_sdry.push_back(s_parm( "samg_iset_cset_read",            &SAMG_ISET_CSET_READ,           ""       ));
  m_sdry.push_back(s_parm( "samg_iset_filnam",               &SAMG_ISET_FILNAM,              "temp"   ));
  m_sdry.push_back(s_parm( "samg_iset_filnam_dump",          &SAMG_ISET_FILNAM_DUMP,         "level"  ));
  m_sdry.push_back(s_parm( "samg_iset_ioform",               &SAMG_ISET_IOFORM,              "f"      ));
  m_sdry.push_back(s_parm( "samg_iset_logfile",              &SAMG_ISET_LOGFILE,             ""       ));
  m_sdry.push_back(s_parm( "samg_set_a_cmplx_agg_default",   &SAMG_SET_A_CMPLX_AGG_DEFAULT,   1.2     ));
  m_sdry.push_back(s_parm( "samg_set_a_cmplx_default",       &SAMG_SET_A_CMPLX_DEFAULT,       2.5     ));
  m_sdry.push_back(s_parm( "samg_set_allow_elim",            &SAMG_SET_ALLOW_ELIM,                 1  ));
  m_sdry.push_back(s_parm( "samg_set_b_cmplx",               &SAMG_SET_B_CMPLX,               0.      ));
  m_sdry.push_back(s_parm( "samg_set_b_cmplx_agg_default",   &SAMG_SET_B_CMPLX_AGG_DEFAULT,   1.2     ));
  m_sdry.push_back(s_parm( "samg_set_b_cmplx_default",       &SAMG_SET_B_CMPLX_DEFAULT,       2.      ));
  m_sdry.push_back(s_parm( "samg_set_blk_fillexp",           &SAMG_SET_BLK_FILLEXP,           0.5     ));
  m_sdry.push_back(s_parm( "samg_set_blk_stab",              &SAMG_SET_BLK_STAB,              0.25    ));
  m_sdry.push_back(s_parm( "samg_set_check_allpnts",         &SAMG_SET_CHECK_ALLPNTS,              1  ));
  m_sdry.push_back(s_parm( "samg_set_check_order",           &SAMG_SET_CHECK_ORDER,                0  ));
  m_sdry.push_back(s_parm( "samg_set_conv_stop",             &SAMG_SET_CONV_STOP,             0.1e-1  ));
  m_sdry.push_back(s_parm( "samg_set_cset_lastpnt",          &SAMG_SET_CSET_LASTPNT,               0  ));
  m_sdry.push_back(s_parm( "samg_set_cset_lessvars",         &SAMG_SET_CSET_LESSVARS,              0  ));
  m_sdry.push_back(s_parm( "samg_set_cset_longrow",          &SAMG_SET_CSET_LONGROW,               0  ));
  m_sdry.push_back(s_parm( "samg_set_cset_zerodiag",         &SAMG_SET_CSET_ZERODIAG,              0  ));
  m_sdry.push_back(s_parm( "samg_set_delta_milu",            &SAMG_SET_DELTA_MILU,            1.0     ));
  m_sdry.push_back(s_parm( "samg_set_densx",                 &SAMG_SET_DENSX,                 0.2     ));
  m_sdry.push_back(s_parm( "samg_set_droptol_cl",            &SAMG_SET_DROPTOL_CL,            0.5e-2  ));
  m_sdry.push_back(s_parm( "samg_set_droptol_smo",           &SAMG_SET_DROPTOL_SMO,           0.5e-2  ));
  m_sdry.push_back(s_parm( "samg_set_dump_correctw",         &SAMG_SET_DUMP_CORRECTW,              0  ));
  m_sdry.push_back(s_parm( "samg_set_ecg",                   &SAMG_SET_ECG,                  21.25    ));
  m_sdry.push_back(s_parm( "samg_set_ecg_default",           &SAMG_SET_ECG_DEFAULT,          21.25    ));
  m_sdry.push_back(s_parm( "samg_set_eps_dd",                &SAMG_SET_EPS_DD,                0.9     ));
  m_sdry.push_back(s_parm( "samg_set_eps_diag",              &SAMG_SET_EPS_DIAG,              1.e-05  ));
  m_sdry.push_back(s_parm( "samg_set_eps_lsq",               &SAMG_SET_EPS_LSQ,               0.5e-02 ));
  m_sdry.push_back(s_parm( "samg_set_etr",                   &SAMG_SET_ETR,                  12.2     ));
  m_sdry.push_back(s_parm( "samg_set_etr_default",           &SAMG_SET_ETR_DEFAULT,          12.2     ));
  m_sdry.push_back(s_parm( "samg_set_ewt",                   &SAMG_SET_EWT,                   0.2     ));
  m_sdry.push_back(s_parm( "samg_set_ewt_default",           &SAMG_SET_EWT_DEFAULT,           0.2     ));
  m_sdry.push_back(s_parm( "samg_set_factor_app_var",        &SAMG_SET_FACTOR_APP_VAR,        1.e-10  ));
  m_sdry.push_back(s_parm( "samg_set_factor_quasi_res",      &SAMG_SET_FACTOR_QUASI_RES,     10.      ));
  m_sdry.push_back(s_parm( "samg_set_full_pivoting",         &SAMG_SET_FULL_PIVOTING,              0  ));
  m_sdry.push_back(s_parm( "samg_set_g_cmplx_agg_default",   &SAMG_SET_G_CMPLX_AGG_DEFAULT,   1.2     ));
  m_sdry.push_back(s_parm( "samg_set_g_cmplx_default",       &SAMG_SET_G_CMPLX_DEFAULT,       2.8     ));
  m_sdry.push_back(s_parm( "samg_set_iauto_stop",            &SAMG_SET_IAUTO_STOP,           1103003  ));
  m_sdry.push_back(s_parm( "samg_set_ib_cmplx",              &SAMG_SET_IB_CMPLX,              0.0     ));
  m_sdry.push_back(s_parm( "samg_set_ib_cmplx_agg_default",  &SAMG_SET_IB_CMPLX_AGG_DEFAULT,  1.2     ));
  m_sdry.push_back(s_parm( "samg_set_ib_cmplx_default",      &SAMG_SET_IB_CMPLX_DEFAULT,      2.0     ));
  m_sdry.push_back(s_parm( "samg_set_ibgs_pivot",            &SAMG_SET_IBGS_PIVOT,                 0  ));
  m_sdry.push_back(s_parm( "samg_set_icrits",                &SAMG_SET_ICRITS,                   110  ));
  m_sdry.push_back(s_parm( "samg_set_ilu_speed",             &SAMG_SET_ILU_SPEED,                  1  ));
  m_sdry.push_back(s_parm( "samg_set_iodump",                &SAMG_SET_IODUMP,                    32  ));
  m_sdry.push_back(s_parm( "samg_set_iogrid",                &SAMG_SET_IOGRID,                    33  ));
  m_sdry.push_back(s_parm( "samg_set_iomovie",               &SAMG_SET_IOMOVIE,                   34  ));
  m_sdry.push_back(s_parm( "samg_set_ioscratch_default",     &SAMG_SET_IOSCRATCH_DEFAULT,         31  ));
  m_sdry.push_back(s_parm( "samg_set_ipass_max_set",         &SAMG_SET_IPASS_MAX_SET,              0  ));
  m_sdry.push_back(s_parm( "samg_set_irestriction_openmp",   &SAMG_SET_IRESTRICTION_OPENMP,        1  ));
  m_sdry.push_back(s_parm( "samg_set_iter_check",            &SAMG_SET_ITER_CHECK,                 3  ));
  m_sdry.push_back(s_parm( "samg_set_iter_pre",              &SAMG_SET_ITER_PRE,                   0  ));
  m_sdry.push_back(s_parm( "samg_set_itmax_conv",            &SAMG_SET_ITMAX_CONV,               200  ));
  m_sdry.push_back(s_parm( "samg_set_lastgrid",              &SAMG_SET_LASTGRID,                   0  ));
  m_sdry.push_back(s_parm( "samg_set_levelx",                &SAMG_SET_LEVELX,                    25  ));
  m_sdry.push_back(s_parm( "samg_set_lfil_cl",               &SAMG_SET_LFIL_CL,                    9  ));
  m_sdry.push_back(s_parm( "samg_set_lfil_smo",              &SAMG_SET_LFIL_SMO,                   9  ));
  m_sdry.push_back(s_parm( "samg_set_logio",                 &SAMG_SET_LOGIO,                      6  ));
  m_sdry.push_back(s_parm( "samg_set_max_level",             &SAMG_SET_MAX_LEVEL,                 25  ));
  m_sdry.push_back(s_parm( "samg_set_maxop_restart",         &SAMG_SET_MAXOP_RESTART,              2  ));
  m_sdry.push_back(s_parm( "samg_set_milu",                  &SAMG_SET_MILU,                       0  ));
  m_sdry.push_back(s_parm( "samg_set_mode_debug",            &SAMG_SET_MODE_DEBUG,                 1  ));
  m_sdry.push_back(s_parm( "samg_set_mode_mess",             &SAMG_SET_MODE_MESS,                  0  ));
  m_sdry.push_back(s_parm( "samg_set_multipass_allcoup",     &SAMG_SET_MULTIPASS_ALLCOUP,          0  ));
  m_sdry.push_back(s_parm( "samg_set_nblk_debug",            &SAMG_SET_NBLK_DEBUG,                 0  ));
  m_sdry.push_back(s_parm( "samg_set_nblk_max",              &SAMG_SET_NBLK_MAX,                1000  ));
  m_sdry.push_back(s_parm( "samg_set_nblk_overlap",          &SAMG_SET_NBLK_OVERLAP,               0  ));
  m_sdry.push_back(s_parm( "samg_set_nblk_resid",            &SAMG_SET_NBLK_RESID,                 0  ));
  m_sdry.push_back(s_parm( "samg_set_nblk_solve",            &SAMG_SET_NBLK_SOLVE,                 1  ));
  m_sdry.push_back(s_parm( "samg_set_nblk_solver",           &SAMG_SET_NBLK_SOLVER,                3  ));
  m_sdry.push_back(s_parm( "samg_set_ncframes",              &SAMG_SET_NCFRAMES,                  -1  ));
  m_sdry.push_back(s_parm( "samg_set_ncg",                   &SAMG_SET_NCG,                    50000  ));
  m_sdry.push_back(s_parm( "samg_set_ncgrad_default",        &SAMG_SET_NCGRAD_DEFAULT,             1  ));
  m_sdry.push_back(s_parm( "samg_set_ncyc_default",          &SAMG_SET_NCYC_DEFAULT,           11030  ));
  m_sdry.push_back(s_parm( "samg_set_ncyc_min",              &SAMG_SET_NCYC_MIN,                   0  ));
  m_sdry.push_back(s_parm( "samg_set_ncyc_start",            &SAMG_SET_NCYC_START,                 2  ));
  m_sdry.push_back(s_parm( "samg_set_neg_diag",              &SAMG_SET_NEG_DIAG,                  20  ));
  m_sdry.push_back(s_parm( "samg_set_neg_diag_brute",        &SAMG_SET_NEG_DIAG_BRUTE,            20  ));
  m_sdry.push_back(s_parm( "samg_set_nint_rowsum1",          &SAMG_SET_NINT_ROWSUM1,               0  ));
  m_sdry.push_back(s_parm( "samg_set_nkdim_default",         &SAMG_SET_NKDIM_DEFAULT,              6  ));
  m_sdry.push_back(s_parm( "samg_set_nmin_matrix",           &SAMG_SET_NMIN_MATRIX,            40000  ));
  m_sdry.push_back(s_parm( "samg_set_nmin_matrix_resc",      &SAMG_SET_NMIN_MATRIX_RESC,      100000  ));
  m_sdry.push_back(s_parm( "samg_set_nmin_vector",           &SAMG_SET_NMIN_VECTOR,           100000  ));
  m_sdry.push_back(s_parm( "samg_set_np_mod1",               &SAMG_SET_NP_MOD1,                    0  ));
  m_sdry.push_back(s_parm( "samg_set_np_mod2",               &SAMG_SET_NP_MOD2,                    0  ));
  m_sdry.push_back(s_parm( "samg_set_np_opt",                &SAMG_SET_NP_OPT,                     0  ));
  m_sdry.push_back(s_parm( "samg_set_nptmax",                &SAMG_SET_NPTMAX,                   300  ));
  m_sdry.push_back(s_parm( "samg_set_nptmn",                 &SAMG_SET_NPTMN,                    100  ));
  m_sdry.push_back(s_parm( "samg_set_nrc",                   &SAMG_SET_NRC,                        0  ));
  m_sdry.push_back(s_parm( "samg_set_nrc_default",           &SAMG_SET_NRC_DEFAULT,                7  ));
  m_sdry.push_back(s_parm( "samg_set_nrc_emergency",         &SAMG_SET_NRC_EMERGENCY,              3  ));
  m_sdry.push_back(s_parm( "samg_set_nrd",                   &SAMG_SET_NRD,                      131  ));
  m_sdry.push_back(s_parm( "samg_set_nru",                   &SAMG_SET_NRU,                      131  ));
  m_sdry.push_back(s_parm( "samg_set_nsolve_default",        &SAMG_SET_NSOLVE_DEFAULT,             2  ));
  m_sdry.push_back(s_parm( "samg_set_ntake_res_in",          &SAMG_SET_NTAKE_RES_IN,               0  ));
  m_sdry.push_back(s_parm( "samg_set_nth_res_scratch",       &SAMG_SET_NTH_RES_SCRATCH,            0  ));
  m_sdry.push_back(s_parm( "samg_set_ntr",                   &SAMG_SET_NTR,                        1  ));
  m_sdry.push_back(s_parm( "samg_set_numtry_max_set",        &SAMG_SET_NUMTRY_MAX_SET,             0  ));
  m_sdry.push_back(s_parm( "samg_set_nwt",                   &SAMG_SET_NWT,                        2  ));
  m_sdry.push_back(s_parm( "samg_set_omega_jacobi",          &SAMG_SET_OMEGA_JACOBI,          0.5     ));
  m_sdry.push_back(s_parm( "samg_set_p_cmplx_agg_default",   &SAMG_SET_P_CMPLX_AGG_DEFAULT,   1.2     ));
  m_sdry.push_back(s_parm( "samg_set_p_cmplx_default",       &SAMG_SET_P_CMPLX_DEFAULT,       2.8     ));
  m_sdry.push_back(s_parm( "samg_set_prim_norm",             &SAMG_SET_PRIM_NORM,                  0  ));
  m_sdry.push_back(s_parm( "samg_set_prim_print",            &SAMG_SET_PRIM_PRINT,                 0  ));
  m_sdry.push_back(s_parm( "samg_set_rcondx",                &SAMG_SET_RCONDX,                1.e-6   ));
  m_sdry.push_back(s_parm( "samg_set_show_un_res",           &SAMG_SET_SHOW_UN_RES,                0  ));
  m_sdry.push_back(s_parm( "samg_set_slow_coarsening",       &SAMG_SET_SLOW_COARSENING,       0.75    ));
  m_sdry.push_back(s_parm( "samg_set_stability",             &SAMG_SET_STABILITY,             0.25    ));
  m_sdry.push_back(s_parm( "samg_set_term_coarsening",       &SAMG_SET_TERM_COARSENING,       0.9     ));
  m_sdry.push_back(s_parm( "samg_set_w_avrge_agg_default",   &SAMG_SET_W_AVRGE_AGG_DEFAULT,   1.2     ));
  m_sdry.push_back(s_parm( "samg_set_w_avrge_default",       &SAMG_SET_W_AVRGE_DEFAULT,       2.5     ));

  // list secondary/hidden parameters function pointers (not in the manual)
  m_sdry.push_back(s_parm( "samg_iset_bnd_partition_file",    &SAMG_ISET_BND_PARTITION_FILE,   ""  ));
  m_sdry.push_back(s_parm( "samg_iset_filnam_uzawa",          &SAMG_ISET_FILNAM_UZAWA,         ""  ));
  m_sdry.push_back(s_parm( "samg_iset_iofile_opta",           &SAMG_ISET_IOFILE_OPTA,          ""  ));
  m_sdry.push_back(s_parm( "samg_iset_omp_partition_file",    &SAMG_ISET_OMP_PARTITION_FILE,   ""  ));
  m_sdry.push_back(s_parm( "samg_iset_partition_file",        &SAMG_ISET_PARTITION_FILE,       ""  ));
  m_sdry.push_back(s_parm( "samg_iset_scratchfile",           &SAMG_ISET_SCRATCHFILE,          ""  ));
  m_sdry.push_back(s_parm( "samg_set_allow_filnam_in",        &SAMG_SET_ALLOW_FILNAM_IN,       0   ));
  m_sdry.push_back(s_parm( "samg_set_alluns_at_allpnts",      &SAMG_SET_ALLUNS_AT_ALLPNTS,     0   ));
  m_sdry.push_back(s_parm( "samg_set_cl_cpl",                 &SAMG_SET_CL_CPL,                0   ));
  m_sdry.push_back(s_parm( "samg_set_cl_order",               &SAMG_SET_CL_ORDER,              0   ));
  m_sdry.push_back(s_parm( "samg_set_droptol",                &SAMG_SET_DROPTOL,               0.  ));
  m_sdry.push_back(s_parm( "samg_set_eps_abs",                &SAMG_SET_EPS_ABS,               0.  ));
  m_sdry.push_back(s_parm( "samg_set_eps_accept_jac_uzawa",   &SAMG_SET_EPS_ACCEPT_JAC_UZAWA,  0.  ));
  m_sdry.push_back(s_parm( "samg_set_eps_converg_jac_uzawa",  &SAMG_SET_EPS_CONVERG_JAC_UZAWA, 0.  ));
  m_sdry.push_back(s_parm( "samg_set_eps_diverg_jac_uzawa",   &SAMG_SET_EPS_DIVERG_JAC_UZAWA,  0.  ));
  m_sdry.push_back(s_parm( "samg_set_factor_res_var",         &SAMG_SET_FACTOR_RES_VAR,        0.  ));
  m_sdry.push_back(s_parm( "samg_set_factor_sol_jac_uzawa",   &SAMG_SET_FACTOR_SOL_JAC_UZAWA,  0.  ));
  m_sdry.push_back(s_parm( "samg_set_force_accel",            &SAMG_SET_FORCE_ACCEL,           0   ));
  m_sdry.push_back(s_parm( "samg_set_gmax_multipass",         &SAMG_SET_GMAX_MULTIPASS,        0.  ));
  m_sdry.push_back(s_parm( "samg_set_icase_jac_uzawa",        &SAMG_SET_ICASE_JAC_UZAWA,       0   ));
  m_sdry.push_back(s_parm( "samg_set_icolor_omp",             &SAMG_SET_ICOLOR_OMP,            0   ));
  m_sdry.push_back(s_parm( "samg_set_idtest_uzawa",           &SAMG_SET_IDTEST_UZAWA,          0   ));
  m_sdry.push_back(s_parm( "samg_set_ijac_uzawa",             &SAMG_SET_IJAC_UZAWA,            0   ));
  m_sdry.push_back(s_parm( "samg_set_ilu4allschwarz",         &SAMG_SET_ILU4ALLSCHWARZ,        0   ));
  m_sdry.push_back(s_parm( "samg_set_impldo_read_write_len",  &SAMG_SET_IMPLDO_READ_WRITE_LEN, 0   ));
  m_sdry.push_back(s_parm( "samg_set_inner_accel",            &SAMG_SET_INNER_ACCEL,           0   ));
  m_sdry.push_back(s_parm( "samg_set_iounit_opta",            &SAMG_SET_IOUNIT_OPTA,           0   ));
  m_sdry.push_back(s_parm( "samg_set_ipressure_uzawa",        &SAMG_SET_IPRESSURE_UZAWA,       0   ));
  m_sdry.push_back(s_parm( "samg_set_isat_uzawa",             &SAMG_SET_ISAT_UZAWA,            0   ));
  m_sdry.push_back(s_parm( "samg_set_iset_vio_dd",            &SAMG_SET_ISET_VIO_DD,           0   ));
  m_sdry.push_back(s_parm( "samg_set_isstep_uzawa",           &SAMG_SET_ISSTEP_UZAWA,          0   ));
  m_sdry.push_back(s_parm( "samg_set_isteering",              &SAMG_SET_ISTEERING,             0   ));
  m_sdry.push_back(s_parm( "samg_set_iswit3_user",            &SAMG_SET_ISWIT3_USER,           0   ));
  m_sdry.push_back(s_parm( "samg_set_itrace_jac_uzawa",       &SAMG_SET_ITRACE_JAC_UZAWA,      0   ));
  m_sdry.push_back(s_parm( "samg_set_itrace_schwarz",         &SAMG_SET_ITRACE_SCHWARZ,        0   ));
  m_sdry.push_back(s_parm( "samg_set_itrace_split",           &SAMG_SET_ITRACE_SPLIT,          0   ));
  m_sdry.push_back(s_parm( "samg_set_itrace_uzawa",           &SAMG_SET_ITRACE_UZAWA,          0   ));
  m_sdry.push_back(s_parm( "samg_set_levels_uzawa",           &SAMG_SET_LEVELS_UZAWA,          0   ));
  m_sdry.push_back(s_parm( "samg_set_max_cl_size",            &SAMG_SET_MAX_CL_SIZE,           0   ));
  m_sdry.push_back(s_parm( "samg_set_max_nb_list_size",       &SAMG_SET_MAX_NB_LIST_SIZE,      0   ));
  m_sdry.push_back(s_parm( "samg_set_maxiter_jac_uzawa",      &SAMG_SET_MAXITER_JAC_UZAWA,     0   ));
  m_sdry.push_back(s_parm( "samg_set_maxretry_jac_uzawa",     &SAMG_SET_MAXRETRY_JAC_UZAWA,    0   ));
  m_sdry.push_back(s_parm( "samg_set_miniter_jac_uzawa",      &SAMG_SET_MINITER_JAC_UZAWA,     0   ));
  m_sdry.push_back(s_parm( "samg_set_modify_mat",             &SAMG_SET_MODIFY_MAT,            0   ));
  m_sdry.push_back(s_parm( "samg_set_nbnd_omega",             &SAMG_SET_NBND_OMEGA,            0   ));
  m_sdry.push_back(s_parm( "samg_set_nbnd_omega_uzawa",       &SAMG_SET_NBND_OMEGA_UZAWA,      0   ));
  m_sdry.push_back(s_parm( "samg_set_nbnd_sweeps",            &SAMG_SET_NBND_SWEEPS,           0   ));
  m_sdry.push_back(s_parm( "samg_set_nbnd_sweeps_uzawa",      &SAMG_SET_NBND_SWEEPS_UZAWA,     0   ));
  m_sdry.push_back(s_parm( "samg_set_ndyn_smo",               &SAMG_SET_NDYN_SMO,              0   ));
  m_sdry.push_back(s_parm( "samg_set_notalluns_cheap",        &SAMG_SET_NOTALLUNS_CHEAP,       0   ));
  m_sdry.push_back(s_parm( "samg_set_nprim_at_allpnts",       &SAMG_SET_NPRIM_AT_ALLPNTS,      0   ));
  m_sdry.push_back(s_parm( "samg_set_nsimple_emergency",      &SAMG_SET_NSIMPLE_EMERGENCY,     0   ));
  m_sdry.push_back(s_parm( "samg_set_nsolve_uzawa",           &SAMG_SET_NSOLVE_UZAWA,          0   ));
  m_sdry.push_back(s_parm( "samg_set_nstar_typ",              &SAMG_SET_NSTAR_TYP,             0   ));
  m_sdry.push_back(s_parm( "samg_set_nsw_omega_uzawa",        &SAMG_SET_NSW_OMEGA_UZAWA,       0   ));
  m_sdry.push_back(s_parm( "samg_set_nsw_uzawa",              &SAMG_SET_NSW_UZAWA,             0   ));
  m_sdry.push_back(s_parm( "samg_set_ntr_prim",               &SAMG_SET_NTR_PRIM,              0   ));
  m_sdry.push_back(s_parm( "samg_set_ntyp_accel",             &SAMG_SET_NTYP_ACCEL,            0   ));
  m_sdry.push_back(s_parm( "samg_set_ntyp_galerkin",          &SAMG_SET_NTYP_GALERKIN,         0   ));
  m_sdry.push_back(s_parm( "samg_set_nxtyp_coarse",           &SAMG_SET_NXTYP_COARSE,          0   ));
  m_sdry.push_back(s_parm( "samg_set_omega_jac_es_uzawa",     &SAMG_SET_OMEGA_JAC_ES_UZAWA,    0.5 ));
  m_sdry.push_back(s_parm( "samg_set_omega_jac_p_uzawa",      &SAMG_SET_OMEGA_JAC_P_UZAWA,     0.5 ));
  m_sdry.push_back(s_parm( "samg_set_omega_uzawa",            &SAMG_SET_OMEGA_UZAWA,           0.5 ));
  m_sdry.push_back(s_parm( "samg_set_perf_meter_enabled",     &SAMG_SET_PERF_METER_ENABLED,    0   ));
  m_sdry.push_back(s_parm( "samg_set_print_nxtyp_coarse_msg", &SAMG_SET_PRINT_NXTYP_COARSE_MSG,0   ));
  m_sdry.push_back(s_parm( "samg_set_read_uzawa_parms",       &SAMG_SET_READ_UZAWA_PARMS,      0   ));
  m_sdry.push_back(s_parm( "samg_set_tau_uzawa",              &SAMG_SET_TAU_UZAWA,             0.  ));
  m_sdry.push_back(s_parm( "samg_set_use_ic",                 &SAMG_SET_USE_IC,                0   ));
  m_sdry.push_back(s_parm( "samg_set_vio_dd",                 &SAMG_SET_VIO_DD,                0   ));
  m_sdry.push_back(s_parm( "samg_set_write_cluster_id",       &SAMG_SET_WRITE_CLUSTER_ID,      0   ));

  // set error descriptions map
  m_ierr.insert(make_pair(   1,"general error (unclassified)"                                               ));
  m_ierr.insert(make_pair(   5,"error in license checking"                                                  ));
  m_ierr.insert(make_pair(  10,"insufficient dimensioning (happens only if memory extension is turned off)" ));
  m_ierr.insert(make_pair(  20,"illegal input parameter"                                                    ));
  m_ierr.insert(make_pair(  30,"undefined or missing input"                                                 ));
  m_ierr.insert(make_pair(  40,"error in input arrays"                                                      ));
  m_ierr.insert(make_pair(  50,"incorrect or inconsistent input"                                            ));
  m_ierr.insert(make_pair(  51,"SAMG re-start terminated: initial call with iswit=4 required"               ));
  m_ierr.insert(make_pair(  52,"SAMG re-start terminated: no AMG decomposition available"                   ));
  m_ierr.insert(make_pair(  53,"SAMG re-start terminated: illegal change in parameters"                     ));
  m_ierr.insert(make_pair(  54,"SAMG re-start terminated: illegal change in matrix"                         ));
  m_ierr.insert(make_pair(  60,"memory management failed (including file I/O to scratch files)"             ));
  m_ierr.insert(make_pair(  70,"allocation and de-allocation errors"                                        ));
  m_ierr.insert(make_pair(  71,"unexpected allocation status of allocatable array"                          ));
  m_ierr.insert(make_pair(  72,"unexpected assoziation status of pointer array"                             ));
  m_ierr.insert(make_pair(  80,"requested AMG component not installed in current release"                   ));
  m_ierr.insert(make_pair(  90,"logfile exists but could not be opened"                                     ));
  m_ierr.insert(make_pair(  91,"logfile already connected to different unit"                                ));
  m_ierr.insert(make_pair(  92,"logfile does not exist and could not be opened"                             ));
  m_ierr.insert(make_pair(  93,"specified unit does not exist"                                              ));
  m_ierr.insert(make_pair(  94,"unit=5: specify a logfile or another unit number"                           ));
  m_ierr.insert(make_pair( 100,"AMG setup: general error (unclassified)"                              ));
  m_ierr.insert(make_pair( 101,"AMG setup: setup failed in CNTRL-routine (automatic setup mechanism)" ));
  m_ierr.insert(make_pair( 110,"AMG setup: error in defining strong connectivity"                     ));
  m_ierr.insert(make_pair( 120,"AMG setup: error in the splitting process"                            ));
  m_ierr.insert(make_pair( 130,"AMG setup: error in the coarsening process"                           ));
  m_ierr.insert(make_pair( 140,"AMG setup: error in defining interpolation"                           ));
  m_ierr.insert(make_pair( 150,"AMG setup: error in computing the coarse-level Galerkin operators"    ));
  m_ierr.insert(make_pair( 160,"AMG setup: error in performing the setup optimization"                ));
  m_ierr.insert(make_pair( 200,"AMG solution phase: general error (unclassified)"             ));
  m_ierr.insert(make_pair( 210,"AMG solution phase: divergence of the method"                 ));
  m_ierr.insert(make_pair( 220,"AMG solution phase: error in relaxation (smoothing)"          ));
  m_ierr.insert(make_pair( 230,"AMG solution phase: error in ILU (smoothing)"                 ));
  m_ierr.insert(make_pair( 240,"AMG solution phase: error in ILUT (smoothing)"                ));
  m_ierr.insert(make_pair( 250,"AMG solution phase: error in inter-grid transfers"            ));
  m_ierr.insert(make_pair( 260,"AMG solution phase: error in alternating \"Schwarz process\"" ));
  m_ierr.insert(make_pair( 310,"300-399 reserved for MPI-parallel SAMG: illegal use of routine in parallel context"   ));
  m_ierr.insert(make_pair( 320,"300-399 reserved for MPI-parallel SAMG: illegal use of routine in sequential context" ));
  m_ierr.insert(make_pair( 800,"auxiliary components: general error (unclassified)"                               ));
  m_ierr.insert(make_pair( 810,"auxiliary components: error in ILU (one-level)"                                   ));
  m_ierr.insert(make_pair( 820,"auxiliary components: error in ILUT (one-level)"                                  ));
  m_ierr.insert(make_pair( 830,"auxiliary components: error in conjugate gradient (CG)"                           ));
  m_ierr.insert(make_pair( 831,"auxiliary components: quasi residual check has been reached (CG)"                 ));
  m_ierr.insert(make_pair( 832,"auxiliary components: residual stagnation check has been reached (CG)"            ));
  m_ierr.insert(make_pair( 833,"auxiliary components: limit of accelerator re-starts has been reached (CG)"       ));
  m_ierr.insert(make_pair( 840,"auxiliary components: error in BiCGstab"                                          ));
  m_ierr.insert(make_pair( 841,"auxiliary components: quasi residual check has been reached (BiCGstab)"           ));
  m_ierr.insert(make_pair( 842,"auxiliary components: residual stagnation check has been reached (BiCGstab)"      ));
  m_ierr.insert(make_pair( 843,"auxiliary components: limit of accelerator re-starts has been reached (BiCGstab)" ));
  m_ierr.insert(make_pair( 855,"auxiliary components: error in GMRES"                                             ));
  m_ierr.insert(make_pair( 851,"auxiliary components: quasi residual check has been reached (GMRES)"              ));
  m_ierr.insert(make_pair( 852,"auxiliary components: residual stagnation check has been reached (GMRES)"         ));
  m_ierr.insert(make_pair( 853,"auxiliary components: limit of accelerator re-starts has been reached (GMRES)"    ));
  m_ierr.insert(make_pair( 900,"solution on coarsest level: general error (unclassified)"                                   ));
  m_ierr.insert(make_pair( 910,"solution on coarsest level: error in method #1 (iterative application of current smoother)" ));
  m_ierr.insert(make_pair( 920,"solution on coarsest level: error in method #2 or #4 (preconditioned CG)"                   ));
  m_ierr.insert(make_pair( 930,"solution on coarsest level: error in method #3 or #5 (preconditioned BiCGstab)"             ));
  m_ierr.insert(make_pair( 960,"solution on coarsest level: error in method #6 (full Gauss elimination)"                    ));
  m_ierr.insert(make_pair( 970,"solution on coarsest level: error in method #7 (sparse Gauss elimination)"                  ));
  m_ierr.insert(make_pair( 980,"solution on coarsest level: error in method #8 (least squares solver)"                      ));
  m_ierr.insert(make_pair( 990,"solution on coarsest level: error in method #9 (user defined solver)"                       ));
}


void ls_samg::init_params()
{
  using std::vector;

  // set primary and secondary parameters
  for (std::vector< aux::s_parm >::iterator i=m_prim.begin(); i!=m_prim.end(); ++i)
    if (xml.isAttributeSet(i->n.c_str()))
      i->set(xml.getAttribute(i->n.c_str()));
    else
      i->def();
  for (std::vector< aux::s_parm >::iterator i=m_sdry.begin(); i!=m_sdry.end(); ++i)
    if (xml.isAttributeSet(i->n.c_str()))
      i->set(xml.getAttribute(i->n.c_str()));

  // show primary (main) parameters
  for (std::vector< aux::s_parm >::iterator i=m_prim.begin(); i!=m_prim.end(); ++i)
    if (i->ismain())
      std::cout << "info: primary parameter \"" << i->n << "\": " << i->get() << std::endl;

  // set solution approach (scalar is a good start)
  p.nsys = 1;
  p.ndiu = 1;
  p.ndip = 1;
  p.iu.assign(p.ndiu,0);
  p.ip.assign(p.ndip,0);
  if      (Nb>1 && digit(p.nsolve,1)==1)  init_block_iu();
  else if (Nb>1 && digit(p.nsolve,1)==2)  init_block_iu();
  else if (digit(p.nsolve,1)>2) {         init_block_iu();  init_block_ip(); }

  // set unknowns scaling
  p.iscale.assign(p.nsys,0);
  if (xml.isAttributeSet("iscale")) {
    const int number = xml.getAttribute< int >("iscale",0);
    for (unsigned i=1; i<=(unsigned)(p.nsys+1); ++i)
      p.iscale[i-1] = (digit(number,i)? 1:0);
  }
}


}  // namespace m

