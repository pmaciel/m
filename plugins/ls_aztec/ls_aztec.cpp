
#include "mfactory.h"
#include "ls_aztec.h"

#include "az_aztec.h"

using namespace std;
using namespace m;


Register< mlinearsystem< double >,ls_aztec > ls_aztec("ls_aztec","linear system solver, using Aztec (double p.)");


namespace m {


ls_aztec::ls_aztec() :
  mlinearsystem< double >(),
  mtype(AZ_MSR_MATRIX)
{
  // allocate options, parameters, status and proc_config
  options = new int[AZ_OPTIONS_SIZE];   for (int i=0; i<AZ_OPTIONS_SIZE; ++i) options[i] = 0;
  params = new double[AZ_PARAMS_SIZE];  for (int i=0; i<AZ_PARAMS_SIZE;  ++i) params[i] = 0.;
  status = new double[AZ_STATUS_SIZE];  for (int i=0; i<AZ_STATUS_SIZE;  ++i) status[i] = 0.;
  proc_config = new int[AZ_PROC_SIZE];  for (int i=0; i<AZ_PROC_SIZE;    ++i) proc_config[i] = 0;

  // setup options/params (at initialize, user option override)
  AZ_defaults(options,params);
  xml.addAttribute("mtype",    "msr");
  xml.addAttribute("output",   "1");
  xml.addAttribute("precond",  "none");
  xml.addAttribute("overlap",  "0");
  xml.addAttribute("solver",   "gmres");
  xml.addAttribute("max_iter", "500");
  xml.addAttribute("kspace",   "30");
  xml.addAttribute("tol",      "1.e-6");
  set_az_options(xml,options,params);
}


ls_aztec::~ls_aztec()
{
  // deallocate proc_config, status, parameters and options
  delete[] status;
  delete[] params;
  delete[] options;
  delete[] proc_config;

  free(update);    free(update_index);
  free(external);  free(extern_index);
}


void ls_aztec::initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb)
{
  cout << "options/params setup..." << endl;
  set_az_options(xml,options,params);
  cout << "options/params setup." << endl;

  cout << "parent/matrix initialize..." << endl;
  mlinearsystem< double >::initialize(_Ne,_Nv,_Nb);
  if (mtype==AZ_MSR_MATRIX)
    m_A_msr.initialize(Ne,Nv,Nb);
  else
    m_A_vbr.initialize(Ne,Nv,Nb);
  cout << "parent/matrix initialize." << endl;
}


void ls_aztec::initialize(const vector< vector< unsigned > >& nz)
{
  cout << "aztec: AZ_set_proc_config..." << endl;
  AZ_set_proc_config(proc_config,AZ_NOT_MPI);
  cout << "aztec: AZ_set_proc_config." << endl;

  cout << "aztec: initialize, AZ_{read_update,transform,check_" << (mtype==AZ_MSR_MATRIX? "msr":"vbr") << "}..." << endl;
  if (mtype==AZ_MSR_MATRIX) {
    m_A_msr.initialize(nz);
    int N_update = 0;
    AZ_read_update(&N_update,&update,proc_config,Ne*Nb,1,AZ_linear);
    AZ_transform( proc_config, &external, m_A_msr.bindx, m_A_msr.val,
      update, &update_index, &extern_index, &data_org, N_update,
      NULL, NULL, NULL, NULL, AZ_MSR_MATRIX );
    AZ_check_msr(m_A_msr.bindx,N_update,0,AZ_GLOBAL,proc_config);
  }
  else {
    m_A_vbr.initialize(nz);
    int N_update = 0;
    AZ_read_update(&N_update,&update,proc_config,Ne,1,AZ_linear);
    AZ_transform( proc_config, &external, m_A_vbr.bindx, m_A_vbr.val,
      update, &update_index, &extern_index, &data_org, N_update,
      m_A_vbr.indx, m_A_vbr.bpntr, m_A_vbr.rpntr, &m_A_vbr.cpntr, AZ_VBR_MATRIX );
    AZ_check_vbr( N_update, 0, AZ_GLOBAL,
      m_A_vbr.bindx, m_A_vbr.bpntr, m_A_vbr.cpntr, m_A_vbr.rpntr, proc_config);
  }
  cout << "aztec: initialize, AZ_{read_update,transform,check_" << (mtype==AZ_MSR_MATRIX? "msr":"vbr") << "}." << endl;

  cout << "aztec: AZ_check_input..." << endl;
  AZ_check_input(data_org,options,params,proc_config);
  cout << "aztec: AZ_check_input." << endl;
}


const double& ls_aztec::A(const unsigned r, const unsigned c) const
{
  return mtype==AZ_MSR_MATRIX? m_A_msr(r,c):m_A_vbr(r,c);
}


double& ls_aztec::A(const unsigned r, const unsigned c)
{
  return mtype==AZ_MSR_MATRIX? m_A_msr(r,c):m_A_vbr(r,c);
}


const double& ls_aztec::A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const
{
  return mtype==AZ_MSR_MATRIX? m_A_msr(R,C,r,c):m_A_vbr(R,C,r,c);
}


double& ls_aztec::A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)
{
  return mtype==AZ_MSR_MATRIX? m_A_msr(R,C,r,c):m_A_vbr(R,C,r,c);
}


void ls_aztec::set_az_options(XMLNode& _xml, int *_options, double *_params)
{
  // set matrix type
  const string mtype_str = _xml.getAttribute< string >("mtype","msr");
  mtype = mtype_str=="vbr"? AZ_VBR_MATRIX : AZ_MSR_MATRIX;

  // set options
  const string solver   = _xml.getAttribute< string >("solver"),
               scaling  = _xml.getAttribute< string >("scaling"),
               precond  = _xml.getAttribute< string >("precond"),
               subdomain_solve = _xml.getAttribute< string >("subdomain_solve"),
               conv     = _xml.getAttribute< string >("conv"),
               output   = _xml.getAttribute< string >("output"),
               pre_calc = _xml.getAttribute< string >("pre_calc"),
               overlap  = _xml.getAttribute< string >("overlap"),
               type_overlap    = _xml.getAttribute< string >("type_overlap"),
               orthog   = _xml.getAttribute< string >("orthog"),
               aux_vec  = _xml.getAttribute< string >("aux_vec");
  const int output_int  = _xml.getAttribute< int >("output",  _options[AZ_output]),
            overlap_int = _xml.getAttribute< int >("overlap", _options[AZ_overlap]);
  _options[AZ_solver] = (solver=="cg"?        AZ_cg       :  // preconditioned conjugate gradient method
                        (solver=="gmres"?     AZ_gmres    :  // preconditioned gmres method
                        (solver=="cgs"?       AZ_cgs      :  // preconditioned cg squared method
                        (solver=="tfqmr"?     AZ_tfqmr    :  // preconditioned transpose-free qmr method
                        (solver=="bicgstab"?  AZ_bicgstab :  // preconditioned stabilized bi-cg method
                        (solver=="slu"?       AZ_slu      :  // super LU direct method
                        (solver=="symmlq"?    AZ_symmlq   :  // indefinite symmetric like symmlq
                        (solver=="GMRESR"?    AZ_GMRESR   :  // recursive GMRES (not supported)
                        (solver=="fixed_pt"?  AZ_fixed_pt :  // fixed point iteration
                        (solver=="analyze"?   AZ_analyze  :  // fixed point iteration
                        (solver=="lu"?        AZ_lu       :  // sparse LU direct method (also used for preconditioning)
                                              _options[AZ_solver] )))))))))));  // (default)

  _options[AZ_scaling] = (scaling=="none"?        AZ_none        :  // no scaling
                         (scaling=="Jacobi"?      AZ_Jacobi      :  // Jacobi scaling
                         (scaling=="BJacobi"?     AZ_BJacobi     :  // block Jacobi scaling
                         (scaling=="row_sum"?     AZ_row_sum     :  // point row-sum scaling
                         (scaling=="sym_diag"?    AZ_sym_diag    :  // symmetric diagonal scaling
                         (scaling=="sym_row_sum"? AZ_sym_row_sum :  // symmetric diagonal scaling
                         (scaling=="equil"?       AZ_equil       :  // equilib scaling
                         (scaling=="sym_BJacobi"? AZ_sym_BJacobi :  // symmetric block Jacobi scaling
                                                  _options[AZ_scaling] ))))))));  // (default)

  _options[AZ_precond] = (precond=="none"?         AZ_none         :  // no preconditioning
                         (precond=="Jacobi"?       AZ_Jacobi       :  // Jacobi (also used for scaling options)
                         (precond=="sym_GS"?       AZ_sym_GS       :  // symmetric Gauss-Siedel
                         (precond=="Neumann"?      AZ_Neumann      :  // Neumann series polynomial
                         (precond=="ls"?           AZ_ls           :  // least-squares polynomial
                         (precond=="ilu"?          AZ_ilu          :  // domain decomposition with ilu in subdomains
                         (precond=="bilu"?         AZ_bilu         :  // domain decomposition with block ilu in subdomains
                         (precond=="lu"?           AZ_lu           :  // domain decomposition with lu in subdomains
                         (precond=="icc"?          AZ_icc          :  // domain decomposition with incomplete Choleski in domains
                         (precond=="ilut"?         AZ_ilut         :  // domain decomposition with ilut in subdomains
                         (precond=="rilu"?         AZ_rilu         :  // domain decomposition with rilu in subdomains
                         (precond=="recursive"?    AZ_recursive    :  // recursive call to AZ_iterate()
                         (precond=="smoother"?     AZ_smoother     :  // recursive call to AZ_iterate()
                         (precond=="dom_decomp"?   AZ_dom_decomp   :  // domain decomposition using subdomain solver given by options[AZ_subdomain_solve]
                         (precond=="multilevel"?   AZ_multilevel   :  // do multiplicative domain decomposition with coarse grid (not supported)
                         (precond=="user_precond"? AZ_user_precond :  // user's preconditioning
                         (precond=="bilu_ifp"?     AZ_bilu_ifp     :  // domain decomposition with bilu using ifpack in subdom
                                                   _options[AZ_precond] )))))))))))))))));  // (default)

  _options[AZ_subdomain_solve] = (subdomain_solve=="lu"?   AZ_lu   :  // domain decomposition with lu in subdomains
                                 (subdomain_solve=="ilut"? AZ_ilut :  // domain decomposition with ilut in subdomains
                                 (subdomain_solve=="ilu"?  AZ_ilu  :  // domain decomposition with ilu in subdomains
                                 (subdomain_solve=="rilu"? AZ_rilu :  // domain decomposition with rilu in subdomains
                                 (subdomain_solve=="bilu"? AZ_bilu :  // domain decomposition with block ilu in subdomains
                                 (subdomain_solve=="icc"?  AZ_icc  :  // domain decomposition with incomplete Choleski in domains
                                                           _options[AZ_subdomain_solve] ))))));  // (default)

  _options[AZ_conv] = (conv=="r0"?              AZ_r0              :  // ||r||_2 / ||r^{(0)}||_2
                      (conv=="rhs"?             AZ_rhs             :  // ||r||_2 / ||b||_2
                      (conv=="Anorm"?           AZ_Anorm           :  // ||r||_2 / ||A||_infty
                      (conv=="sol"?             AZ_sol             :  // ||r||_infty/(||A||_infty ||x||_1+||b||_infty)
  //                  (conv=="weighted"?        AZ_weighted        :  // ||r||_WRMS
                      (conv=="expected_values"? AZ_expected_values :  // ||r||_WRMS with weights taken as |A||x0|
                      (conv=="noscaled"?        AZ_noscaled        :  // ||r||_2
                                                _options[AZ_conv] ))))));  // (default)

  _options[AZ_output] = (output=="all"?      AZ_all      :  // print out everything including matrix
                        (output=="none"?     AZ_none     :  // print out no results (not even warnings)
                        (output=="last"?     AZ_last     :  // print out final residual and warnings
                        (output=="warnings"? AZ_warnings :  // print out only warning messages
                        (output_int>0?       output_int  :  // print residual every output_int iteration
                                             _options[AZ_output] )))));  // (default)

  _options[AZ_pre_calc] = (pre_calc=="calc"?      AZ_calc       :  // use no previous information
                          (pre_calc=="recalc"?    AZ_recalc     :  // use last symbolic information
                          (pre_calc=="reuse"?     AZ_reuse      :  // use a previous factorization to precondition
                          (pre_calc=="sys_reuse"? AZ_sys_reuse  :  // use last factorization to precondition
                                                  _options[AZ_pre_calc] ))));  // (default)

  _options[AZ_overlap] = (overlap=="none"? AZ_none     :  // no overlap
                         (overlap=="diag"? AZ_diag     :  // use diagonal blocks for overlapping
                         (overlap=="full"? AZ_full     :  // use external rows for overlapping
                         (overlap_int>=0?  overlap_int :  // overlapping steps
                                           _options[AZ_overlap] ))));  // (default)

  _options[AZ_type_overlap] = (type_overlap=="standard"?  AZ_standard  :
                              (type_overlap=="symmetric"? AZ_symmetric :
                                                          _options[AZ_type_overlap] ));  // (default)

  _options[AZ_orthog] = (orthog=="classic"?  AZ_classic  :  // 2 steps of classical Gram-Schmidt orthogonalization
                        (orthog=="modified"? AZ_modified :  // modified Gram-Schmidt orthogonalization
                                             _options[AZ_orthog] ));  // (default)

  _options[AZ_aux_vec] = (aux_vec=="resid"? AZ_resid :  // r is set to the initial residal vector
                         (aux_vec=="rand"?  AZ_rand  :  // r is set to the random numbers between -1. and 1.
                                            _options[AZ_aux_vec] ));  // (default)

  _options[AZ_graph_fill] = _xml.getAttribute< int >("graph_fill", _options[AZ_graph_fill]);
  _options[AZ_max_iter]   = _xml.getAttribute< int >("max_iter",   _options[AZ_max_iter]);
  _options[AZ_poly_ord]   = _xml.getAttribute< int >("poly_ord",   _options[AZ_poly_ord]);
  _options[AZ_kspace]     = _xml.getAttribute< int >("kspace",     _options[AZ_kspace]);
  _options[AZ_reorder]    = _xml.getAttribute< int >("reorder",    _options[AZ_reorder]);
  _options[AZ_keep_info]  = _xml.getAttribute< int >("keep_info",  _options[AZ_keep_info]);

  // set params
  _params[AZ_tol]       = _xml.getAttribute< double >("tol",       _params[AZ_tol]);
  _params[AZ_drop]      = _xml.getAttribute< double >("drop",      _params[AZ_drop]);
  _params[AZ_ilut_fill] = _xml.getAttribute< double >("ilut_fill", _params[AZ_ilut_fill]);
  _params[AZ_omega]     = _xml.getAttribute< double >("omega",     _params[AZ_omega]);
  // _params[AZ_weights] = // harder to implement, not done
}


void ls_aztec::reset(const double &v)
{
  mlinearsystem< double >::reset(v);
  if (mtype==AZ_MSR_MATRIX)
    m_A_msr.reset(v);
  else
    m_A_vbr.reset(v);
}


void ls_aztec::zerorow(const unsigned r)
{
  zerorow(r/Nb,r%Nb);
}


void ls_aztec::zerorow(const unsigned R, const unsigned r)
{
  B(R,r) = 0.;
  if (mtype==AZ_MSR_MATRIX)
    m_A_msr.zerorow(R,r);
  else
    m_A_vbr.zerorow(R,r);
}


void ls_aztec::print(std::ostream& o, bool pmatrix)
{
/*  mlinearsystem< double >::print(o,false);return;*/
  mlinearsystem< double >::print(o,false);
  if (!pmatrix)
    return;
  if (mtype==AZ_MSR_MATRIX)
    AZ_print_out(
      update_index, extern_index, update, external,
      m_A_msr.val, NULL, m_A_msr.bindx, NULL, NULL, NULL,
      proc_config, AZ_explicit, mtype, Ne, 0, 0 );
  else
    AZ_print_out(
      update_index, extern_index, update, external,
      m_A_vbr.val, m_A_vbr.indx, m_A_vbr.bindx, m_A_vbr.rpntr, m_A_vbr.cpntr, m_A_vbr.bpntr,
      proc_config, AZ_explicit, mtype, Ne, 0, 0 );
}


void ls_aztec::solve()
{
  status[AZ_its] = 0.;
  /*for (int i=0; i<5 && status[AZ_its]<5.; ++i)*/ {
    AZ_solve( &m_X[0], &m_B[0], options, params,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.indx,
      mtype==AZ_MSR_MATRIX? m_A_msr.bindx : m_A_vbr.bindx,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.rpntr,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.cpntr,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.bpntr,
      mtype==AZ_MSR_MATRIX? m_A_msr.val   : m_A_vbr.val,
      data_org, status, proc_config );
    cout << "aztec status:"
         << "  [its]=" << status[AZ_its]
         << "  [why]=" << (status[AZ_why]==AZ_normal?    "normal"    :
                          (status[AZ_why]==AZ_param?     "param"     :
                          (status[AZ_why]==AZ_breakdown? "breakdown" :
                          (status[AZ_why]==AZ_loss?      "loss"      :
                          (status[AZ_why]==AZ_ill_cond?  "ill_cond"  :
                          (status[AZ_why]==AZ_maxits?    "maxits"    :
                                                         "(unknown)" ))))))
         << "  [r]="             << status[AZ_r]
         << "  [scaled_r]="      << status[AZ_scaled_r]
         << "  [rec_r]="         << status[AZ_rec_r]
         << "  [solve_time]="    << status[AZ_solve_time]
    //   << "  [Aztec_version]=" << status[AZ_Aztec_version]
         << endl;
  }
}


}  // namespace m

