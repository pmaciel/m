
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
  // allocate options, parameters, status, proc_config and update
  proc_config = new int[AZ_PROC_SIZE];  for (int i=0; i<AZ_PROC_SIZE;    ++i) proc_config[i] = 0;
  options = new int[AZ_OPTIONS_SIZE];   for (int i=0; i<AZ_OPTIONS_SIZE; ++i) options[i] = 0;
  params = new double[AZ_PARAMS_SIZE];  for (int i=0; i<AZ_PARAMS_SIZE;  ++i) params[i] = 0.;
  status = new double[AZ_STATUS_SIZE];  for (int i=0; i<AZ_STATUS_SIZE;  ++i) status[i] = 0.;

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
  cout << "aztec: AZ_processor_info..." << endl;
  AZ_processor_info(proc_config);
  cout << "aztec: AZ_processor_info." << endl;

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


void ls_aztec::set_az_options(XMLNode& xml, int *options, double *params)
{
  // set matrix type
  const string mtype_str = xml.getAttribute< string >("mtype","msr");
  mtype = mtype_str=="vbr"? AZ_VBR_MATRIX : AZ_MSR_MATRIX;

  // set maximum iterations, number of Krylov subspaces and tolerance
  options[AZ_max_iter] = xml.getAttribute< int    >("max_iter",options[AZ_max_iter]);
  options[AZ_kspace]   = xml.getAttribute< int    >("kspace",  options[AZ_kspace]);
  params[AZ_tol]       = xml.getAttribute< double >("tol",     params[AZ_tol]);

  // set preconditioner, overlap, solver and output
  const string precond = xml.getAttribute< string >("precond"),
               overlap = xml.getAttribute< string >("overlap"),
               solver  = xml.getAttribute< string >("solver"),
               output  = xml.getAttribute< string >("output");
  const int output_int = xml.getAttribute< int    >("output",options[AZ_output]);
  options[AZ_precond] = (precond=="Jacobi"?   AZ_Jacobi  :  // Jacobi
                        (precond=="sym_GS"?   AZ_sym_GS  :  // symmetric Gauss-Siedel
                        (precond=="Neumann"?  AZ_Neumann :  // Neumann series polynomial
                        (precond=="ls"?       AZ_ls      :  // least-squares polynomial
                        (precond=="ilu"?      AZ_ilu     :  // domain decomp with ilu in subdomains
                        (precond=="bilu"?     AZ_bilu    :  // domain decomp with block ilu in subdomains
                        (precond=="bmilu"?    AZ_bmilu   :  // domain decomp with block milu in subdomains
                        (precond=="icc"?      AZ_icc     :  // domain decomp with incomp Choleski in domains
                        (precond=="none"?     AZ_none    :  // no preconditioning
                        options[AZ_precond] )))))))));  // (default)

  options[AZ_overlap] = (overlap=="diag"?     AZ_diag     :  // use diagonal blocks for overlapping
                        (overlap=="sym_full"? AZ_sym_full :  // use external rows (symmetric) for overlapping
                        (overlap=="full"?     AZ_full     :  // use external rows for overlapping
                        (overlap=="none"?     AZ_none     :  // no overlap
                                              options[AZ_overlap] ))));  // (default)

  options[AZ_solver] = (solver=="cg"?        AZ_cg       :  // preconditioned conjugate gradient method
                       (solver=="gmres"?     AZ_gmres    :  // preconditioned gmres method
                       (solver=="cgs"?       AZ_cgs      :  // preconditioned cg squared method
                       (solver=="tfqmr"?     AZ_tfqmr    :  // preconditioned transpose-free qmr method
                       (solver=="bicgstab"?  AZ_bicgstab :  // preconditioned stabilized bi-cg method
                       (solver=="slu"?       AZ_slu      :  // super LU direct method
                       (solver=="symmlq"?    AZ_symmlq   :  // indefinite symmetric like symmlq
                       (solver=="lu"?        AZ_lu       :  // sparse LU direct method
                                             options[AZ_solver] ))))))));  // (default)

  options[AZ_output] = (output=="all"?      AZ_all      :  // print matrix/vectors and residuals
                       (output=="none"?     AZ_none     :  // no intermediate results
                       (output=="warnings"? AZ_warnings :  // print only warnings
                       (output=="last"?     AZ_last     :  // print only last residual
                       (output_int>0?       output_int  :  // print residual every output_int iteration
                                            options[AZ_output] )))));  // (default)

  // override polynomial order if Jacobi preconditioner is used
  options[AZ_poly_ord] = options[AZ_precond]==AZ_Jacobi? 3:0;
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
  for (int i=0; i<1/*5 && status[AZ_its]<5.*/; ++i) {
    AZ_solve( &m_X[0], &m_B[0], options, params,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.indx,
      mtype==AZ_MSR_MATRIX? m_A_msr.bindx : m_A_vbr.bindx,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.rpntr,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.cpntr,
      mtype==AZ_MSR_MATRIX? NULL          : m_A_vbr.bpntr,
      mtype==AZ_MSR_MATRIX? m_A_msr.val   : m_A_vbr.val,
      data_org, status, proc_config );
    cout << "aztec status:"
         << "  [AZ_its]=" << status[AZ_its]
         << "  [AZ_r]="   << status[AZ_r]
         << "  [AZ_why]=" << (status[AZ_why]==AZ_normal?    "AZ_normal"    :
                             (status[AZ_why]==AZ_param?     "AZ_param"     :
                             (status[AZ_why]==AZ_breakdown? "AZ_breakdown" :
                             (status[AZ_why]==AZ_loss?      "AZ_loss"      :
                             (status[AZ_why]==AZ_ill_cond?  "AZ_ill_cond"  :
                             (status[AZ_why]==AZ_maxits?    "AZ_maxits"    :
                                                            "unknown" )))))) << endl;
  }
}


}  // namespace m

