
#include "mfactory.h"
#include "ls_aztec.h"

#include "az_aztec.h"

using namespace std;
using namespace m;


Register< mlinearsystem< double >,ls_aztec > ls_aztec("ls_aztec","linear system solver, using Aztec (double p.)");


namespace m {


ls_aztec::ls_aztec() : mlinearsystem< double >()
{
  // allocate options, parameters, status and proc_config
  proc_config = new int[AZ_PROC_SIZE];
  options = new int[AZ_OPTIONS_SIZE];
  params = new double[AZ_PARAMS_SIZE];
  status = new double[AZ_STATUS_SIZE];
}


ls_aztec::~ls_aztec()
{
  // deallocate proc_config, status, parameters and options
  delete[] status;
  delete[] params;
  delete[] options;
}


void ls_aztec::initialize(const vector< vector< unsigned > >& nz)
{
  // set sparse matrix
  m_A.initialize(nz);

  // set aztec
  setup();
  AZ_processor_info(proc_config);
  AZ_check_input(data_org,options,params,proc_config);
  AZ_check_msr(m_A.BINDX,m_A.NNU,0,AZ_LOCAL,proc_config);
}


void ls_aztec::setup()
{
  options[AZ_solver]   = AZ_cgs;
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_ls;
  options[AZ_output]   = 1;
  options[AZ_max_iter] = 640;
  options[AZ_poly_ord] = 7;
  params[AZ_tol]   = 0.0000001;
  params[AZ_drop]  = 0.;
//params[AZ_omega] = 1.;
}


void ls_aztec::solve()
{
  int error = 0;
  error = error;
  AZ_solve(&m_X[0], &m_B[0], options, params, NULL, m_A.BINDX,
    NULL, NULL, NULL, m_A.VAL, data_org, status, proc_config);
  cout << "Aztec status: " << endl
       << "  AZ_its: " << status[AZ_its] << endl
       << "  AZ_why: " << (status[AZ_why]==AZ_normal?    "AZ_normal"    :
                          (status[AZ_why]==AZ_param?     "AZ_param"     :
                          (status[AZ_why]==AZ_breakdown? "AZ_breakdown" :
                          (status[AZ_why]==AZ_loss?      "AZ_loss"      :
                          (status[AZ_why]==AZ_ill_cond?  "AZ_ill_cond"  :
                          (status[AZ_why]==AZ_maxits?    "AZ_maxits"    :
                                                         "unknown" )))))) << endl
       << "  AZ_r:   " << status[AZ_r]             << endl
       << "  AZ_scaled_r:      " << status[AZ_scaled_r]      << endl
       << "  AZ_rec_r:         " << status[AZ_rec_r]         << endl;
}


}  // namespace m

