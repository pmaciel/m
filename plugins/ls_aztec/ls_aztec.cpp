
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


void ls_aztec::solve()
{
  int error = 0;
  error = error;
  AZ_solve(&m_X[0], &m_B[0], options, params, NULL, m_A.BINDX,
    NULL, NULL, NULL, m_A.VAL, data_org, status, proc_config);
  std::cout << "Aztec summary: " << std::endl
            << "  number of iterations: " << status[AZ_its] << std::endl
            << "  termination reason:   " << (status[AZ_why]==AZ_normal?    "user requested convergence criteria is satisfied" :
                                             (status[AZ_why]==AZ_param?     "user requested option is not available" :
                                             (status[AZ_why]==AZ_breakdown? "numerical breakdown occured" :
                                             (status[AZ_why]==AZ_loss?      "loss of precision occured" :
                                             (status[AZ_why]==AZ_ill_cond?  "GMRES Hessenberg matrix is ill-conditioned" :
                                             (status[AZ_why]==AZ_maxits?    "maximum iterations taken without convergence" :
                                                                            "unknown" )))))) << std::endl
            << "  true residual norm:   " << status[AZ_r] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl
            << "  " << status[AZ_its] << std::endl;
}


}  // namespace m

