
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_config.h"                    // required by ML
#include "ml_MultiLevelPreconditioner.h"  // ...
#include "ml_include.h"                   // ...

#if defined(HAVE_MPI)
  #include "mpi.h"
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "mfactory.h"
#include "ls_trilinos_aztecoo.h"


using namespace std;
using namespace m;


Register< mlinearsystem< double >,ls_trilinos_aztecoo > ls_trilinos_aztecoo("ls_trilinos_aztecoo","linear system solver, using Trilinos AztecOO (double p.)");


namespace m {


ls_trilinos_aztecoo::ls_trilinos_aztecoo() :
  mlinearsystem< double >(),
  issetup(false)
{
  // set solver options (surelly there is a way without arrays and Teuchos...)
#if 1
  // - options array first
  m_solver.SetAztecDefaults();
  az_options[AZ_solver]          = AZ_gmres;
  az_options[AZ_scaling]         = AZ_none;
  az_options[AZ_precond]         = AZ_dom_decomp;
  az_options[AZ_subdomain_solve] = AZ_ilut;
  az_options[AZ_conv]            = AZ_r0;
  az_options[AZ_output]          = AZ_all;
  az_options[AZ_pre_calc]        = AZ_calc;
  az_options[AZ_max_iter]        = 15;
  az_options[AZ_aux_vec]         = AZ_resid;
  az_options[AZ_keep_info]       = 0;
  m_solver.SetAllAztecOptions(az_options);
  // (alternative) see Aztec 2.1 User Guide for a complete list
  // m_solver.SetAztecDefaults();
  m_solver.SetAztecOption(AZ_solver,AZ_gmres);

  // - params array next
  az_params[AZ_tol]   = 1.e-60;
  az_params[AZ_drop]  = 0.;
  az_params[AZ_omega] = 1.;
  m_solver.SetAllAztecParams(az_params);
#endif

#if 1
  // - now Teuchos
  Teuchos::ParameterList MLList;
  ML_Epetra::SetDefaults("DD", MLList);
  // MLList.set("max levels", 5);
  // MLList.set("increasing or decreasing", "increasing");
  // MLList.set("aggregation: type", "METIS");
  // MLList.set("aggregation: nodes per aggregate", 16);
  MLList.set("smoother: type", "Gauss-Seidel");
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: type","Amesos-Superlu");
/*
  ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
  m_solver.SetPrecOperator(MLPrec);  
*/
#endif
}


ls_trilinos_aztecoo::~ls_trilinos_aztecoo()
{
  if (issetup) {
  }
}


void ls_trilinos_aztecoo::solve()
{
  vector< int > vNumEntriesPerRow(Ne,0);
  for (unsigned i=0; i<Ne; ++i)
    vNumEntriesPerRow[i] = m_A.IA[i+1] - m_A.IA[i];

  Epetra_Map EpetraMap((int) Ne,0,Epetra_SerialComm());

  Epetra_CrsMatrix A(Copy,EpetraMap,&vNumEntriesPerRow[0]);

  for (int i=0; i<m_A.NNU; ++i)
    A.InsertGlobalValues(i,m_A.IA[i+1]-m_A.IA[i],&m_A.A[ m_A.IA[i] ],&m_A.JA[ m_A.IA[i] ]);
  A.FillComplete();

  // set linear system vectors and AztecOO solver
  Epetra_Vector* X = new Epetra_Vector(View, A.OperatorDomainMap(), &m_X[0]);
  Epetra_Vector* B = new Epetra_Vector(View, A.OperatorRangeMap(), &m_B[0]);
  m_solver.SetUserMatrix(&A);
  m_solver.SetLHS(X);
  m_solver.SetRHS(B);

  m_solver.Iterate(az_options[AZ_max_iter], az_params[AZ_tol]);

  delete X;
  delete B;
}


}  // namespace m

