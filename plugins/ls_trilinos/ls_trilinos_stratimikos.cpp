
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_MsrMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#if defined(HAVE_MPI)
  #include "mpi.h"
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "mfactory.h"
#include "ls_trilinos_stratimikos.h"


using namespace std;
using namespace m;


Register< mlinearsystem< double >,ls_trilinos_stratimikos > ls_trilinos_stratimikos("ls_trilinos_stratimikos","linear system solver, using Trilinos Stratimikos (double p.)");


namespace m {


ls_trilinos_stratimikos::ls_trilinos_stratimikos() : mlinearsystem< double >()
{
}


void ls_trilinos_stratimikos::solve()
{
#if 0
  /*********************
  * general setup phase
  *********************/
  using Teuchos::rcp;
  using Teuchos::RCP;

  fflush(stdout); cout << flush;
  string xmlfilename(opts->trilinos.XMLfile);
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder(xmlfilename); // the most important in general setup

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream(); // TODO: decouple from fancyostream to ostream or to C stdout when possible
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
  Teuchos::CommandLineProcessor  clp(false); // false: don't throw exceptions
  linearSolverBuilder.setupCLP(&clp); // not used, TODO: see if can be removed safely since not really used
  double tol=1.234e-5; clp.setOption( "tol",            &tol,            "Tolerance to check against the scaled residual norm." );
  int argc=2; char* argv[]={"morpheus","--tol=5.678e-9"}; Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) _MERROR_(-1,"Emulated command line parsing for stratimikos failed");

  /**************
  * epetra setup
  **************/
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_MsrMatrix> epetra_A=rcp(new Epetra_MsrMatrix(opts->trilinos.proc_config, opts->trilinos.Amat));
  RCP<Epetra_Vector>    epetra_x=rcp(new Epetra_Vector(View, epetra_A->OperatorDomainMap(), opts->trilinos.x));
  RCP<Epetra_Vector>    epetra_b=rcp(new Epetra_Vector(View, epetra_A->OperatorRangeMap(), opts->trilinos.b));

  /**************************
  * wrapping epetra to thyra
  **************************/
  RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( epetra_A );
  RCP<Thyra::VectorBase<double> >         x = Thyra::create_Vector( epetra_x, A->domain() );
  RCP<const Thyra::VectorBase<double> >   b = Thyra::create_Vector( epetra_b, A->range() );

  /*****************************
  * solving the system actually
  *****************************/

  // r = b - A*x, initial L2 norm
  double nrm_r=0.;
  opts->systemResidual=-1.;
  {
    Epetra_Vector epetra_r(*epetra_b);
    Epetra_Vector epetra_A_x(epetra_A->OperatorRangeMap());
    epetra_A->Apply(*epetra_x,epetra_A_x);
    epetra_r.Update(-1.0,epetra_A_x,1.0);
    epetra_r.Norm2(&nrm_r);
  }

  // Reading in the solver parameters from the parameters file and/or from
  // the command line.  This was setup by the command-line options
  // set by the setupCLP(...) function above.
  linearSolverBuilder.readParameters(0); /* out.get() if want confirmation about the xml file within trilinos */
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = linearSolverBuilder.createLinearSolveStrategy(""); // create linear solver strategy
  lowsFactory->setVerbLevel((Teuchos::EVerbosityLevel)opts->trilinos.verbosity); // set verbosity

  // print back default and current settings
  if (opts->trilinos.dumpDefault!=0) {
    fflush(stdout); cout << flush;
    _MMESSAGE_(0,1,"Dumping Trilinos/Stratimikos solver defaults to files: 'trilinos_default.txt' and 'trilinos_default.xml'...\n");
    fflush(stdout); cout << flush;
    std::ofstream ofs("./trilinos_default.txt");
    linearSolverBuilder.getValidParameters()->print(ofs,PLPrintOptions().indent(2).showTypes(true).showDoc(true)); // the last true-false is the deciding about whether printing documentation to option or not
    ofs.flush();ofs.close();
    ofs.open("./trilinos_default.xml");
    Teuchos::writeParameterListToXmlOStream(*linearSolverBuilder.getValidParameters(),ofs);
    ofs.flush();ofs.close();
  }
  if (opts->trilinos.dumpCurrXML!=0) {
    fflush(stdout); cout << flush;
    _MMESSAGE_(0,1,"Dumping Trilinos/Stratimikos current settings to file: 'trilinos_current.xml'...\n");
    fflush(stdout); cout << flush;
    linearSolverBuilder.writeParamsFile(*lowsFactory,"./trilinos_current.xml");
  }

  // solve the matrix
  RCP<Thyra::LinearOpWithSolveBase<double> > lows = Thyra::linearOpWithSolve(*lowsFactory, A); // create solver
  Thyra::solve(*lows, Thyra::NOTRANS, *b, &*x); // solve

  // r = b - A*x, final L2 norm
  {
    Epetra_Vector epetra_r(*epetra_b);
    Epetra_Vector epetra_A_x(epetra_A->OperatorRangeMap());
    epetra_A->Apply(*epetra_x,epetra_A_x);
    epetra_r.Update(-1.0,epetra_A_x,1.0);
    opts->systemResidual=1./nrm_r;
    nrm_r=0.;
    epetra_r.Norm2(&nrm_r);
    opts->systemResidual*=nrm_r;
  }
  
  // print relative residual
  fflush(stdout); cout << flush;
  _MMESSAGE_(0,0,"L2 norm of linear system: % 8.8e\n",opts->systemResidual);
#endif
}


}  // namespace m

