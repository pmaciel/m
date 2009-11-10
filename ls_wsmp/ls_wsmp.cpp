
#include "mfactory.h"
#include "ls_wsmp.h"

using namespace std;
using namespace m;


Register< mlinearsystem< double >,ls_wsmp > ls_wsmp("ls_wsmp","linear system solver, using WSMP (double p.)");


// WSMP prototypes
extern "C" {
  void wsetmaxthrds_(int *);
  void wgsmp_(int*, int*, int*, double*, double*, int*, int*, double*, int*, double*);
}


namespace m {


ls_wsmp::ls_wsmp() : mlinearsystem< double >()
{
  int nthd = 1;
  char* nthreads = getenv("WSMP_NUM_THREADS");
  sscanf(nthreads? nthreads:"1","%d",&nthd);
  wsetmaxthrds_(&nthd);
  cout << "info: number of threads: " << nthd << " (WSMP_NUM_THREADS: " << (nthreads? "":"not ") << "set)" << endl;

  char* licpath  = getenv("WSMPLICPATH");
  cout << "info: license file path: " << (licpath? licpath:"(WSMPLICPATH: not set)") << endl;

  iparm[ 0] = 0;  // iparm/dparm defaults
  call_wsmp(0);   // ...
  iparm[ 4] = 0;  // + C-style numbering
  iparm[19] = 2;  // + ordering option 5
}


void ls_wsmp::solve()
{
  call_wsmp(1);  // analysis and reordering
  call_wsmp(2);  // LU factorization
  call_wsmp(3);  // forward and backward elimination
  call_wsmp(4);  // iterative refinement
  m_X.swap(m_B); // swap solution from B vector to X
  cout << "info: summary [it/res]: " << iparm[25] << '/' << dparm[25] << endl;
}


int ls_wsmp::call_wsmp(int task)
{
  cout << "info: wsmp task " << task << "..." << endl;
  iparm[ 1] = task;
  iparm[ 2] = task;
  int nrhs = 1;
  boost::timer t;
  wgsmp_(&m_A.NNU,m_A.IA,m_A.JA,m_A.A,&m_B[0],&m_A.NNU,&nrhs,NULL,iparm,dparm);
  if (iparm[63]) {
    cerr << "error: SolverWSMP task/error: " << task << '/' << iparm[63] << endl;
    throw 42;
  }
  switch (task) {
    case 1: {
        cout << "info: number of nonzeros in LU factors: " << iparm[23] << endl;
        cout << "info: number of FLOPS in factorization: " << dparm[23] << endl;
      } break;
    case 2: {
        cout << "info: factorization MFLOPS: " << dparm[22]*1.e-6/t.elapsed() << endl;
      } break;
    case 4: {
        cout << "info: maximum relative error: " << dparm[6] << endl;
      } break;
    default: break;
  }
  cout << "info: wsmp task/time: " << task << '/' << t.elapsed() << 's' << endl;
  return iparm[63];
}


}  // namespace m

