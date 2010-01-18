
#include <cstdio>  // for sscanf
#include "mfactory.h"
#include "ls_pardiso.h"

using namespace std;
using namespace m;


Register< mlinearsystem< double >,ls_pardiso > ls_pardiso("ls_pardiso","linear system solver, using Pardiso (double p.)");


// Pardiso prototypes
extern "C" {
  int pardisoinit_(
    void *, int *, int *, int *, double *, int *);
  int pardiso_(
    void *, int *, int *, int *, int *, int *,
    double *, int *, int *, int *, int *, int *,
    int *, double *, double *, int *, double *);
}


namespace m {


ls_pardiso::ls_pardiso() : mlinearsystem< double >()
{
  for (int i=0; i<64; ++i)  iparm[i] = 0;
  for (int i=0; i<64; ++i)  dparm[i] = 0.;

  char* nthreads = getenv("OMP_NUM_THREADS");
  sscanf(nthreads? nthreads:"1","%d",&iparm[2]);
  cout << "info: number of threads: " << iparm[2] << " (OMP_NUM_THREADS: " << (nthreads? "":"not ") << "set)" << endl;

  iparm[ 7] = 0;  // max numbers of iterative refinement steps
  iparm[31] = 0;  // [0|1] sparse direct solver or multi-recursive iterative solver
  maxfct    = 1,  // maximum number of numerical factorizations
  mnum      = 1,  // which factorization to use
  mtype     = 1,  // real structurally symmetric matrix
  nrhs      = 1;  // number of right hand sides
}


void ls_pardiso::solve()
{
  call_pardisoinit();  // setup
  call_pardiso(11,1);  // 11: reordering and symbolic factorization
  call_pardiso(22,0);  // 22: numerical factorization and
  call_pardiso(33,0);  // 33: back substitution and iterative refinement
  call_pardiso(-1,0);  // -1: termination and release of memory
}


int ls_pardiso::call_pardisoinit() {
  int error = 0;
  pardisoinit_(pt,&mtype,&iparm[31],iparm,dparm,&error);
  if (error) {
    cerr << "error: pardiso init error: "
         << (error==-10? "no license file pardiso.lic found"      :
            (error==-11? "license is expired"         :
            (error==-12? "wrong username or hostname" : "unknown error" ))) << endl;
    throw 42;
  }
  return error;
}


int ls_pardiso::call_pardiso(int phase, int msglvl)
{
  int error = 0;
  pardiso_( pt, &maxfct, &mnum, &mtype, &phase, &m_A.NNU,
            m_A.A, m_A.IA, m_A.JA, NULL, &nrhs, iparm,
            &msglvl, &m_B[0], &m_X[0], &error, dparm );
  if (error) {
    cerr << "error: pardiso phase/error: " << phase << '/' << error << ": "
         << (error==  -1? "input inconsistent" :
            (error==  -2? "not enough memory"  :
            (error==  -3? "reordering problem" :
            (error==  -4? "zero pivot, numerical fact. or iterative refinement problem" :
            (error==  -5? "unclassified (internal) error"     :
            (error==  -6? "preordering failed (matrix types 11, 13 only)" :
            (error==  -7? "diagonal matrix problem"           :
            (error==  -8? "32-bit integer overflow problem"   :
            (error== -10? "no license file pardiso.lic found" :
            (error== -11? "license is expired"                :
            (error== -12? "wrong username or hostname"        :
            (error==-100? "reached maximum number of Krylov-subspace iteration in iterative solver"     :
            (error==-101? "no sufficient convergence in Krylov-subspace iteration within 25 iterations" :
            (error==-102? "error in Krylov-subspace iteration"      :
            (error==-103? "break-down in Krylov-subspace iteration" : "unknown error" ))))))))))))))) << endl;
    throw 42;
  }
  return error;
}


}  // namespace m

