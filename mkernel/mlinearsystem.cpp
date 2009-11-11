
#include <cmath>
#include <algorithm>
#include "mfactory.h"
#include "mlinearsystem.h"


using namespace m;

Register< mlinearsystem< double >,ls_gauss > __ls_gauss("ls_gauss","linear system solver, using gaussian elimination (double p.)");
Register< mlinearsystem< double >,ls_bandlu > __ls_bandlu("ls_bandlu","linear system solver, using band LU (double p.)");


namespace m {


void ls_gauss::solve()
{
  const unsigned N = Ne;
  m_X = m_B;

  double C;
  boost::progress_display pbar(N-1);
  for (unsigned m=0; m<N-1; ++m, ++pbar) {

    // put row with highest diagonal element on top
    C = A(m,m);
    for (unsigned n=m+1; n<N; n++) {
      if (abs(A(n,m)) > abs(C)) {
        for (unsigned p=m; p<N; p++) {
          C = A(m,p);
          A(m,p) = A(n,p);
          A(n,p) = C;
        }
        C    = X(m);
        X(m) = X(n);
        X(n) = C;
        C    = A(m,m);
      }
    }

    // check if diagonal element is (close to) zero
    if (std::abs(C)<1.e-32) {
      std::cerr << "error: matrix is singular (line:" << m << ",C:" << std::abs(C) << ")" << std::endl;
      throw 42;
    }

    // normalize row m
    for (unsigned n=m+1; n<N; n++)
      A(m,n) /= C;
    X(m) /= C;

    // subtract row m from subsequent rows
    for (unsigned n=m+1; n<N; n++) {
      C = A(n,m);
      for (unsigned p=m+1; p<N; p++)
        A(n,p) -= C*A(m,p);
      X(n) -= C*X(m);
    }
  }

  // solve by back substitution
  X(N-1) /= A(N-1,N-1);
  for (unsigned p=0; p<N-1; p++) {
    unsigned m = N-p-2;
    for (unsigned n=m+1; n<N; n++)
      X(m) -= A(m,n)*X(n);
  }
}


void ls_bandlu::solve()
{
  // check if matrix is not already in upper triangular form
  if (m_ld>0)
    LUDecomposition();

  // back substitution
  for (int i=Ne-1; i>=0; --i) {
    const unsigned kend = std::min(m_ud+1,Ne-i);
    for (unsigned k=1; k<kend; ++k)
      B(i) -= m_A(i,k+m_ld) * B(i+k);
    B(i) /= m_A(i,m_ld);
  }

  // swap solution from B vector to X
  m_X.swap(m_B);
}


void ls_bandlu::initialize(const std::vector< std::vector< unsigned > >& nz)
{
  // determine number of lower and upper co-diagonals
  for (unsigned i=0; i<Ne; ++i) {
    const unsigned min = *std::min_element(nz[i].begin(),nz[i].end());
    const unsigned max = *std::max_element(nz[i].begin(),nz[i].end());
    std::max(m_ld,max-i);
    std::max(m_ud,i-min);
  }
  m_bandwidth = 1 + m_ld + m_ud;

  // construct matrix
  m_A.initialize(Ne,m_bandwidth);
}


/*
 *  LUDecomposition factors a condensed banded matrix
 *  (uses Gauss algorithm without column pivot search)
 */
void ls_bandlu::LUDecomposition()
{
  // loop over all rows
  boost::progress_display pbar(Ne-1);
  for (unsigned i=0; i<Ne-1; ++i, ++pbar) {
    const unsigned kend  = std::min(m_ld+1,Ne-i);
    const unsigned kjend = std::min(m_ud+1,Ne-i);

    // check if matrix is singular
    if (std::abs(A(i,m_ld))<1.e-32) {
      std::cerr << "error: matrix is singular (line:" << i << ")" << std::endl;
      throw 42;
    }

    // loop over all rows below row i
    for (unsigned k=1; k<kend; ++k) {
      A(k+i,m_ld-k) /= A(i,m_ld);
      B(k+i) -= A(k+i,m_ld-k) * B(i);

      // loop over columns
      for (unsigned j=1; j!=kjend; ++j) {
        const unsigned jk = j + m_ld - k;
        const unsigned jm = j + m_ld;
        A(k+i,jk) -= A(k+i,m_ld-k) * A(i,jm);
      }
    }
  }
}


}  // namespace m
