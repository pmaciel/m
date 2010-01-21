#ifndef ls_pardiso_h
#define ls_pardiso_h

#include "mlinearsystem.h"


namespace m {


// implementation of a linear system solver, using Pardiso (double p.)
class ls_pardiso : public mlinearsystem< double > {
 public:
  // constructor
  ls_pardiso();
  // interfacing functions
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  void solve();
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb=1) {
    mlinearsystem< double >::initialize(_Ne,_Nv,_Nb);
    m_A.initialize(Ne,Nv,Nb);
  }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {
    m_A.initialize(nz);
  }
  // indexing functions (absolute indexing)
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
  // indexing functions (block indexing)
  const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return m_A(R,C,r,c); }
        double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return m_A(R,C,r,c); }
  // auxiliary functions
 private:
  int call_pardisoinit();
  int call_pardiso(int phase, int msglvl);
  // members
 private:
  mmatrix_csr< double,1 > m_A;
  void*  pt[64];  // internal memory pointer (void* for both 32/64-bit)
  int    iparm[64];
  double dparm[64];
  int    maxfct,
         mnum,
         mtype,
         nrhs;
};


}  // namespace m


#endif

