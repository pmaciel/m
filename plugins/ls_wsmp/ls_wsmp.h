#ifndef ls_wsmp_h
#define ls_wsmp_h

#include "mlinearsystem.h"


namespace m {


// implementation of a linear system solver, using WSMP (double p.)
class ls_wsmp : public mlinearsystem< double > {
 public:
  // constructor
  ls_wsmp();
  // interfacing functions
  void reset(const double& v=0.) { mlinearsystem< double >::reset(v); m_A.reset(v); }
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
 private:
  // auxiliary functions
  int call_wsmp(int task);
  // members
 private:
  mmatrix_csr< double,0 > m_A;
  int    iparm[64];
  double dparm[64];
};


}  // namespace m


#endif

