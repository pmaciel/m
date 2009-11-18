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
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  void solve();
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv) {
    mlinearsystem< double >::initialize(_Ne,_Nv);
    m_A.initialize(Ne,Nv);
  }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {
    m_A.initialize(nz);
  }
  // indexing functions
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
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

