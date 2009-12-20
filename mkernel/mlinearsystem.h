#ifndef mlinearsystem_h
#define mlinearsystem_h

#include <vector>
#include "boost/progress.hpp"
#include "mmatrix.h"


namespace m {


/*
 * description of a linear system solver
 * (solution and right-hand side vectors are included, matrix should be
 * included in the derived class)
 */
template< typename T >
class mlinearsystem {
 public:
  // cons/destructor
  mlinearsystem() : issparse(true), Ne(0), Nv(0) {}
  virtual ~mlinearsystem()                       {}
  // interfacing functions
  virtual void zerorow(const unsigned r) = 0;
  virtual void solve()                   = 0;
  void print(std::ostream& o) {
    o << "m::mlinearsystem::X(Nv:" << Nv << "):" << std::endl;
    for (unsigned j=0; j<Nv; ++j)
      o << ' ' << m_X[j] << std::endl;
    o << "m::mlinearsystem::B(Ne:" << Nv << "):" << std::endl;
    for (unsigned i=0; i<Ne; ++i)
      o << ' ' << m_B[i] << std::endl;
    o << "m::mlinearsystem::A(" << (issparse? "true":"false") << "):" << std::endl;
    for (unsigned i=0; i<Ne; ++i) {
      for (unsigned j=0; j<Nv; ++j)
        o << ' ' << A(i,j);
      o << std::endl;
    }
  }
  // initialize methods for dense/sparse variations
  virtual void initialize(unsigned _Ne, unsigned _Nv) {
    Ne = _Ne;
    Nv = _Nv;
    m_X.assign(Ne,T());
    m_B.assign(Nv,T());
  }
  virtual void initialize(const std::vector< std::vector< unsigned > >& nz) = 0;
  // indexing functions
  virtual const T& A(const unsigned r, const unsigned c) const = 0;
  virtual       T& A(const unsigned r, const unsigned c)       = 0;
          const T& X(const unsigned r) const { return m_X[r]; }
                T& X(const unsigned r)       { return m_X[r]; }
          const T& B(const unsigned c) const { return m_B[c]; }
                T& B(const unsigned c)       { return m_B[c]; }
  // members
 public:
  bool issparse;
 protected:
  unsigned Ne;
  unsigned Nv;
  std::vector< T > m_X;
  std::vector< T > m_B;
};


// implementation of a linear system solver, using gaussian elimination (double p.)
class ls_gauss : public mlinearsystem< double > {
 public:
  // constructor
  ls_gauss() : mlinearsystem< double >() { issparse=false; }
  // interfacing functions
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  void solve();
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv) {
    mlinearsystem< double >::initialize(_Ne,_Nv);
    m_A.initialize(Ne,Nv);
  }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {}
  // indexing functions
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
  // members
 private:
  mmatrix_aa< double > m_A;
};


// implementation of a linear system solver, using band LU (double p.)
class ls_bandlu : public mlinearsystem< double > {
 public:
  // constructor
  ls_bandlu() : mlinearsystem< double >(), m_ld(0), m_ud(0) { issparse=true; }
  // interfacing functions
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  void solve();
  // initialize methods for sparse variation
  void initialize(const std::vector< std::vector< unsigned > >& nz);
  // indexing functions
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,m_ld+c-r); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,m_ld+c-r); }
 private:
  // auxiliary functions
  void LUDecomposition();
  // members
 private:
  mmatrix_aa< double > m_A;  // matrix (dense)
  unsigned m_bandwidth;      // bandwidth
  unsigned m_ld;             // number of lower co-diagonals
  unsigned m_ud;             // number of upper co-diagonals
};


}  // namespace m


#endif

