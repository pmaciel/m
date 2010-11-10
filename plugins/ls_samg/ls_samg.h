#ifndef ls_samg_h
#define ls_samg_h

#include <map>
#include "mlinearsystem.h"
#include "ls_samg_parm.h"


namespace m {


// implementation of a linear system solver, using SAMG (double p.)
class ls_samg : public mlinearsystem< double > {
 public:
  // interfacing functions
  void reset(const double& v=0.) { mlinearsystem< double >::reset(v); m_A.reset(v); }
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  void solve();
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb=1) {
    mlinearsystem< double >::initialize(_Ne,_Nv,_Nb);
    m_A.initialize(Ne,Nv,Nb);
  }
  void initialize(const std::vector< std::vector< unsigned > >& nz);
  // indexing functions (absolute indexing)
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
  // indexing functions (block indexing)
  const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return m_A(R,C,r,c); }
        double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return m_A(R,C,r,c); }
 private:
  // auxiliary functions
  int digit(const int &number, const unsigned &place, unsigned radix=10) const;
  void init_block_iu();
  void init_block_ip();
  void init_maps();
  void init_params();
 private:
  // members:
  // linear system matrix (1-based indexed CSR)
  mmatrix_csr< double,1 > m_A;
  // primary parameters
  struct {

    // configurable
    int matrix, nsolve, ifirst, ncyc, iswtch, idump, iout;
    double eps, a_cmplx, g_cmplx, p_cmplx, w_avrge, chktol;

    // non-configurable (set with other options)
    int nsys, ndiu, ndip;
    std::vector< int > iu, ip, iscale;

    // non-configurable (output)
    int ncyc_done, ierr;
    double res_in, res_out;

  } p;
  // maps for primary and secondary parameters
  std::vector< aux::s_parm > m_prim;
  std::vector< aux::s_parm > m_sdry;
  // map for error descriptions
  std::map< int,std::string > m_ierr;
};


}  // namespace m


#endif

