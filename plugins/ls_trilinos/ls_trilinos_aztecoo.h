#ifndef ls_trilinos_aztecoo_h
#define ls_trilinos_aztecoo_h

#include "Epetra_CrsMatrix.h"
#include "AztecOO.h"

#include "mlinearsystem.h"


namespace m {


// implementation of a linear system solver, using Trilinos AztecOO (double p.)
class ls_trilinos_aztecoo : public mlinearsystem< double > {
 public:
  // cons/destructor
  ls_trilinos_aztecoo();
  ~ls_trilinos_aztecoo();
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
  // members
 private:
  bool issetup;
  Epetra_CrsMatrix *myA;
  Epetra_Vector    *myX;
  Epetra_Vector    *myB;
  AztecOO m_solver;
  double az_params[AZ_PARAMS_SIZE];
  int    az_options[AZ_OPTIONS_SIZE];
  mmatrix_csr< double,0 > m_A;
};


}  // namespace m


#endif

