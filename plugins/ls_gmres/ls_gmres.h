#ifndef ls_gmres_h
#define ls_gmres_h

#include "mlinearsystem.h"


namespace m {


// implementation of a linear system solver, using GMRES (double p.)
class ls_gmres : public mlinearsystem< double > {
 public:
  // constructor
  ls_gmres() : mlinearsystem< double >(), c__1(1) {}
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
  int iluk( int *n, double *a, int *ja, int *ia, int *lfil,
            double*& aluold, int*& jluold, int *ju, int*& levsold, int *iwk,
            double *w, int *jw, int *ierr );
  double ddot(int *n, double *dx, int *incx, double *dy, int *incy);
  int daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  double dnrm2(int *n, double *dx, int *incx);
  void amux(int *n, double *x, double *y, double *a, int *ja, int *ia);
  void lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju);
  void pgmres( int *n, int *im, double *rhs, double *sol, double *vv,
               double *eps, int *maxits, int*iout,
               double *aa, int *ja, int *ia,
               double *alu, int *jlu, int *ju,
               int *ierr );
 private:
  // members
 private:
  mmatrix_csr< double,1 > m_A;
  int c__1;
};


}  // namespace m


#endif

