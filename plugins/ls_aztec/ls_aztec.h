#ifndef ls_aztec_h
#define ls_aztec_h

#include "mlinearsystem.h"

/**/
namespace m {
// implementation of a sparse matrix, MBR format
template< typename T >
struct mmatrix_mbr : mmatrix< T,mmatrix_mbr< T > > {
  typedef mmatrix< T,mmatrix_mbr< T > > P;
  // cons/destructor and indexing operators
  mmatrix_mbr() : P(), zero(T()) { P::issparse = true; }
  ~mmatrix_mbr() {
    if (NNZ) {
      delete[] VAL;
      delete[] BINDX;
    }
  }
  const T& operator()(const unsigned r, const unsigned c) const { const int i=getindex(r,c); if (i<0) return zero; return VAL[i]; }
        T& operator()(const unsigned r, const unsigned c)       { const int i=getindex(r,c); if (i<0) return zero; return VAL[i]; }
  void zerorow(const unsigned r) {
    for (int k=BINDX[r]; k<BINDX[r+1]; ++k)
      VAL[k] = T();
  }
  // initialize method for base/sparse variation
  void initialize(unsigned _Nr, unsigned _Nc) { P::initialize(_Nr,_Nc); }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {
    // set number of rows/non-zero entries and allocate data structure
    NNU = (int) nz.size();
    NNZ = 0;
    for (int i=0; i<NNU; ++i)
      NNZ += (int) nz[i].size();
    BINDX = new int[NNZ + 1];

    // set row number of off-diagonal non-zeros and their entries
    BINDX[0] = NNU + 1;
    for (int i=0; i<NNU; ++i) {
      BINDX[i+1] = BINDX[i] + (int) nz[i].size() - 1;
      int k = BINDX[i];
      for (int j=0; j<(int) nz[i].size(); ++j)
        if ((int) nz[i][j]!=i)
          BINDX[k++] = nz[i][j];
    }

    // set entries
    VAL = new T[BINDX[NNU]/*NNZ*/];
    for (int i=0; i<BINDX[NNU]/*NNZ*/; ++i)
      VAL[i] = T();
  }
  // utilities
  int getindex(const unsigned r, const unsigned c) const {
    if (r==c)
      return VAL[r];
    for (int k=BINDX[r]; k<BINDX[r+1]; ++k)
      if (BINDX[k]==(int) c)
        return k;
    return -1;
  }
  // members
  T   zero;
  T   *VAL;
  int *BINDX;
  int NNU;
  int NNZ;
};


}  // namespace m
/**/
namespace m {


// implementation of a linear system solver, using Aztec (double p.)
class ls_aztec : public mlinearsystem< double > {
 public:
  // cons/destructor
  ls_aztec();
  ~ls_aztec();
  // interfacing functions
  void zerorow(const unsigned r) { B(r)=0.; m_A.zerorow(r); }
  void solve();
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv) {
    mlinearsystem< double >::initialize(_Ne,_Nv);
    m_A.initialize(Ne,Nv);
  }
  void initialize(const std::vector< std::vector< unsigned > >& nz);
  // indexing functions
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
 private:
  // auxiliary functions
  void setup();
 private:
  // members
  mmatrix_mbr< double > m_A;
  int *proc_config;  // for communication? (set by AZ_set_proc_config)
  int *data_org;     // for communication? (allocated and set by AZ_transform)
  int *options;      // options
  double *params;    // parameters
  double *status;    // result of AZ_solve
};


}  // namespace m


#endif

