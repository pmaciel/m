#ifndef ls_aztec_h
#define ls_aztec_h

#include "mlinearsystem.h"

/**/
namespace m {


// implementation of a sparse matrix, VBR format
template< typename T >
struct mmatrix_vbr : mmatrix< T,mmatrix_vbr< T > > {
  typedef mmatrix< T,mmatrix_vbr< T > > P;
  // cons/destructor and indexing operators
  mmatrix_vbr() : P(), zero(T()) { P::issparse = true; }
  ~mmatrix_vbr() {
    if (bnnz) {
      delete[] val;
      delete[] indx;
      delete[] bindx;
      delete[] rpntr;
      delete[] bpntr;
    }
  }
  const T& operator()(const unsigned r, const unsigned c) const { const int i=getindex(r,c); if (i<0) return zero; return val[i]; }
        T& operator()(const unsigned r, const unsigned c)       { const int i=getindex(r,c); if (i<0) return zero; return val[i]; }
  void zerorow(const unsigned r) {
/*
    for (int k=bindx[r]; k<bindx[r+1]; ++k)
      val[k] = T();
*/
  }
  // initialize method for base/sparse variation
  void initialize(unsigned _Nr, unsigned _Nc, unsigned _Nb=1) { P::initialize(_Nr,_Nc,_Nb); }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {

    bnnu = (int) nz.size();
    bnnz = 0;
    for (int i=0; i<bnnu; ++i)
      bnnz += (int) nz[i].size();

    const unsigned Nblock = 2;

    // allocate data structure (cpntr=rpntr in a struct. symmetric matrix)
    cpntr = rpntr = new int[bnnu+1];
    bpntr = new int[bnnu+1];
    bindx = new int[bnnz];
    indx  = new int[bnnz+1];

    // set it up (diagonal blocks come first)
    rpntr[0] = bpntr[0] = bindx[0] = indx[0] = 0;
    for (int i=0; i<bnnu; ++i) {
      rpntr[i+1] = rpntr[i] + (int) Nblock;
      bpntr[i+1] = bpntr[i] + (int) nz[i].size();
      int k = bpntr[i];
      bindx[k++] = i;
      for (int j=0; j<(int) nz[i].size(); ++j)
        if (nz[i][j]!=i)
          bindx[k++] = nz[i][j];
      for (int j=0, k=bpntr[i]; j<(int) nz[i].size(); ++j, ++k)
        indx[k+1] = indx[k] + (int) (Nblock*Nblock);
    }

    // set entries
    val = new double[indx[bpntr[bnnu]]];  // or bnnz*Nblock^2+1
    for (int i=0; i<indx[bpntr[bnnu]]; ++i)
      val[i] = T();
  }
  // utilities
  int getindex(const unsigned r, const unsigned c) const {
    for (int i=bpntr[r]; i<bpntr[r+1]; ++i)
      if (bindx[i]==c)
        return indx[i];
    return -1;
  }
  // members
  T   zero;
  T   *val;
  int *indx, *bindx, *rpntr, *cpntr, *bpntr;
  int bnnu;
  int bnnz;
};


}  // namespace m
/**/
class XMLNode;
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
  void initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb);
  void initialize(const std::vector< std::vector< unsigned > >& nz);
  // indexing functions (absolute indexing)
  const double& A(const unsigned r, const unsigned c) const { return m_A(r,c); }
        double& A(const unsigned r, const unsigned c)       { return m_A(r,c); }
  // indexing functions (block indexing)
  const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const { return m_A(R,C,r,c); }
        double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c)       { return m_A(R,C,r,c); }
 private:
  // utilities
  void set_az_options(XMLNode& xml, int *options, double *params);
 private:
  // members
  mmatrix_msr< double > m_A;     // linear system matrix
  int *options;                  // solver options
  double *params;                // ... parameters
  double *status;                // ... status
  int *proc_config, *data_org,   // for communication?
      *update,   *update_index,  // update/external vectors
      *external, *extern_index;  // ...
};


}  // namespace m


#endif

