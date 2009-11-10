#ifndef mmatrix_h
#define mmatrix_h

#include <vector>


namespace m {


/*
 * this is a struct describing a matrix; since it is based on parent class
 * template extension to simulate virtual functions, it is not a polymorphic
 * type and thus can't be used with factories.
 */

// description of a matrix (T: storage type, P: parent class)
template< typename T, class P >
struct mmatrix {
  // constructor and indexing operators
  mmatrix() : Nr(0), Nc(0), issparse(false) {}
  const T& operator()(const unsigned r, const unsigned c) const { return P::operator()(r,c); };
        T& operator()(const unsigned r, const unsigned c)       { return P::operator()(r,c); };
  void zerorow(const unsigned r)                                {        P::zerorow(r);      };
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Nr, unsigned _Nc)                       { Nr=_Nr; Nc=_Nc; }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {}
  // members
  unsigned Nr;
  unsigned Nc;
  bool issparse;
};


// implementation of a dense matrix, std::vector< std::vector< T > > based
template< typename T >
struct mmatrix_vv : mmatrix< T,mmatrix_vv< T > > {
  typedef mmatrix< T,mmatrix_vv< T > > P;
  // indexing operators
  const T& operator()(const unsigned r, const unsigned c) const { return A[r][c]; }
        T& operator()(const unsigned r, const unsigned c)       { return A[r][c]; }
  void zerorow(const unsigned r) { A[r].assign(P::Nc,T()); }
  // initialize method for dense variation
  void initialize(unsigned _Nr, unsigned _Nc) {
    P::initialize(_Nr,_Nc);
    A.assign(P::Nr,std::vector< T >(P::Nc,T()));
  }
  // members
  std::vector< std::vector< T > > A;
};


// implementation of a dense matrix, T[][] based
template< typename T >
struct mmatrix_aa : mmatrix< T,mmatrix_aa< T > > {
  typedef mmatrix< T,mmatrix_aa< T > > P;
  // destructor and indexing operators
  ~mmatrix_aa() {
    if (P::Nr || P::Nc) {
      delete[] A[0];
      delete[] A;
    }
  }
  const T& operator()(const unsigned r, const unsigned c) const { return A[r][c]; }
        T& operator()(const unsigned r, const unsigned c)       { return A[r][c]; }
  void zerorow(const unsigned r) {
    for (unsigned c=0; c<P::Nc; ++c) A[r][c] = T();
  }
  // initialize method for dense variation
  void initialize(unsigned _Nr, unsigned _Nc) {
    P::initialize(_Nr,_Nc);
    A    = new T*[ P::Nr ];
    A[0] = new T [ P::Nr * P::Nc ];
    for (unsigned r=1; r<P::Nr; ++r)
      A[r] = A[r-1] + P::Nc;
    for (unsigned r=0; r<P::Nr; ++r)
      zerorow(r);
  }
  // members
  T **A;
};


// implementation of a sparse matrix, compressed sparse rows format
// note: BASE={0,1}: {0,1}-based indexing (other values probably don't work)
template< typename T, int BASE >
struct mmatrix_csr : mmatrix< T,mmatrix_csr< T,BASE > > {
  typedef mmatrix< T,mmatrix_csr< T,BASE > > P;
  // cons/destructor and indexing operators
  mmatrix_csr() : P(), zero(T()) { P::issparse = true; }
  ~mmatrix_csr() {
    if (NNZ) {
      delete[] A;
      delete[] JA;
      delete[] IA;
    }
  }
  const T& operator()(const unsigned r, const unsigned c) const { const int i=getindex(r,c); if (i<0) return zero; return A[i]; }
        T& operator()(const unsigned r, const unsigned c)       { const int i=getindex(r,c); if (i<0) return zero; return A[i]; }
  void zerorow(const unsigned r) {
    for (int k=IA[r]-BASE; k<IA[r+1]-BASE; ++k)
      A[k] = T();
  }
  // initialize method for base/sparse variation
  void initialize(unsigned _Nr, unsigned _Nc) { P::initialize(_Nr,_Nc); }
  void initialize(const std::vector< std::vector< unsigned > >& nz) {
    // set base class members
    initialize((unsigned) nz.size(), (unsigned) nz.size());
    // set row indices
    NNU = (int) nz.size();
    IA = new int[NNU+1];
    IA[0] = BASE;
    for (int r=0; r<NNU; ++r)
      IA[r+1] = IA[r] + nz[r].size();
    NNZ = IA[NNU]-BASE;
    // set column indices
    JA = new int[NNZ];
    for (int r=0; r<NNU; ++r) {
      int j = IA[r]-BASE;
      for (unsigned n=0; n<(unsigned) nz[r].size(); ++n)
        JA[j++] = (int) nz[r][n] + BASE;
    }
    // set entries
    A = new T[NNZ];
    for (int i=0; i<NNZ; ++i)
      A[i] = T();
  }
  int getindex(const unsigned r, const unsigned c) const {
    for (int k=IA[r]-BASE; k<IA[r+1]-BASE; ++k)
      if (JA[k]-BASE==(int) c)
        return k;
    return -1;
  }
  // members
  T zero;
  int *IA;
  int *JA;
  T   *A;
  int NNZ;
  int NNU;
};


}  // namespace m


#endif

