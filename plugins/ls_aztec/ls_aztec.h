#ifndef ls_aztec_h
#define ls_aztec_h

#include "mlinearsystem.h"


namespace m {


// implementation of a linear system solver, using Aztec (double p.)
class ls_aztec : public mlinearsystem< double > {
 public:
  // cons/destructor
  ls_aztec();
  ~ls_aztec();
  // interfacing functions
  void reset(const double& v=0.);
  void zerorow(const unsigned r);
  void zerorow(const unsigned R, const unsigned r);
  void print(std::ostream& o, bool pmatrix=false);
  void solve();
  // initialize methods for dense/sparse variations
  void initialize(unsigned _Ne, unsigned _Nv, unsigned _Nb);
  void initialize(const std::vector< std::vector< unsigned > >& nz);
  // indexing functions (absolute indexing)
  const double& A(const unsigned r, const unsigned c) const;
        double& A(const unsigned r, const unsigned c);
  // indexing functions (block indexing)
  const double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c) const;
        double& A(const unsigned R, const unsigned C, const unsigned r, const unsigned c);
 private:
  // utilities
  void set_az_options(XMLNode& xml, int *options, double *params);
 private:
  // members
  int mtype;                      // matrix type
  mmatrix_msr< double > m_A_msr;  // ... in MSR format
  mmatrix_vbr< double > m_A_vbr;  // ... in VBR format
  int *options;                   // solver options
  double *params;                 // ... parameters
  double *status;                 // ... status
  int *proc_config, *data_org,    // for communication?
      *update,   *update_index,   // update/external vectors
      *external, *extern_index;   // ...
};


}  // namespace m


#endif

