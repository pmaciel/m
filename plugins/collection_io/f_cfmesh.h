#ifndef f_cfmesh_h
#define f_cfmesh_h

#include "mkernel.h"

// module for .CFmesh
class f_cfmesh : public m::mfinput,
                 public m::mfoutput {
 public:
  void read(GetPot& o, m::mmesh& m, const XMLNode& x);
  void write(GetPot& o, const m::mmesh& m, const XMLNode& x);
 private:
  // reading functions
  void setvarnames(std::vector< std::string >& vn, unsigned Ndim, unsigned Neqs);
  void readlnodes(std::ifstream& f, m::mmesh& m, unsigned Ndim, unsigned Neqs, unsigned N);
  void readlstates1(std::ifstream& f, m::mmesh& m, unsigned Ndim, unsigned Neqs, unsigned N);
  void readlstates0(std::ifstream& f, m::mmesh& m, unsigned Ndim, unsigned Neqs, unsigned N);
  void readlelems(std::ifstream& f, std::vector< m::melem >& e2n, unsigned Nnodes, unsigned N);
};

#endif

