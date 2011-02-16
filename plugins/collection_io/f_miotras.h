#ifndef f_miotras_h
#define f_miotras_h

#include "mkernel.h"

// module for Miotras
class f_miotras : public m::mfinput,
                  public m::mfoutput {
 public:
  void read(GetPot& o, m::mmesh& m);
  void write(GetPot& o, const m::mmesh& m);
 private:
  std::vector< unsigned > findneighbours(const std::vector< m::melem >& ve, const std::vector< unsigned >& n);
};

#endif

