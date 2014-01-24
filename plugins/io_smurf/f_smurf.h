#ifndef f_smurf_h
#define f_smurf_h

#include "mkernel.h"

// module for .smurf
class f_smurf : public m::mfinput,
                public m::mfoutput {
 public:
  void read(GetPot& o, m::mmesh& m, const XMLNode& x);
  void write(GetPot& o, const m::mmesh& m, const XMLNode& x);
 private:
  const std::vector< std::vector< unsigned > > convert_from_vtelem(const std::vector< m::melem >& ve1);
  const std::vector< m::melem > convert_to_vtelem(const std::vector< std::vector< unsigned > >& ve1);
};

#endif

