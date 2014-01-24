#ifndef f_cgns_h
#define f_cgns_h

#include "mkernel.h"


// module for .cgns
class f_cgns : public m::mfinput,
               public m::mfoutput {
 public:
  void read(GetPot& o, m::mmesh& m, const XMLNode& x);
  void write(GetPot& o, const m::mmesh& m, const XMLNode& x);
};


#endif

