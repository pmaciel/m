#ifndef f_muphys_h
#define f_muphys_h

#include "mkernel.h"

// module for MuPhyS files
class f_muphys : public m::mfoutput {
 public:
  void write(GetPot& o, const m::mmesh& m, const XMLNode& x);
};

#endif

