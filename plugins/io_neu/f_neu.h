#ifndef f_neu_h
#define f_neu_h

#include "mkernel.h"

// module for .neu
class f_neu : public m::mfinput {
 public:
  void read(GetPot& o, m::mmesh& m, const XMLNode& x);
};

#endif

