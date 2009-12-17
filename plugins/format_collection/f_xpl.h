#ifndef f_xpl_h
#define f_xpl_h

#include "mkernel.h"

// module for xplot
class f_xpl : public m::mfinput {
 public:
  void read(GetPot& o, m::mmesh& m);
};

#endif

