#ifndef f_msh_h
#define f_msh_h

#include "mkernel.h"

// module for .msh
class f_msh : public m::mfinput {
 public:
  void read(GetPot& o, m::mmesh& m);
};

#endif

