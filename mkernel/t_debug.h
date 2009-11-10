#ifndef t_debug_h
#define t_debug_h

#include "mkernel.h"

// module to debug
class t_debug : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
};

#endif

