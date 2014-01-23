#ifndef t_probe_h
#define t_probe_h

#include "mkernel.h"

// module to check variable values at given coordinates/node index
class t_probe : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
};

#endif

