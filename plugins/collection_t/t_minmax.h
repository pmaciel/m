#ifndef t_minmax_h
#define t_minmax_h

#include "mkernel.h"

// module to check minimum and maximum variable values at given zones
class t_minmax : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
};

#endif

