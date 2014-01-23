#ifndef t_info_h
#define t_info_h

#include "mkernel.h"

// module to get information about the current mesh
class t_info : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
};

#endif

