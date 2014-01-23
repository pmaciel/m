#ifndef t_ncompress_h
#define t_ncompress_h

#include "mkernel.h"

// module to remove unconnected nodes
class t_ncompress : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
};

#endif

