#ifndef t_zmerge_h
#define t_zmerge_h

#include "mkernel.h"

// zone merge transform module
class t_zmerge : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
};

#endif

