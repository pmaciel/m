#ifndef t_zmerge_h
#define t_zmerge_h

#include "mkernel.h"

// zone boolean operations transform module
class t_zbool : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
};

#endif

