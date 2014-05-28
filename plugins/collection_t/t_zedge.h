#ifndef t_zedge_h
#define t_zedge_h

#include "mkernel.h"

// zone edge extraction transform module
class t_zedge : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
};

#endif

