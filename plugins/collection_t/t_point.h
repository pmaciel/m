#ifndef t_point_h
#define t_point_h

#include "mkernel.h"

// point insertion transform module
class t_point : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
};

#endif


