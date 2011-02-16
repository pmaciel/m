#ifndef t_walldistance_h
#define t_walldistance_h

#include "mkernel.h"

// module to calculate wall distance to given boundary zones
class t_walldistance : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
};

#endif

