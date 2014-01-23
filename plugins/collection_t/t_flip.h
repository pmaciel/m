#ifndef t_flip_h
#define t_flip_h

#include "mkernel.h"


// fix negative element volumes by flipping nodes
class t_flip : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
};


#endif

