#ifndef t_solution_h
#define t_solution_h

#include "mkernel.h"

// Read solution from another file transform module
class t_solution : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m2, const XMLNode& x);
};

#endif

