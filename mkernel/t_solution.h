#ifndef t_solution_h
#define t_solution_h

#include "mkernel.h"

// Solution from Tecplot file transform module
class t_solution : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m2);
};

#endif

