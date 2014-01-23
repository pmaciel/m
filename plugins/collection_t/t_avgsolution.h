#ifndef t_avgsolution_h
#define t_avgsolution_h

#include "mkernel.h"

// averaging of solution field transform module
class t_avgsolution : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m1, const XMLNode& x);
 private:
  static unsigned Navg;
};

unsigned t_avgsolution::Navg = 0;

#endif


