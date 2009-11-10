#ifndef f_steven_h
#define f_steven_h

#include "mkernel.h"

// module for Steven files
class f_steven : public m::mfoutput {
 public:
  void write(GetPot& o, const m::mmesh& m);
};

#endif

