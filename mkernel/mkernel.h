#ifndef mkernel_h
#define mkernel_h

#include "ext/GetPot.h"
#include "mconfig.h"
#include "mmesh.h"


namespace m {


// transformation plugin
class mtransform {
 public:
  virtual ~mtransform() {}
  virtual void transform(GetPot& o, mmesh& m) {}
};


// file input plugin
class mfinput {
 public:
  virtual ~mfinput() {}
  virtual void read(GetPot& o, mmesh& m) {}
};


// file output plugin
class mfoutput {
 public:
  virtual ~mfoutput() {}
  virtual void write(GetPot& o, const mmesh& m) {}
};


}  // namespace m


#endif

