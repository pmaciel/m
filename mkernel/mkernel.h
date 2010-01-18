#ifndef mkernel_h
#define mkernel_h

// imp/export symbols for shared libraries
#ifdef M_CSHLIB
#define M_SYMBOL M_SYMBOL_IMPORT
#else
#define M_SYMBOL M_SYMBOL_EXPORT
#endif

#include "ext/GetPot.h"
#include "mconfig.h"
#include "mmesh.h"

namespace m {


// transformation plugin
class M_SYMBOL mtransform {
 public:
  virtual ~mtransform() {}
  virtual void transform(GetPot& o, mmesh& m) {}
};


// file input plugin
class M_SYMBOL mfinput {
 public:
  virtual ~mfinput() {}
  virtual void read(GetPot& o, mmesh& m) {}
};


// file output plugin
class M_SYMBOL mfoutput {
 public:
  virtual ~mfoutput() {}
  virtual void write(GetPot& o, const mmesh& m) {}
};


}  // namespace m


#endif

