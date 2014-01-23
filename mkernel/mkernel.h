#ifndef mkernel_h
#define mkernel_h

// imp/export symbols for shared libraries
#ifdef M_CSHLIB
#define M_SYMBOL M_SYMBOL_IMPORT
#else
#define M_SYMBOL M_SYMBOL_EXPORT
#endif


#include <string>
#include <utility>  // (for std::pair)
#include <vector>

#include "ext/GetPot.h"  //TODO: remove
#include "ext/xmlParser.h"
#include "mconfig.h"
#include "mmesh.h"


namespace m {


// transformation plugin
class M_SYMBOL mtransform {
 public:
  virtual ~mtransform() {}
  virtual void transform(GetPot& o, mmesh& m /*, const XMLNode& x */) {}
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


// utilities
namespace utils {


  // std::string splitting/cleaning
  std::vector< std::pair< std::string,std::string > > getoperands(const std::string& s);
  std::vector< std::string >& split(const std::string& s, char delim, std::vector<std::string>& elems);
  std::vector< std::string >  split(const std::string& s, char delim);
  std::string trim(const std::string& s, const std::string& t=" ");
  std::string upper(const std::string& s);


}


}  // namespace m


#endif

