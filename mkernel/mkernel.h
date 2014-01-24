#ifndef mkernel_h
#define mkernel_h

// imp/export symbols for shared libraries
#ifdef M_CSHLIB
#define M_SYMBOL M_SYMBOL_IMPORT
#else
#define M_SYMBOL M_SYMBOL_EXPORT
#endif


#include <utility>  // (for std::pair)

#include "ext/GetPot.h"  //TODO: remove
#include "ext/xmlParser.h"
#include "mconfig.h"
#include "mmesh.h"


namespace m {


// transformation plugin
class M_SYMBOL mtransform {
 public:
  virtual ~mtransform() {}
  virtual void transform(GetPot& o, mmesh& m, const XMLNode& x) {}
};


// file input plugin
class M_SYMBOL mfinput {
 public:
  virtual ~mfinput() {}
  virtual void read(GetPot& o, mmesh& m, const XMLNode& x) {}
};


// file output plugin
class M_SYMBOL mfoutput {
 public:
  virtual ~mfoutput() {}
  virtual void write(GetPot& o, const mmesh& m, const XMLNode& x) {}
};


// utilities
namespace utils {
  using std::pair;
  using std::string;
  using std::vector;


  // string splitting/cleaning
  vector< string >& split(const string& s, char delim, vector<string>& elems);
  vector< string >  split(const string& s, char delim);
  string trim(const string& s, const string& t=" ");
  string upper(const string& s);
  string get_file_extension(const string& f);


  // internal/commands options as vector/xml string
  vector< pair< string,string > > get_operands(const string& s);
  XMLNode                         get_operands_xml(const string& s);   //FIXME remove


}


}  // namespace m


#endif

