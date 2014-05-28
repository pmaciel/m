#ifndef utils_h
#define utils_h


#include "mmesh.h"


// useful common functions
namespace m {
namespace utils {


  unsigned getvindex(const m::mmesh& m, const std::string& n);
  unsigned getzindex(const m::mmesh& m, const std::string& n);

  std::vector< std::pair< std::string, std::string > > getoperands(const std::string& s);
  std::vector< std::string >& split(const std::string& s, char delim, std::vector< std::string >& elems);
  std::vector< std::string >  split(const std::string& s, char delim);


}  // namespace utils
}  // namespace m


#endif
