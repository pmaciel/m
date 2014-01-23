#ifndef t_fix_h
#define t_fix_h

#include "mkernel.h"

// module to fix a mesh
class t_fix : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
 private:
  template< int D, int T >
  void check(const m::mzone& z, const std::vector< std::string >& vn, const std::vector< std::vector< double > >& vv);
};


#endif

