#ifndef t_simplex_h
#define t_simplex_h

#include "mkernel.h"

// module to transform elements into simplexes
class t_simplex : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
  static std::vector< m::melem > simplex(const m::mtype& t, const std::vector< unsigned >& en);
};

#endif

