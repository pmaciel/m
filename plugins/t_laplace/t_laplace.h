#ifndef t_laplace_h
#define t_laplace_h

#include "mkernel.h"


// module with a 2D Laplace equation solver
class t_laplace : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
 private:
  // mesh utilities
  std::vector< double > getvvalues(const std::string& s);
  std::vector< m::mzone >::const_iterator getzoneit  (const std::string& n, const m::mmesh& m);
  unsigned                                getzoneidx (const std::string& n, const m::mmesh& m);
};

#endif

