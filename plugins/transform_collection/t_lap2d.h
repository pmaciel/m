#ifndef t_lap2d_h
#define t_lap2d_h

#include "mkernel.h"


// module with a 2d laplace equation solver
class t_lap2d : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
 private:
  // mesh utilities
  std::vector< double > getvvalues(const std::string& s);
  std::vector< m::mzone >::const_iterator getzoneit  (const std::string& n, const m::mmesh& m);
  unsigned                                getzoneidx (const std::string& n, const m::mmesh& m);
};

#endif

