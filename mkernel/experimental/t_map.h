#ifndef t_map_h
#define t_map_h

#include "m2m.h"
#include "mpoint.h"

// module to  map 3D surface mesh onto another
class t_map : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);

 private:
  // parametrization utilities
  std::pair< m::mpoint,m::mpoint > findpoints(const m::mzone& z, const std::vector< std::vector< double > >& vv);
  m::mpoint findnormal(const m::mzone& z, const std::vector< std::vector< double > >& vv);

 private:
  // utilities
  void dump2d(std::ostream& out, const m::mzone& z, const std::vector< std::vector< double > >& vv);
  void dump3d(std::ostream& out, const m::mzone& z, const std::vector< std::vector< double > >& vv);
  void dumpuv(std::ostream& out, const m::mzone& z, const std::vector< std::vector< double > >& vv, const m::mpoint& c, const m::mpoint& u, const m::mpoint& v);
};

#endif

