#ifndef t_surfmap_h
#define t_surfmap_h

#include "mkernel.h"
#include "mpoint.h"

// module to map a structured grid on a 3D surface mesh
class t_surfmap : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
 private:
  // triangulation utilities
  double trisplit(const std::vector< double >& vu, const double& u1, const double& u2);
  bool tricheck( const double& x, const double& xa, const double& xb, const double& xc,
                 const double& y, const double& ya, const double& yb, const double& yc );
};

#endif

