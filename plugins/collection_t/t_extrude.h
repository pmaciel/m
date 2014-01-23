#ifndef t_extrude_h
#define t_extrude_h

#include "mkernel.h"

// extrusion transform module
class t_extrude : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m2, const XMLNode& x);
 private:
  void zsteps_idf(const std::string& str);
  void zsteps_geomp(const std::string& str);
  std::vector< double > m_zsteps;
};

#endif

