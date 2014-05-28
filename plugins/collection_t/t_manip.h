#ifndef t_manip_h
#define t_manip_h

#include "mkernel.h"

// manipulation of variables and zones module
class t_manip : public m::mtransform {

 public:
  void transform(GetPot& o, m::mmesh& m);

 private:
  // variables operations, per se
  void vsort(m::mmesh& m);
  void vkeep(m::mmesh& m, const std::string& s);
  void vrm(m::mmesh& m, const unsigned i);
  void vmv(m::mmesh& m, const unsigned i, const unsigned j);
  void vren(m::mmesh& m, const unsigned i, const std::string& n);
  void vadd(m::mmesh& m, const std::string& n, const std::string& f, const std::vector< unsigned > zindex=std::vector< unsigned >());
  void vaxiz(m::mmesh& m);

  // zone operations, per se
  void zsort(m::mmesh& m);
  void zkeep(m::mmesh& m, const std::string& s);
  void zrm(m::mmesh& m, const unsigned i);
  void zmv(m::mmesh& m, const unsigned i, const unsigned j);
  void zren(m::mmesh& m, const unsigned i, const std::string& n);

};

#endif


