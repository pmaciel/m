#ifndef t_manip_h
#define t_manip_h

#include "mkernel.h"

// manipulation of variables and zones module
class t_manip : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);

 private:
  // variables operations, per se
  void vrm(m::mmesh& m, const unsigned i);
  void vmv(m::mmesh& m, const unsigned i, const unsigned j);
  void vren(m::mmesh& m, const unsigned i, const std::string& n);
  void vadd(m::mmesh& m, const std::string& n, const std::string& f);

  // zone operations, per se
  void zrm(m::mmesh& m, const unsigned i);
  void zmv(m::mmesh& m, const unsigned i, const unsigned j);
  void zren(m::mmesh& m, const unsigned i, const std::string& n);

  // utilities
  unsigned getvindex(const m::mmesh& m, const std::string& n);
  unsigned getzindex(const m::mmesh& m, const std::string& n);
  std::vector< std::pair< std::string,std::string > > getoperations(const std::string& s);
};

#endif


