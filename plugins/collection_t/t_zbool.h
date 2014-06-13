#ifndef t_zmerge_h
#define t_zmerge_h

#include "mkernel.h"

// zone boolean operations transform module
class t_zbool : public m::mtransform {

public:
  void transform(GetPot& o, m::mmesh& m);

private:

  // zone union (merge)
  void zunion(m::mmesh& m, const std::vector< unsigned >& zidx);

  // zone intersection
  void zinter(m::mmesh& m, const std::vector< unsigned >& zidx);

  // zone splitting, neighboring node-sharing zone based
  void zsplit_sharednodes(m::mmesh& m, const std::vector< unsigned >& zidx);

};

#endif

