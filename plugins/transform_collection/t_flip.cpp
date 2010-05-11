
#include <list>
#include <sstream>
#include "boost/progress.hpp"
#include "mfactory.h"
#include "t_flip.h"
#include "t_geo.h"

using namespace std;
using namespace m;


Register< mtransform,t_flip > mt_flip("-tflip","fix negative element volumes by flipping nodes");


void t_flip::transform(GetPot& o, mmesh& m)
{
  const unsigned Ndim   = m.d();
  if (!Ndim || !m.v() || !m.n())
    return;

  // set statistics variables and progress bar
  unsigned Nelem   = 0,
           Nerrors = 0,
           Ncorr = 0,
           Nswap = 0;
  for (unsigned i=0; i<m.z(); ++i)
    Nelem += m.e(i);
  boost::progress_display pbar(Nelem);

  // for all zones and all elements, check (and flip)
  for (unsigned i=0; i<m.z(); ++i) {
    const unsigned Nenodes = (m.e(i)? m.vz[i].e2n[0].n.size():0);
    vector< vector< double > > c(Nenodes,vector< double >(Ndim));
    for (unsigned j=0; j<m.e(i); ++j, ++pbar) {
      vector< unsigned >& enodes = m.vz[i].e2n[j].n;

      // check if ok
      for (unsigned i=0; i<Nenodes; ++i)
        for (unsigned j=0; j<Ndim; ++j)
          c[i][j] = m.vv[j][ enodes[i] ];
      if (t_geo::cellgeom_ok(c))
        continue;

      // not ok, flip until ok
      ++Nerrors;
      for (unsigned l=0; l<Nenodes; ++l) {
        ++Nswap;
        swap( enodes[l], enodes[(l+1)%Nenodes] );
        for (unsigned i=0; i<Nenodes; ++i)
          for (unsigned j=0; j<Ndim; ++j)
            c[i][j] = m.vv[j][ enodes[i] ];
        if (t_geo::cellgeom_ok(c)) {
          ++Ncorr;
          break;
        }
      }

    }
  }
  cout << "info: elements/errors: " << Nelem << '/' << Nerrors << endl;
  cout << "info: swaps/corrections: " << Nswap << '/' << Ncorr << endl;
}


