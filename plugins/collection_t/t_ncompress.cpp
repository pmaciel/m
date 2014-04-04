
#include "mfactory.h"
#include "t_ncompress.h"

using namespace std;
using namespace m;


Register< mtransform,t_ncompress > mt_ncompress("-tncompress","remove unconnected nodes");


void t_ncompress::transform(GetPot& o, mmesh& m)
{
  const unsigned nold = m.n();
  m.compress();
  cout << "info: n=" << nold << '>' << m.n() << endl;
}

