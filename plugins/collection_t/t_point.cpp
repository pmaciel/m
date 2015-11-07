
#include <algorithm>
#include <sstream>
#include "mfactory.h"

#include "t_point.h"

using namespace std;
using namespace m;


Register< mtransform,t_point > mt_point("-tpoint","[d[:d[:d]]] insert a x[,y[,z]] point");


void t_point::transform(GetPot& o, mmesh& m)
{
  // get intended point
  double x[3] = {0.,0.,0.};

  string str = o[o.inc_cursor()];
  const unsigned ndim = std::min<unsigned>(3,1+count(str.begin(),str.end(),':'));
  replace(str.begin(),str.end(),':',' ');
  istringstream(str) >> x[0] >> x[1] >> x[2];


  // get dimensions:
  // - if nothing is there yet, in the beggining there was point
  // - if something is there already, conform (or else)
  const unsigned d = m.d();
  if (!d || d==ndim) {

    // possibly append new variable names
    if (!d) {
      m.vn = vector<string>(ndim,"x");
      for (unsigned i=1; i<ndim; ++i)
        m.vn[i] = char('x'+i);
    }

    // append new point coordinates
    const unsigned v = m.v();
    if (!m.vv.size())
      m.vv.assign(v,vector<double>());
    for (unsigned i=0; i<v; ++i)
      m.vv[i].push_back(i<ndim? x[i] : 0.);

    // append new point zone (getting to a unique name is a pain)
    vector<string> badnames;
    for (vector<mzone>::const_iterator iz=m.vz.begin(); iz!=m.vz.end(); ++iz)
      badnames.push_back(iz->n);
    unsigned c = 1;
    string n = "point";
    while (find(badnames.begin(),badnames.end(),n) != badnames.end()) {
      stringstream ss;
      ss << "point_" << (c++);
      n = ss.str();
    }
    m.vz.push_back(mzone(n));
    m.vz.back().e2n.push_back(melem(vector<unsigned>(1,m.n()-1)));

  }
  else {
    cerr << "error: d=" << d << " != #point=" << ndim << endl;
    throw 42;
  }
}


