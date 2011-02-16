
#include "mfactory.h"
#include "../io_smurf/f_plt.h"
#include "t_avgsolution.h"

using namespace std;
using namespace m;


Register< mtransform,t_avgsolution > mt_avgsolution(2,
  "-tavgsolution","[str] Tecplot (nodal) solution file, to average with current",
  "-tavg",        "[str] ...");


void t_avgsolution::transform(GetPot& o, mmesh& m1)
{
  const unsigned D = m1.d();
  const unsigned V = m1.v();
  const unsigned N = m1.n();
  if (V<=D) {
    cerr << "::avgsolution: nothing to do, exiting" << endl;
  }
  ++Navg;
  const double dNavg = (double) Navg;


  cout << "::avgsolution: read..." << endl;
  const string fn(o.get(o.inc_cursor(),""));
  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }
  mmesh m2;
  string s;
  while (getline(f,s)) {
    const string key(f_plt::upper(f_plt::splitstring(s)[0]));
    if (key.find("VARIABLES")==0) {

      m2.vn = f_plt::getVariables(s);
      if (V!=m2.v()) {
        cerr << "error: there are " << V << " variables, file has " << m2.v() << " instead" << endl;
        throw 42;
      }

    }
    else if (key=="ZONE") {

      string zn;
      mtype zt;
      TecZone tz = f_plt::getZoneHeader(s,zn,zt);
      f_plt::readZoneNodeValues(f,m2.vv,tz.i,m2.v(),tz.isblock);
      if (N!=m2.n()) {
        cerr << "error: there are " << N << " nodes, file has " << m2.n() << " instead" << endl;
        throw 42;
      }

    }
  }
  f.close();
  cout << "::avgsolution: read." << endl;


  cout << "::avgsolution: average..." << endl;
  for (unsigned i=D; i<V; ++i)
    for (unsigned j=0; j<N; ++j)
      m1.vv[i][j] = (m1.vv[i][j] * dNavg + m2.vv[i][j])/(dNavg+1.);
  cout << "::avgsolution: average." << endl;
}


