
#include "mfactory.h"
#include "f_plt.h"
#include "t_solution.h"

using namespace std;
using namespace m;


Register< mtransform,t_solution > mt_solution("-tsol","[str]: Tecplot (nodal) solution file");


void t_solution::transform(GetPot& o, mmesh& m)
{
  if (m.d()!=2 && m.d()!=3)
    return;
  const string fn = o.get(o.inc_cursor(),"");
  cout << "::solution [d]: " << m.d() << "..." << endl;


  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }
  string s;
  while (getline(f,s)) {
    const string key(f_plt::upper(f_plt::splitstring(s)[0]));
    if (key.find("VARIABLES")==0) {

      m.vn = f_plt::getVariables(s);

    }
    else if (key=="ZONE") {

      string   zn = "";
      mtype    zt = ORDERED;
      TecZone tz = f_plt::getZoneHeader(s,zn,zt);
      f_plt::readZoneNodeValues(f,m.vv,tz.i,m.v(),tz.isblock);

    }
  }


  f.close();
  cout << "::solution [d]: " << m.d() << "." << endl;
}


