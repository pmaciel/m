
#include <sstream>
#include <limits>
#include "mfactory.h"
#include "t_probe.h"

using namespace m;


Register< mtransform,t_probe > mt_probe(2,"-tprobe","[real:real[:real]] check variable values at given coordinates",
                                          "",       "[int]              ... or given node index");


namespace aux {
  template< typename T >
  std::vector< T > getVValues(const std::string& l)
  {
      using std::string;
      using std::vector;
      using std::istringstream;

      // split string by ':'
      vector< T > r;
      if (l.length() && l!=":") {
        string::size_type p1 = 0,
                          p2 = 0;
        while (p2!=string::npos) {
          p2 = l.find(":",p1);
          istringstream ss( l.substr(p1,(p2==string::npos? p2:p2-p1)) );
          r.push_back(T());
          ss >> r.back();
          p1 = p2+1;
        }
      }
      return r;
  }
}


void t_probe::transform(GetPot& o, mmesh& m)
{
  using std::cout;
  using std::cerr;
  using std::endl;


  // interpret option (get "tentative" position and node index)
  const std::string option(o.get(o.inc_cursor(),""));
  const std::vector< double   > p(aux::getVValues< double   >(option));
  const std::vector< unsigned > i(aux::getVValues< unsigned >(option));


  // set probing node index and distance to given coordinates
  int idx = -1;
  double dmin(std::numeric_limits< double >::max());
  if ((unsigned) p.size()==m.d()) {

    // search closest node by coordinates proximity
    for (unsigned n=0; n<m.n(); ++n) {
      double d = 0.;
      for (unsigned i=0; i<m.d(); ++i)
        d += (p[i]-m.vv[i][n]) * (p[i]-m.vv[i][n]);
      if (d<dmin) {
        dmin = d;
        idx  = (int) n;
      }
    }
    dmin = sqrt(dmin);
    if (idx<0) {
      cerr << "error: probe node index not found!" << endl;
      throw 42;
    }

  }
  else if ((unsigned) i.size()==1 && i[0]<m.n()) {

    // set from option
    idx = (int) i[0];
    dmin = 0.;

  } else {

    return;

  }


  // display probing
  cout << "info: node: " << idx << "  (d=" << dmin << ')' << endl;
  for (unsigned t=0; t<m.v(); ++t)
    cout << "info: variable \"" << m.vn[t] << "\": " << m.vv[t][idx] << endl;
}

