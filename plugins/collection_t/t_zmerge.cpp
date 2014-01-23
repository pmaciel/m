
#include "mfactory.h"
#include "t_zmerge.h"
#include <sstream>

using namespace std;
using namespace m;


Register< mtransform,t_zmerge > mt_zmerge("-tzm","[str:...] zone merge, for zones of the same type");


void t_zmerge::transform(GetPot& o, mmesh& m, const XMLNode& x)
{
  cout << "::zone merge..." << endl;


  // find zones, asserting same type
  vector< unsigned > zm;  // where old zones are
  string str = o.get(o.inc_cursor(),"");
  replace(str.begin(),str.end(),':',' ');
  istringstream is(str);
  string zname;
  while (is >> zname) {
    for (unsigned j=0; j!=m.vz.size(); ++j)
      if (m.vz[j].n==zname && (zm.size()? m.vz[j].t==m.vz[ zm[0] ].t:true)) {
        zm.push_back(j);
        break;
      }
  }
  if (zm.size()<2) {
    cerr << "error: incorrect zone names definition (different types, missing ':'?)" << endl;
    throw 42;
  }
  mzone& z = m.vz[ zm[0] ];  // this is the zone that merges the others


  // reserve new connectivity size
  cout << "info [e]: \"" << z.n << "\" [" << z.e2n.size() << "]";
  const unsigned Nelemi = z.e2n.size();  // original nb. elements
  unsigned Nelemf = Nelemi;              // count new elements
  for (unsigned i=1; i<zm.size(); ++i) {
    z.n.append("_m_" + m.vz[ zm[i] ].n);
    Nelemf += m.vz[ zm[i] ].e2n.size();
    cout << ", \"" << m.vz[ zm[i] ].n << "\" [" << m.vz[ zm[i] ].e2n.size() << "]";
  }
  cout << "..." << endl;
  z.e2n.resize(Nelemf);


  // do an hostile take-over
  for (unsigned i=1, k=Nelemi; i<zm.size(); ++i) {
    for (unsigned j=0; j<m.vz[ zm[i] ].e2n.size(); ++j) {
      vector< unsigned >& en1 = z.e2n[k++].n;
      vector< unsigned >& en2 = m.vz[ zm[i] ].e2n[j].n;
      en1.swap(en2);
    }
  }


  // delete useless zones
  for (unsigned i=zm.size()-1; i>=1; --i)
    m.vz.erase( m.vz.begin()+zm[i] );
  cout << "info [e]: \"" << z.n << "\" [" << z.e2n.size() << "]." << endl;


  cout << "::zone merge." << endl;
}

