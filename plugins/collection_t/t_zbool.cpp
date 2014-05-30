
#include <sstream>

#include "mfactory.h"
#include "t_zbool.h"


using namespace m;


Register< mtransform,t_zbool > mt_zbool( 2,
  "-tzunion", "[str:...] zone union (merging of zones of the same type)",
  "-tzinter", "[str:...] zone intersection" );


void t_zbool::transform(GetPot& o, mmesh& m)
{
  using namespace std;


#if 0
  const string k = o[o.get_cursor()],
               v = (k=="-tvsort"||k=="-tvaxiz"||k=="-tzsort"? "" : o.get(o.inc_cursor(),""));

  // operations that apply in one shot
       if (k=="-tvsort") { vsort(m);   return; }
  else if (k=="-tzsort") { zsort(m);   return; }
#endif


  cout << "::zone union..." << endl;


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
    z.n.append("_union_" + m.vz[ zm[i] ].n);
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


  cout << "::zone union." << endl;
}

