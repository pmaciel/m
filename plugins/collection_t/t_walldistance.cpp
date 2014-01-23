
#include <sstream>
#include "boost/progress.hpp"
#include "ext/Vec.h"
#include "mfactory.h"
#include "t_walldistance.h"

using namespace std;
using namespace m;


Register< mtransform,t_walldistance > mt_walldistance("-twd","[str:...] calculate distance to given boundary (wall) zones");


void t_walldistance::transform(GetPot& o, mmesh& m, const XMLNode& x)
{
  cout << "::walldistance..." << endl;

  const unsigned Ndim = m.d();
  if (Ndim!=2 && Ndim!=3) {
    cerr << "error: dimension can be either 2 or 3!" << endl;
    throw 42;
  }


  cout << "info: mark nodes to calculate distance to..." << endl;
  bool wallnodes_present = false;
  vector< unsigned > wallnodes(m.n(),0);
  string str = o.get(o.inc_cursor(),"");
  replace(str.begin(),str.end(),':',' ');
  istringstream is(str);
  string zname;
  while (is >> zname) {
    for (unsigned j=0; j<m.vz.size(); ++j)
      if (m.vz[j].n==zname) {
        for (unsigned e=0; e<m.e(j); ++e) {
          const vector< unsigned >& en = m.vz[j].e2n[e].n;
          for (vector< unsigned >::const_iterator n=en.begin(); n!=en.end(); ++n)
            wallnodes[*n] = 1;
        }
        wallnodes_present = true;
        break;
      }
  }
  if (!wallnodes_present) {
    cerr << "error: no wall nodes present, specify at least one wall" << endl;
    throw 42;
  }
  cout << "info: mark nodes to calculate distance to." << endl;


  cout << "info: substitute/append \"wd\" variable..." << endl;
  for (unsigned i=0; i<m.v(); ++i)
    if (m.vn[i]=="wd") {
      m.vn.erase(m.vn.begin()+i);
      m.vv.erase(m.vv.begin()+i);
      --i;  // make sure no variables are skipped
    }
  m.vn.insert(m.vn.end(),1,"wd");
  m.vv.insert(m.vv.end(),1,vector< double >(m.n(),1.e99));
  cout << "info: substitute/append \"wd\" variable." << endl;


  cout << "info: calculate distances..." << endl;
  boost::progress_display pbar(m.n());
  for (unsigned i=0; i<m.n(); ++i, ++pbar) {
    double &d = m.vv.back()[i];
    const Vec< 3,double > a(m.vv[0][i], m.vv[1][i], Ndim>2? m.vv[2][i]:0.);
          Vec< 3,double > b;
    bool found_nonwall_node = false;
    for (unsigned j=0; !wallnodes[i] && j<m.n(); ++j) {
      if (wallnodes[j]) {
        b[0] = m.vv[0][j];
        b[1] = m.vv[1][j];
        b[2] = Ndim>2? m.vv[2][j]:0.;
        const double newd = dist2(a,b);
        d = newd<d? newd:d;
        found_nonwall_node = true;
      }
    }
    for (unsigned j=0; !found_nonwall_node && j<m.n(); ++j) {
      if (i!=j) {
        b[0] = m.vv[0][j];
        b[1] = m.vv[1][j];
        b[2] = Ndim>2? m.vv[2][j]:0.;
        const double newd = dist2(a,b);
        d = newd<d? newd:d;
      }
    }
    d = sqrt(d);
  }
  cout << "info: calculate distances." << endl;


  cout << "::walldistance." << endl;
}

