
#include <sstream>
#include <numeric>
#include "mfactory.h"
#include "f_steven.h"

using namespace std;
using namespace m;


Register< mfoutput,f_steven > mf_steven(3,".steven","Steven output format",
                                          "","--steven-vector-v [str]: variable for velocity",
                                          "","--steven-vector-b [str]: variable for magnetic field");


void f_steven::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  const unsigned dim = m.d();
  if (dim!=2 && dim!=3) {
    cerr << "number of dimensions should be 2 or more!" << endl;
    throw 42;
  }
  const string s_vvec(o.follow("","--steven-vector-v"));
  const string s_bvec(o.follow("","--steven-vector-b"));


  // -- find vector variables indices for velocity and magnetic field --

  int i_vvec = -1;
  int i_bvec = -1;
  for (unsigned i=0; i<m.v(); ++i) {
    i_vvec = (s_vvec.length() && s_vvec==m.vn[i] && i_vvec<0? i:i_vvec);
    i_bvec = (s_bvec.length() && s_bvec==m.vn[i] && i_bvec<0? i:i_bvec);
  }
  if (s_vvec.length() && i_vvec<0) {
    cerr << "variable for velocity \"" << s_vvec << "\" not found!" << endl;
    throw 42;
  }
  if (s_bvec.length() && i_bvec<0) {
    cerr << "variable for magnetic field \"" << s_bvec << "\" not found!" << endl;
    throw 42;
  }


  // -- reorder zones, setting the first one as the volume --

  vector< const mzone* > z;
  for (unsigned i=0; i<m.vz.size(); ++i) {
    if (m.d(i)==dim && (m.vz[i].t==FETRIANGLE || m.vz[i].t==FETETRAHEDRON)) {
      z.push_back(&m.vz[i]);
      break;  // only one
    }
  }
  if (z.size()!=1) {
    cerr << "couldn't find a volume zone!" << endl;
    throw 42;
  }
  for (unsigned i=0; i<m.vz.size(); ++i) {
    if (m.d(i)==dim-1)
      z.push_back(&m.vz[i]);
  }

  vector< unsigned > vnbe;  // nb. boundary elements p. zone
  for (unsigned t=1; t<z.size(); ++t)
    vnbe.push_back(z[t]->e2n.size());


  // -- open files --

  ofstream f1((fn+".boundaries").c_str(),ios_base::trunc);
  ofstream f2((fn+".boundaryelements").c_str(),ios_base::trunc);
  ofstream f3((fn+".elements").c_str(),ios_base::trunc);
  ofstream f4((fn+".flowfield").c_str(),ios_base::trunc);
  ofstream f5((fn+".magneticfield").c_str(),ios_base::trunc);
  ofstream f6((fn+".nodes").c_str(),ios_base::trunc);
  if (!f1 || !f2 || !f3 || !f4 || !f5 || !f6) {
    cerr << "error creating files" << endl;
    throw 42;
  }
  f4.precision(15);  // only these file write doubles
  f5.precision(15);  // ...
  f6.precision(15);  // ...


  // -- f1: .boundaries --

  f1 << "Versie 1.0" << endl
     << "[nBoundaries] = " << vnbe.size() << endl;
  for (unsigned t=0; t<vnbe.size(); ++t)
    f1 << "\t<type> = Insulator" << endl;
  f1 << endl;
  for (unsigned t=0; t<vnbe.size(); ++t) {
    const unsigned begin = accumulate(vnbe.begin(),vnbe.begin()+t,0);
    f1 << "[nBoundaryElements" << t << "] = " << vnbe[t] << endl;
    for (unsigned i=0; i<z[t+1]->e2n.size(); ++i)
      f1 << "\t<boundaryElement> = " << begin+i << endl;
    f1 << endl;
  }


  // -- f2: .boundaryelements --

  f2 << "Versie 1.0" << endl
     << accumulate(vnbe.begin(),vnbe.end(),0) << endl;
  for (unsigned t=1; t<=vnbe.size(); ++t)
    for (unsigned i=0; i<z[t]->e2n.size(); ++i) {
      const vector< unsigned >& en = z[t]->e2n[i].n;
      for (unsigned j=0; j<en.size(); ++j)
        f2 << en[j] << '\t';
      f2 << endl;
    }


  // -- f3: .elements --

  f3 << "Versie 1.0" << endl
     << z[0]->e2n.size() << endl;
  for (unsigned i=0; i<z[0]->e2n.size(); ++i) {
    const vector< unsigned >& en = z[0]->e2n[i].n;
    for (unsigned j=0; j<en.size(); ++j)
      f3 << en[j] << '\t';
    f3 << endl;
  }


  // -- f4: .flowfield --

  f4 << "Versie 1.0" << endl
     << m.n() << endl;
  for (unsigned j=0; j<m.n(); ++j) {
    for (int i=0; i<(int) dim; ++i)
      f4 << (i? "\t":"") << (i_vvec<0? 0.:m.vv[i+i_vvec][j]);
    f4 << endl;
  }

  // -- f5: .magneticfield --

  f5 << "Versie 1.0" << endl
     << m.n() << endl;
  for (unsigned j=0; j<m.n(); ++j) {
    for (int i=0; i<(int) dim; ++i)
      f5 << (i? "\t":"") << (i_bvec<0? 0.:m.vv[i+i_bvec][j]);
    f5 << endl;
  }


  // -- f6: .nodes --

  f6 << "Versie 1.0" << endl
     << m.n() << endl;
  for (unsigned j=0; j<m.n(); ++j) {
    for (unsigned i=0; i<dim; ++i)
      f6 << (i? "\t":"") << m.vv[i][j];
    f6 << endl;
  }


  // -- close files --

  f1.close();
  f2.close();
  f3.close();
  f4.close();
  f5.close();
  f6.close();
}


