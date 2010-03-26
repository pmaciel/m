
#include "mfactory.h"
#include "t_extrude.h"
#include <sstream>

using namespace std;
using namespace m;


Register< mtransform,t_extrude > mt_extrude(2,"-tz","[z1:n:z2] mesh extrusion from constant steps, or",
                                              "",   "[z1:r1:n:r2:z2]  ... from geometric progression");


void t_extrude::transform(GetPot& o, mmesh& m2)
{
  const unsigned dim = m2.d();
  if (dim>2)
    return;
  cout << "::extrusion [n/e]: " << m2.n() << " / " << m2.e() << "..." << endl;


  // -- point cloud --


  cout << "info: z-steps..." << endl;
  const string str = o.get(o.inc_cursor(),"");
  const unsigned ncolon = count(str.begin(),str.end(),':');
  m_zsteps.clear();
  if      (ncolon==2)  zsteps_idf(str);
  else if (ncolon==4)  zsteps_geomp(str);
  else {
    cerr << "error: incorrect z-steps definition" << endl;
    throw 42;
  }
  const unsigned Nzsteps = m_zsteps.size();
  if (Nzsteps<=1) {
    cerr << "error: not enough z-steps" << endl;
    throw 42;
  }
  sort(m_zsteps.begin(),m_zsteps.end());
  cout << "info: z-steps." << endl;


  cout << "info: point cloud variable names..." << endl;
  mmesh m3;
  const vector< bool > vectors2d = m2.vvectors();
  m3.vn.clear();
  for (unsigned i=0; i<vectors2d.size(); ++i) {
    if (vectors2d[i]) {
      const char x = *(m2.vn[i].end()-1);
      const string n = m2.vn[i].substr(0,m2.vn[i].length()-1);
      m3.vn.push_back(n + char(x));
      m3.vn.push_back(n + char(x+1));
      m3.vn.push_back(n + char(x+2));
      i+=1;
    }
    else
      m3.vn.push_back(m2.vn[i]);
  }
  cout << "info: point cloud variable names." << endl;


  cout << "info: point cloud variable values copy..." << endl;
  const vector< bool > vectors3d = m3.vvectors();
  m3.vv.resize(m3.vn.size());
  for (unsigned j=0, i=0; j<vectors3d.size(); ++j, ++i) {
    m3.vv[j].reserve(Nzsteps*m2.n());
    if (vectors3d[j]) {

      // add vector field
      for (unsigned d=0; d<2; ++d)
        for (unsigned l=0; l<Nzsteps; ++l)
          m3.vv[j+d].insert(m3.vv[j+d].end(),m2.vv[i+d].begin(),m2.vv[i+d].end());
      m3.vv[j+2].assign(Nzsteps*m2.n(),0.);
      i+=1;  // skip y-component
      j+=2;  // skip y and z-component

    }
    else {

      // add scalar field
      for (unsigned l=0; l<Nzsteps; ++l)
        m3.vv[j].insert(m3.vv[j].end(),m2.vv[i].begin(),m2.vv[i].end());

    }
  }
  cout << "info: point cloud variable values copy." << endl;


  cout << "info: point cloud extrusion layer values..." << endl;
  m3.vv[dim].clear();
  for (unsigned l=0; l<Nzsteps; ++l)
    m3.vv[dim].insert(m3.vv[dim].end(),m2.n(),m_zsteps[l]);
  cout << "info: point cloud extrusion layer values." << endl;


  // -- zones --


  cout << "info: allocate new zones..." << endl;
  vector< mzone* > vzbottomtop;
  for (unsigned i=0; i<m2.vz.size(); ++i) {
    if (m2.d(i)==dim) {
      vzbottomtop.push_back(&(m2.vz[i]));
      cout << "info: zone to replicate: \"" << m2.vz[i].n << "\"" << endl;
    }
  }

  m3.vz.resize(m2.vz.size()+vzbottomtop.size()*2);
  for (vector< mzone >::iterator z2=m2.vz.begin(), z3=m3.vz.begin(); z2!=m2.vz.end(); ++z2, ++z3) {
    z3->n = z2->n;
    z3->t = (z2->t==ORDERED?         FELINESEG       :
            (z2->t==FELINESEG?       FEQUADRILATERAL :
            (z2->t==FETRIANGLE?      PRISM3          :
            (z2->t==FEQUADRILATERAL? FEBRICK         :
                                     ORDERED ))));
    z3->e2n.reserve((z2->e2n.size())*(Nzsteps-1));
  }
  cout << "info: allocate new zones." << endl;


  cout << "info: top/bottom replication..." << endl;
  // (it's here because extrusion is destructive)
  for (unsigned i=0; i<vzbottomtop.size(); ++i) {
    mzone& zb = m3.vz[ m2.vz.size() + i*2+0 ];
    mzone& zt = m3.vz[ m2.vz.size() + i*2+1 ];
    zb = *vzbottomtop[i];  zb.n.append("_extruded_bottom");
    zt = *vzbottomtop[i];  zt.n.append("_extruded_top");

    // shift top zone element indices and flip orientation
    const unsigned Nshift = (Nzsteps-1)*m2.n();
    for (unsigned c=0; c<zt.e2n.size(); ++c) {
      vector< unsigned >& en = zt.e2n[c].n;
      for (unsigned j=0; j<en.size(); ++j)
        en[j] += Nshift;
      switch (zt.t) {
        case FELINESEG:
        case FETRIANGLE:      { swap(en[0],en[1]); } break;
        case FEQUADRILATERAL: { swap(en[0],en[2]); } break;
        case FEPOLYGON:  // not extrudible (yet)
        case ORDERED:
        default: break;
      }
    }
  }
  cout << "info: top/bottom replication." << endl;


  cout << "info: connect extrusion layers..." << endl;
  // (this destroys m2 zones connectivities)
  for (unsigned l=0; l<Nzsteps; ++l) {
    cout << "info: extrusion layer z=" << m_zsteps[l] << endl;
    for (unsigned t=0; t<m2.vz.size(); ++t) {
      for (unsigned i=0; l>0 && i<m2.e(t); ++i) {
        melem& ei = m2.vz[t].e2n[i];
        melem  ef = ei;
        for (unsigned n=0; n<ei.n.size(); ++n)
          ei.n[n] += m2.n();
        ef.n.insert(ef.n.end(),ei.n.begin(),ei.n.end());
        switch (m2.vz[t].t) {
          case FEQUADRILATERAL: { swap(ef.n[6],ef.n[7]); }
          case FELINESEG:       { swap(ef.n[2],ef.n[3]); }
          default: break;
        }
        m3.vz[t].e2n.push_back(ef);

      }
    }
  }
  cout << "info: connect extrusion layers." << endl;


  m2 = m3;
  cout << "::extrusion [n/e]: " << m2.n() << " / " << m2.e() << "." << endl;
}


void t_extrude::zsteps_idf(const string& str)
{
  double z1=0., z2=1.;
  unsigned n=2;
  string s = str;
  replace(s.begin(),s.end(),':',' ');
  istringstream is(s);
  is >> z1 >> n >> z2;
  if (n<2) {
    cerr << "extrusion error: z-steps n should be greater than 1" << endl;
    throw 42;
  }

  const double d = (z2-z1)/((double) n-1);
  m_zsteps.assign(1,z1);
  while (m_zsteps.size()<n)
    m_zsteps.push_back(m_zsteps.back()+d);
}


void t_extrude::zsteps_geomp(const string& str)
{
  double z1=0., r1=1., r2=1., z2=1.;
  unsigned n=2;
  string s = str;
  replace(s.begin(),s.end(),':',' ');
  istringstream is(s);
  is >> z1 >> r1 >> n >> r2 >> z2;
  if (n<2) {
    cerr << "extrusion error: z-steps n should be greater than 1" << endl;
    throw 42;
  }

  // create a double-sided progression from 0. to 1.
  vector< double > g;
  {
    // number of intervals is number of layers minus one
    const unsigned ni = n-1;

    // first do two geom. progressions starting with 1.; if ni is odd,
    // last element of g1 is averaged with last value of g2
    vector< double > g1(1,1.),
                     g2(1,1.);
    for (unsigned i=1; g1.size()+g2.size()<ni; ++i) {
      g1.push_back(r1*g1[i-1]);
      g2.push_back(r2*g2[i-1]);
    }
    if (ni%2) {
     g1.back() = (g1.back()+g2.back())*.5;
     g2.pop_back();
    }

    // introduce the reference entry (thus 1 + ni intervals = n layers),
    // append the two progressions and accumulate
    g.assign(1,0.);
    for (vector< double >::iterator i=g1.begin(); i!=g1.end(); ++i)
      g.push_back(*i);
    for (vector< double >::reverse_iterator i=g2.rbegin(); i!=g2.rend(); ++i)
      g.push_back(*i);

    // accumulate and bound up to 1.
    for (unsigned i=1; i<g.size(); ++i)
      g[i] += g[i-1];
    const double gmin = 0.;  // by definition
    const double gmax = *max_element(g.begin(),g.end());
    for (unsigned i=0; i<g.size(); ++i)
      g[i] = (g[i]-gmin)/(gmax-gmin);
  }

  // set zsteps, bounded between z1:z2
  m_zsteps.resize(g.size());
  for (unsigned i=0; i<g.size(); ++i)
    m_zsteps[i] = g[i]*(z2-z1)+z1;
}

