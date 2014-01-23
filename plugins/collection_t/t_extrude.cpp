
#include "mfactory.h"
#include "t_extrude.h"
#include <sstream>

using namespace m;


Register< mtransform,t_extrude > mt_extrude(3,"-tz","[z1:n:z2] mesh extrusion from constant steps, or",
                                              "",   "[z1:r1:n:r2:z2]  ... from geometric progression",
                                              "","--extrude-wires: connect extrusion layers with \"wires\" (default: no)");


void t_extrude::transform(GetPot& o, mmesh& m2, const XMLNode& x)
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::vector;
  using std::string;

  cout << "info: [d/n/e]: " << m2.d() << '/' << m2.n() << '/' << m2.e() << "..." << endl;
  const unsigned D2 = m2.d(),
                 N2 = m2.n();
  if (D2>2)
    return;
  const string str = o.get(o.inc_cursor(),"");
  const bool wires = o.search("--extrude-wires");


  cout << "info: z-steps..." << endl;
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


  // mesh structure, extruded
  mmesh m3;


  // -- variables --


  cout << "info: variable names and values..." << endl;
  {
    // resize (new) variables names and values
    const vector< bool > vectors2d = m2.vvectors();
    const size_t vextra = count(vectors2d.begin(),vectors2d.end(),true);
    m3.vn.assign(m2.vn.size()+vextra,"");
    m3.vv.assign(m2.vn.size()+vextra,vector< double >(Nzsteps*N2,0.));

    // variable names and values iterators
    vector< string           >::iterator n3 = m3.vn.begin();
    vector< vector< double > >::iterator v3 = m3.vv.begin();
    for (unsigned i=0; i<vectors2d.size(); ++i) {
      vector< string           >::const_iterator n2 = m2.vn.begin() + i;
      vector< vector< double > >::const_iterator v2 = m2.vv.begin() + i;
      if (vectors2d[i]) {

        // add vector field names
        const char   x = *(n2->end()-1);
        const string n = n2->substr(0,n2->length()-1);
        for (unsigned j=0; j<=D2; ++j, ++n3)
          *n3 = n + char(x+j);

        // add vector field values
        for (unsigned j=0; j<D2; ++j, ++v2, ++v3)
          for (unsigned l=0; l<Nzsteps; ++l)
            copy(v2->begin(),v2->end(),(v3->begin()) + (l*N2));
        ++v3;

        ++i;  // one variable was added
      }
      else {

        // add scalar field name/values
        *n3 = *n2;
        for (unsigned l=0; l<Nzsteps; ++l)
          copy(v2->begin(),v2->end(),(v3->begin()) + (l*N2));
        ++n3;
        ++v3;

      }
    }

    // set extrusion layer values
    v3 = m3.vv.begin() + D2;  // extruded coordinate
    for (unsigned l=0; l<Nzsteps; ++l)
      fill_n((v3->begin()) + (l*N2),N2,m_zsteps[l]);
  }
  cout << "info: variable names and values." << endl;


  // -- zones --


  cout << "info: allocate new zones..." << endl;
  vector< mzone* > vzbottomtop;
  for (unsigned i=0; i<m2.vz.size(); ++i) {
    if (m2.d(i)==D2) {
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
    const unsigned Nshift = (Nzsteps-1)*N2;
    for (unsigned c=0; c<zt.e2n.size(); ++c) {
      vector< unsigned >& en = zt.e2n[c].n;
      for (unsigned j=0; j<en.size(); ++j)
        en[j] += Nshift;
      switch (zt.t) {
        case FELINESEG:
        case FETRIANGLE:      { std::swap(en[0],en[1]); } break;
        case FEQUADRILATERAL: { std::swap(en[0],en[2]); } break;
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
          ei.n[n] += N2;
        ef.n.insert(ef.n.end(),ei.n.begin(),ei.n.end());
        switch (m2.vz[t].t) {
          case FEQUADRILATERAL: { std::swap(ef.n[6],ef.n[7]); }
          case FELINESEG:       { std::swap(ef.n[2],ef.n[3]); }
          default: break;
        }
        m3.vz[t].e2n.push_back(ef);

      }
    }
  }
  cout << "info: connect extrusion layers." << endl;


  if (wires) {
    cout << "info: connect wires..." << endl;

    // deafult element type for wires is a line segment
    melem default_e;
    default_e.n.assign(2,0);

    // create new wire zones (separate)
    vector< mzone > w(N2+2*N2+1);
    for (unsigned i=0; i<N2; ++i) {

      // set name and type
      std::ostringstream s;
      s << "wire_" << i+1;
      w[N2*0+i].n = w[N2*1+i].n = w[N2*2+i].n = s.str();
      w[N2*1+i].n += "_bottom";
      w[N2*2+i].n += "_top";
      w[N2*0+i].t = w[N2*1+i].t = w[N2*2+i].t = FELINESEG;

      // allocate wire connectivity
      w[N2*0+i].e2n.assign(Nzsteps-1,default_e);
      w[N2*1+i].e2n.assign(1,default_e);
      w[N2*2+i].e2n.assign(1,default_e);

      // set wire connectivity
      for (unsigned j=0; j<Nzsteps-1; ++j) {
        w[N2*0+i].e2n[j].n[0] = i+N2*(j+0);
        w[N2*0+i].e2n[j].n[1] = i+N2*(j+1);
      }
      w[N2*1+i].e2n[0].n[0] = w[N2*1+i].e2n[0].n[1] = i;
      w[N2*2+i].e2n[0].n[0] = w[N2*2+i].e2n[0].n[1] = i+N2*(Nzsteps-1);
    }

    // create a zone with all the wires together
    mzone& allw = w.back();
    allw.n = "wires";
    allw.t = FELINESEG;
    allw.e2n.assign(N2*(Nzsteps-1),default_e);
    for (unsigned i=0; i<N2; ++i) {
      for (unsigned j=0; j<Nzsteps-1; ++j) {
        allw.e2n[i*(Nzsteps-1)+j].n[0] = i+N2*(j+0);
        allw.e2n[i*(Nzsteps-1)+j].n[1] = i+N2*(j+1);
      }
    }

    // append to mesh zones
    m3.vz.insert(m3.vz.end(),w.begin(),w.end());

    cout << "info: connect wires." << endl;
  }


  m2 = m3;
  cout << "info: [d/n/e]: " << m2.d() << '/' << m2.n() << '/' << m2.e() << "." << endl;
}


void t_extrude::zsteps_idf(const std::string& str)
{
  using namespace std;

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


void t_extrude::zsteps_geomp(const std::string& str)
{
  using namespace std;

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

