
#include <sstream>
#include <memory>
#include "boost/progress.hpp"
#include "mfactory.h"
#include "mlinearsystem.h"
#include "t_lap2d.h"

using namespace std;
using namespace m;


Register< mtransform,t_lap2d > mt_lap2d(6,"-tlap2d","[str[=real]:...] [str] 2d laplace equation solver,",
                                          "",       "with [str[=real]:...]:",
                                          "",       "  zone name[=conductivity] for inner elements, or",
                                          "",       "  zone name=value pairs for boundary elements, and",
                                          "",       "with [str] linear system solver:",
                                          "",       "  ls_gauss|ls_pardiso|ls_trilinos|ls_wsmp");


namespace aux {


// element definition
struct AElement {
  unsigned S;           // number of elem. nodes
  vector< double >  x,  // coord. x copy
                    y;  // ...    y copy
  AElement(const vector< unsigned >& en, const vector< double >& vx, const vector< double >& vy) : S(en.size()), x(S), y(S) {
    for (unsigned i=0; i<en.size(); ++i) {
      x[i] = vx[en[i]];
      y[i] = vy[en[i]];
    }
  }
  double Size() const {
    return (S==3? .5*(x[1]*y[2]-x[2]*y[1] - x[0]*y[2]+x[2]*y[0] + x[0]*y[1]-x[1]*y[0]) :
                  0. );
  }
  double NX(unsigned i) const { return y[(i+1)%S] - y[(i+2)%S]; }
  double NY(unsigned i) const { return x[(i+2)%S] - x[(i+1)%S]; }
};


}


void t_lap2d::transform(GetPot& o, mmesh& m)
{
  if (m.d()!=2)
    return;
  const string o_zones  = o.get(o.inc_cursor(),""),
               o_method = o.get(o.inc_cursor(),"");
  using namespace aux;


  cout << "info: setup zones and mesh..." << endl;
  vector< pair< string,double > > vzones = getzones(o_zones,m);

  // minimize point cloud
  m.vn.resize(3);
  m.vv.resize(3);
  m.vn.back() = "T";
  m.vv.back().assign(m.n(),0.);

  // minimize connectivities
  {
    vector< bool > vused(m.z(),false);
    for (vector< pair< string,double > >::const_iterator p=vzones.begin(); p!=vzones.end(); ++p)
      vused[ getzoneidx(p->first,m) ] = true;
    for (unsigned i=m.z(); i>0; --i)
      if (!vused[i-1])
        m.vz.erase( m.vz.begin()+i-1 );
  }
  m.compress();
  cout << "info: setup zones and mesh." << endl;


  // linear system
  auto_ptr< mlinearsystem< double > > ls(Create< mlinearsystem< double > >(o_method));
  ls->initialize(m.n(),m.n());
  if (ls->issparse) {
    cout << "info: setup linear system sparsity..." << endl;

    vector< vector< unsigned > > nz(m.n());
    for (unsigned i=0; i<m.z(); ++i)
      for (unsigned j=0; j<m.e(i); ++j) {
        const vector< unsigned >& en = m.vz[i].e2n[j].n;
        for (vector< unsigned >::const_iterator n1=en.begin(); n1!=en.end(); ++n1)
          for (vector< unsigned >::const_iterator n2=en.begin(); n2!=en.end(); ++n2)
            nz[*n1].push_back(*n2);
      }
    for (vector< vector< unsigned > >::iterator n=nz.begin(); n!=nz.end(); ++n) {
      sort(n->begin(),n->end());
      n->erase(unique(n->begin(),n->end()),n->end());
    }
    ls->initialize(nz);

    cout << "info: setup linear system sparsity." << endl;
  }


  cout << "info: assemble linear system..." << endl;
  {
    // timing and progress utilities
    unsigned Nelem = 0;
    for (vector< pair< string,double > >::iterator p=vzones.begin(); p!=vzones.end(); ++p) {
      vector< mzone >::const_iterator z = getzoneit(p->first,m);
      Nelem += (unsigned) z->e2n.size();
    }
    boost::progress_display pbar(Nelem);
    boost::progress_timer t(cout);

    for (vector< pair< string,double > >::iterator p=vzones.begin(); p!=vzones.end(); ++p) {
      vector< mzone >::const_iterator z = getzoneit(p->first,m);
      for (vector< melem >::const_iterator e=z->e2n.begin(); e!=z->e2n.end(); ++e, ++pbar) {
        if (z->d()==2) {
          const vector< unsigned >& en = e->n;
          const AElement E(en,m.vv[0],m.vv[1]);
          for (unsigned i=0; i<(unsigned) en.size(); ++i)
            for (unsigned j=0; j<(unsigned) en.size(); ++j)
              (ls->A)(en[i],en[j]) += 0.25 * p->second / E.Size() *
                 (E.NX(i)*E.NX(j) + E.NY(i)*E.NY(j));
        }
        else if (z->d()==1) {
          for (vector< unsigned >::const_iterator n=e->n.begin(); n!=e->n.end(); ++n) {
            ls->zerorow(*n);
            (ls->A)(*n,*n) = 1.;
            (ls->B)(*n)    = p->second;
          }
        }
      }
    }
    cout << "info: timer: ";
  }
  cout << "info: assemble linear system." << endl;


  cout << "info: solve linear system..." << endl;
  {
    boost::progress_timer t(cout);
    ls->solve();
    cout << "info: timer: ";
  }
  cout << "info: solve linear system." << endl;


  // copy solution
  for (unsigned i=0; i<m.n(); ++i)
    m.vv[2][i] = (ls->X)(i);
}


vector< pair< string,double > > t_lap2d::getzones(const string& s, const mmesh& m)
{
  vector< pair< string,double > > r;

  // split string find a '=' then a ':'
  string::size_type p1 = 0,
                    p2 = 0;
  while (p2!=string::npos) {
    pair< string,double > p("",0.);

    p2 = min(s.find(":",p1),s.find("=",p1));
    p.first = s.substr(p1,(p2==string::npos? p2:p2-p1));
    if (s.find("=",p1)<s.find(":",p1)) {
      p1 = p2+1;
      p2 = s.find(":",p1);
      string s2 = s.substr(p1,(p2==string::npos? p2:p2-p1));
      istringstream ss(s2);
      ss >> p.second;
    }

    p1 = p2+1;
    r.push_back(p);
  }

  // check zones exist and comply with dimensions
  for (vector< pair< string,double > >::iterator p=r.begin(); p!=r.end(); ++p) {
    vector< mzone >::const_iterator z = getzoneit(p->first,m);
    if (z->d()==1) {
      cout << "info: zone \"" << z->n << "\": impose = " << p->second << endl;
    }
    else if (z->d()==2) {
      p->second = max(1.e-20,p->second);
      cout << "info: zone \"" << z->n << "\": sigma = " << p->second << endl;
    }
    else {
      cerr << "error: zone \"" << z->n << "\" d!=1 && d!=2" << endl;
      throw 42;
    }
  }

  return r;
}


vector< mzone >::const_iterator t_lap2d::getzoneit(const string& n, const mmesh& m)
{
  return m.vz.begin() + getzoneidx(n,m);
}


unsigned t_lap2d::getzoneidx(const string& n, const mmesh& m)
{
  for (unsigned r=0; r!=m.z(); ++r)
    if (m.vz[r].n==n)
      return r;
  cerr << "error: zone \"" << n << "\" not present!" << endl;
  throw 42;
  return 0;
}

