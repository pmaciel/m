
#include <sstream>
#include <memory>
#include "boost/progress.hpp"
#include "mfactory.h"
#include "mlinearsystem.h"
#include "t_laplace.h"

using namespace std;
using namespace m;


Register< mtransform,t_laplace > mt_laplace( 19,
  "-tlaplace", "[str] 2D Laplace equation solver, with",
  "", "[str] filename or string with xml formatted as:",
  "", "<laplace",
  "", " ls=\"\" linear system solver (default ls_gauss)",
  "", " Nb=\"\" block size (default 1)",
  "", ">",
  "", " <ls",
  "", "  mtype=\"\"    (default msr)",
  "", "  output=\"\"   (default 1)",
  "", "  precond=\"\"  (default none)",
  "", "  overlap=\"\"  (default 0)",
  "", "  solver=\"\"   (default gmres)",
  "", "  max_iter=\"\" (default 500)",
  "", "  kspace=\"\"   (default 30)",
  "", "  tol=\"\"      (default 1.e-6)",
  "", " />",
  "", " <i zone=\"\" (no default) conductivity=\"\" (default 1.)/> (inner zone, 1 or more)",
  "", " <b zone=\"\" (no default) value=\"\" (default 0)/> (boundary zone, 1 or more)",
  "", "</laplace>" );


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


void t_laplace::transform(GetPot& o, mmesh& m)
{
  if (m.d()!=2 || !m.z())
    return;
  using namespace aux;

  cout << "info: setup laplace xml..." << endl;
  const string o_xml = o.get(o.inc_cursor(),"");
  XMLNode x = ((o_xml.size()? o_xml[0]:'?')=='<')? XMLNode::parseString(o_xml.c_str(),"laplace")
                                                 : XMLNode::openFileHelper(o_xml.c_str(),"laplace");
  cout << "info: setup laplace xml." << endl;


  // set block size
  const unsigned Nb = x.getAttribute< unsigned >("Nb",1);


  cout << "info: minimize mesh connectivity..." << endl;
  {
    vector< bool > vused(m.z(),false);
    for (int i=0; i<x.nChildNode("i"); ++i)
      vused[ getzoneidx( x.getChildNode("i",i).getAttribute< string >("zone") ,m) ] = true;
    for (int i=0; i<x.nChildNode("b"); ++i)
      vused[ getzoneidx( x.getChildNode("b",i).getAttribute< string >("zone") ,m) ] = true;
    for (int i=m.z()-1; i>=0; --i)
      if (!vused[i] || (m.d(i)!=2 && m.d(i)!=1))
        m.vz.erase( m.vz.begin()+i );
    m.compress();
  }
  cout << "info: minimize mesh connectivity." << endl;


  cout << "info: minimize mesh point cloud..." << endl;
  m.vn.resize(2+Nb);
  m.vv.resize(2+Nb);
  for (unsigned i=0; i<Nb; ++i) {
    m.vn[2+i] = string("T") + (char) ('0'+i);
    m.vv[2+i].assign(m.n(),0.);
  }
  cout << "info: minimize mesh point cloud." << endl;


  cout << "info: setup linear system..." << endl;

  // create pointer to linear system and update its options from lap2d xml
  auto_ptr< mlinearsystem< double > > ls(
    Create< mlinearsystem< double > >(x.getAttribute< string >("ls","ls_gauss")) );
  XMLNode x_ls  = x.getChildNode("ls");
  for (int a=0; a<x_ls.nAttribute(); ++a)
    ls->xml.updateAttribute(x_ls.getAttribute(a).lpszValue,NULL,x_ls.getAttribute(a).lpszName);
  ls->initialize(m.n(),m.n(),Nb);
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
  cout << "info: setup linear system." << endl;


  cout << "info: assemble linear system..." << endl;
  {
    // timing and progress utilities
    unsigned Nelem = 0;
    for (vector< mzone >::const_iterator z=m.vz.begin(); z!=m.vz.end(); ++z)
      Nelem += (unsigned) z->e2n.size();
    boost::progress_display pbar(Nelem);
    boost::progress_timer t(cout);

    // assembly
    for (int i=0; i<x.nChildNode("i"); ++i) {
      vector< mzone >::const_iterator z = getzoneit(x.getChildNode("i",i).getAttribute< string >("zone"),m);
      vector< double > conductivity = getvvalues(x.getChildNode("i",i).getAttribute< string >("conductivity","1."));
      for (vector< melem >::const_iterator e=z->e2n.begin(); e!=z->e2n.end(); ++e, ++pbar) {

        const vector< unsigned >& en = e->n;
        const AElement E(en,m.vv[0],m.vv[1]);
        for (unsigned j=0; j<(unsigned) en.size(); ++j)
          for (unsigned k=0; k<(unsigned) en.size(); ++k)
            for (unsigned b=0; b<Nb; ++b)
              (ls->A)(en[j],en[k],b,b) += 0.25 * conductivity[b] / E.Size() *
                 (E.NX(j)*E.NX(k) + E.NY(j)*E.NY(k));

      }
    }
    for (int i=0; i<x.nChildNode("b"); ++i) {
      vector< mzone >::const_iterator z = getzoneit(x.getChildNode("b",i).getAttribute< string >("zone"),m);
      vector< double > value = getvvalues(x.getChildNode("b",i).getAttribute< string >("value","1."));
      for (vector< melem >::const_iterator e=z->e2n.begin(); e!=z->e2n.end(); ++e, ++pbar) {

        for (vector< unsigned >::const_iterator n=e->n.begin(); n!=e->n.end(); ++n)
          for (unsigned b=0; b<Nb; ++b) {
            ls->zerorow(*n,b);
            (ls->A)(*n,*n,b,b) = 1.;
            (ls->B)(*n,b)      = value[b];
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
  for (unsigned i=0; i<Nb; ++i)
    for (unsigned n=0; n<m.n(); ++n)
      m.vv[2+i][n] = (ls->X)(n,i);
}


vector< double > t_laplace::getvvalues(const string& s)
{
  vector< double > r;

  // split string find a '=' then a ':'
  string::size_type p1 = 0,
                    p2 = 0;
  while (p2!=string::npos) {
    double v;
    p2 = s.find(":",p1);
    istringstream ss( s.substr(p1,(p2==string::npos? p2:p2-p1)) );
    ss >> v;
    p1 = p2+1;
    r.push_back(v);
  }

  return r;
}


vector< mzone >::const_iterator t_laplace::getzoneit(const string& n, const mmesh& m)
{
  return m.vz.begin() + getzoneidx(n,m);
}


unsigned t_laplace::getzoneidx(const string& n, const mmesh& m)
{
  for (unsigned r=0; r!=m.z(); ++r)
    if (m.vz[r].n==n)
      return r;
  cerr << "error: zone \"" << n << "\" not present!" << endl;
  throw 42;
  return 0;
}

