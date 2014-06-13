
#include <sstream>

#include "mfactory.h"

#include "utils.h"
#include "t_zbool.h"


using namespace m;


Register< mtransform,t_zbool > mt_zbool( 3,
  "-tzunion", "[str:...] zone union (only for zones of the same type)",
  "-tzinter", "[str:...] zone intersection",
  "-tzsplit-sharednodes", "[str:str:...] ..., based on other zone(s) sharing edges" );


void t_zbool::transform(GetPot& o, mmesh& m)
{
  using namespace std;
  cout << "::zone boolean operation..." << endl;


  // get operation, zone names and indices
  const string k = o[o.get_cursor()];
  const vector< string > znames(utils::split(o.get(o.inc_cursor(),""),':'));
  vector< unsigned > zidx;
  for (vector< string >::const_iterator i=znames.begin(); i!=znames.end(); ++i) {
    try {
      const unsigned idx = utils::getzindex(m,*i);
      zidx.push_back(idx);
    }
    catch (const int&) { continue; }
  }

  if (zidx.size()>1) {
    if      (k=="-tzunion") { zunion(m,zidx); }
    else if (k=="-tzinter") { zinter(m,zidx); }
    else if (k=="-tzsplit-sharednodes") { zsplit_sharednodes(m,zidx); }
  }
  else {
    cerr << "warn: nothing to do.";
  }

  cout << "::zone boolean operation." << endl;
}


void t_zbool::zunion(mmesh& m, const std::vector< unsigned >& zidx)
{
  using namespace std;

  // assert same type and reserve new connectivity size
  mzone& z = m.vz[ zidx[0] ];  // (this zone merges the others)
  unsigned
      Nelemi = z.e2n.size(),  // original nb. elements
      Nelemf = Nelemi;        // final nb. elements
  for (vector< unsigned >::const_iterator i=zidx.begin()+1; i!=zidx.end(); ++i) {
    if (z.t!=m.vz[*i].t) {
      cerr << "error: zones must be of same type" << endl;
      throw 42;
    }
    else{
      z.n.append("_union_" + m.vz[*i].n);
      Nelemf += m.vz[*i].e2n.size();
    }
  }
  z.e2n.resize(Nelemf);

  // do an hostile take-over
  for (unsigned i=1, k=Nelemi; i<zidx.size(); ++i) {
    for (unsigned j=0; j<m.vz[ zidx[i] ].e2n.size(); ++j) {
      vector< unsigned >
       &en1 = z.e2n[k++].n,
       &en2 = m.vz[ zidx[i] ].e2n[j].n;
      en1.swap(en2);
    }
  }

  // delete useless zones
  for (unsigned i=zidx.size()-1; i>=1; --i)
    m.vz.erase( m.vz.begin()+zidx[i] );

  cout << "info: zone \"" << z.n << "\" [e]: " << Nelemi << '>' << Nelemf << endl;
}


void t_zbool::zinter(mmesh& m, const std::vector< unsigned >& zidx)
{
  using namespace std;

  // create list of usable nodes for the intersected zone
  // (also, generate a name for the new zone with the intersected elements)
  vector< bool > oknodes(m.n(),false);
  string n(m.vz[ zidx[0] ].n);
  for (vector< unsigned >::const_iterator i=zidx.begin()+1; i!=zidx.end(); ++i) {
    n += "_inter_" + m.vz[*i].n;
    for (vector< melem >::const_iterator j=m.vz[*i].e2n.begin(); j!=m.vz[*i].e2n.end(); ++j)
      for (vector< unsigned >::const_iterator k=j->n.begin(); k!=j->n.end(); ++k)
        oknodes[*k] = true;
  }

  // create new zone, intersecting the first with a union of the rest
  mzone
   &zi(m.vz[ zidx[0] ]),  // reference
    zf(n,zi.t);           // (copy!)
  for (vector< melem >::const_iterator j=zi.e2n.begin(); j!=zi.e2n.end(); ++j) {
    bool okelem = true;
    for (vector< unsigned >::const_iterator k=j->n.begin(); okelem && k!=j->n.end(); ++k)
      okelem = oknodes[*k];
    if (okelem)
      zf.e2n.push_back(*j);
  }

  // append to mesh if elements are found, otherwise issue warning
  if (zf.e2n.size()) {
    m.vz.push_back(zf);
    cout << "info: zone \"" << zf.n << "\" [e]: " << m.vz[ zidx[0] ].e2n.size() << '>' << m.vz.back().e2n.size() << endl;
  }
  else {
    cerr << "warn: zone intersection generated no valid elements." << endl;
  }
}


void t_zbool::zsplit_sharednodes(mmesh &m, const std::vector< unsigned > &zidx)
{
  using namespace std;
  typedef vector< unsigned >::const_iterator vui;

  // build list of used nodes by the zones other than the first
  vector< bool > shared_nodes(m.n(),false);
  for (vui i=zidx.begin()+1; i!=zidx.end() && m.vz[zidx[0]].e2n.size(); ++i)
    for (vector< melem >::const_iterator j=m.vz[*i].e2n.begin(); j!=m.vz[*i].e2n.end(); ++j)
      for (vui k=j->n.begin(); k!=j->n.end(); ++k)
        shared_nodes[*k] = true;

  // append empty zone to contain the splitted elements, and create & move list
  // of all elements to be tested
  string n(m.vz[zidx[0]].n);
  while (utils::zexists(m,n))
    n += "_split";
  m.vz.push_back(mzone(n,m.vz[zidx[0]].t));
  mzone
   &zi(m.vz[zidx[0]]),  // reference!
   &zf(m.vz.back());    // reference (possibly deleted if empty)
  vector< melem > e2n;
  e2n.swap(zi.e2n);

  // check if element in original zone shares any node with the other zones,
  // moving element to appropriate zone according to this
  for (vector< melem >::iterator i=e2n.begin(); i!=e2n.end(); ++i) {
    bool move = false;
    for (vui j=i->n.begin(); j!=i->n.end() && !(move=shared_nodes[*j]); ++j) {}
    (move? zf:zi).e2n.push_back(melem());
    (move? zf:zi).e2n.back().n.swap(i->n);
  }

  // check if splitted zone contains any element, otherwise erase it
  if (zf.e2n.size()) {
    cout << "info: zone \"" << zi.n << "\" [e]: " << e2n.size() << '=' << zi.e2n.size() << '+' << zf.e2n.size() << endl;
  }
  else {
    cerr << "warn: zone splitting generated no valid elements." << endl;
    m.vz.pop_back();
  }
}
