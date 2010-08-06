
#include <sstream>
#include "mfactory.h"
#include "t_ren.h"
#include "t_geo.h"

#include "boost/graph/bandwidth.hpp"
#include "boost/graph/connected_components.hpp"

#include "boost/graph/cuthill_mckee_ordering.hpp"
#include "boost/graph/king_ordering.hpp"


using namespace std;
using namespace m;


//reverse Cuthill-McKee algorithm
Register< mtransform,t_ren > mt_ren(10,"-tren","[str] [int] renumber nodes using boost graph library,",
                                       "",     "with [str] method:",
                                       "",     "  \"rcm\"/\"cm\", for reversed/normal Cuthill-McKee, or",
                                       "",     "  \"rking\"/\"king\", for reversed/normal King, or",
                                       "",     "  \"dump\", for just dump sparsity (\"dump.sparsity.plt\"), and",
                                       "",     "with [int] node:",
                                       "",     "  [i] for start node i, or",
                                       "",     "  [-i] for automatic start node, or",
                                       "",     "  \"b\" for minimum bandwidth nodes, or",
                                       "",     "  \"d\" for minimum degree node");


// auxiliary definitions
namespace aux {


// renumbering method skeleton
struct renumber {
  virtual ~renumber() {}
  virtual void calc(
    // outputs: renumbering maps new-to-old and old-to-new,
    // inputs: start node (negative for automatic) and mesh graph
    vector< unsigned >& B2A, vector< unsigned >& A2B,
    const int start, Graph& G ) = 0;
  virtual void adjustperm(vector< unsigned >& perm) {}
};


// boost graph library shortcuts
using boost::get;
using boost::vertex_color;
using boost::vertex_degree;
using boost::vertex_index;
using boost::vertex_index_t;


// renumbering method: Cuthill-McKee
struct ren_cm : renumber {
  void calc( vector< unsigned >& B2A, vector< unsigned >& A2B,
    const int start, Graph& G )
  {
    // calculate permutation (B2A)
    for (unsigned i=0; i<B2A.size(); ++i)
      B2A[i] = i;
    if (start<0) {
      boost::cuthill_mckee_ordering( G,
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G) );
    }
    else {
      boost::cuthill_mckee_ordering( G,boost::vertex((unsigned) start,G),
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G) );
    }
    adjustperm(B2A);

    // calculate inverse permutation (A2B)
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;
  }
};


// renumbering method: King
struct ren_king : renumber {
  void calc( vector< unsigned >& B2A, vector< unsigned >& A2B,
    const int start, Graph& G )
  {
    // calculate permutation (B2A)
    for (unsigned i=0; i<B2A.size(); ++i)
      B2A[i] = i;
    boost::property_map< Graph,vertex_index_t >::type index_map =
      get(vertex_index,G);
    if (start<0) {
      boost::king_ordering( G,
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }
    else {
      boost::king_ordering( G,vertex((unsigned) start,G),
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }
    adjustperm(B2A);

    // calculate inverse permutation (A2B)
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;
  }
};


// renumbering methods (reversed): Cuthill-McKee, King
struct ren_rcm    : ren_cm    { void adjustperm(vector< unsigned >& perm) { std::reverse(perm.begin(),perm.end()); } };
struct ren_rking  : ren_king  { void adjustperm(vector< unsigned >& perm) { std::reverse(perm.begin(),perm.end()); } };


}  // namespace aux
// auxiliary definitions


void t_ren::transform(GetPot& o, mmesh& m)
{
#if 1
  const unsigned Nnodes = m.n();
  Graph G(Nnodes);

  cout << "info: building graph..." << endl;
  {
    vector< vector< unsigned > > M(Nnodes);

    // use all available zones (hopefully connecting all nodes in mesh)
    for (unsigned j=0; j<m.z(); ++j) {
      const vector< melem >& e2n = m.vz[j].e2n;

      // get all nodes neighbors
      for (unsigned j=0; j<e2n.size(); ++j) {
        const vector< unsigned >& nodes = e2n[j].n;
        const unsigned Nnodes = nodes.size();
        for (unsigned k=0; k<Nnodes; ++k)
          for (unsigned l=1; l<Nnodes; ++l) {  // guarantee "except itself"
            const unsigned i = nodes[k];
            const unsigned j = nodes[(k+l)%Nnodes];
            M[i].push_back(j);
            M[j].push_back(i);
          }
      }
    }

    //  remove duplicate neighbor entries and add graph edges
    unsigned i=0;
    for (vector< vector< unsigned > >::iterator n=M.begin(); n!=M.end(); ++n, ++i) {
      sort(n->begin(),n->end());
      n->erase(unique(n->begin(),n->end()),n->end());
      for (vector< unsigned >::const_iterator j=n->begin(); j<n->end(); ++j)
        if (*j>i)  // guarantee no parallel edges are introduced
          boost::add_edge(i,*j,G);
    }
  }
  if (!boost::num_edges(G)) {
    cerr << "error: couldn't build graph!" << endl;
    throw 42;
  }
  cout << "info: building graph." << endl;
#else

#if 1 // test graph 1: small graph
  enum NODES {N1,N2,N3,N4,N5,N6,N7,N8,Nnodes};
  Graph G(Nnodes);
  boost::add_edge(N1,N5,G);
  boost::add_edge(N5,N3,G);
  boost::add_edge(N3,N2,G);
  boost::add_edge(N2,N6,G);
  boost::add_edge(N2,N8,G);
  boost::add_edge(N6,N8,G);
  boost::add_edge(N4,N7,G);
  boost::add_edge(N4,N5,G); // connects its two components
#else // test graph 2: small graph, three disconnected components
  enum NODES {N1,N2,N3,N4,N5,N6,Nnodes};
  Graph G(Nnodes);
  boost::add_edge(N1,N2,G);
  boost::add_edge(N2,N5,G);
  boost::add_edge(N5,N1,G);
  boost::add_edge(N3,N6,G);
#endif

#endif


  vector< unsigned > new2old(Nnodes,0);
  vector< unsigned > old2new(Nnodes,0);
  for (unsigned i=0; i<Nnodes; ++i)
    new2old[i] = old2new[i] = i;
  const unsigned orig_b = boost::bandwidth(G);
  cout << "info: original bandwidth: " << orig_b << endl;

  ofstream f("dump.sparsity.plt",ios_base::trunc);
  f << "VARIABLES = i j" << endl;
  dump(f,"original",G,old2new);

  const string o_method = o.get(o.inc_cursor(),"");
  const string o_node   = o.get(o.inc_cursor(),"");
  if (o_method=="dump")
    return;


  cout << "info: set method..." << endl;
  aux::renumber* renumber = NULL;
  {
    mfactory< aux::renumber >* f = mfactory< aux::renumber >::instance();
    f->Register(new ProductionLine< aux::renumber,aux::ren_cm     >("cm",    "Cuthill-McKee"));
    f->Register(new ProductionLine< aux::renumber,aux::ren_king   >("king",  "King"));
    f->Register(new ProductionLine< aux::renumber,aux::ren_rcm    >("rcm",   "reversed Cuthill-McKee"));
    f->Register(new ProductionLine< aux::renumber,aux::ren_rking  >("rking", "reversed King"));
    renumber = f->Create(o_method);
  }
  if (renumber==NULL) {
    cerr << "error: incorrect method" << endl;
    throw 42;
  }
  cout << "info: set method." << endl;


  cout << "info: check graph components..." << endl;
  int Ncomp = -1;
  {
    vector< unsigned > vcomp(Nnodes);
    Ncomp = boost::connected_components(G,&vcomp[0]);
    cout << "info: number of graph components: " << Ncomp << endl;
  }
  cout << "info: check graph components." << endl;


  int start = (int) Nnodes;
  if (Ncomp>1) {
    cout << "info: multiple components, automatic start node forced" << endl;
    start = -1;
  }
  else if (o_node=="b") {
    cout << "info: finding minimum bandwidth..." << endl;
    unsigned min_b = orig_b;
    string info;
    for (int s=0; s<(int) Nnodes; ++s) {

      // calculate renumber mappings and bandwidth
      renumber->calc(new2old,old2new,s,G);
      boost::property_map< Graph,boost::vertex_index_t >::type index_map =
        boost::get(boost::vertex_index,G);
      const unsigned calc_b = boost::bandwidth( G,
        boost::make_iterator_property_map(&old2new[0],index_map,old2new[0]) );
      cout << "info: renumber s:" << s << " bw:" << calc_b << (calc_b<min_b? "*":"") << endl;

      // save minimum-bandwidth start node
      if (calc_b<min_b) {
        min_b = calc_b;
        start = s;
      }
    }
    cout << "info: found s:" << start << " bw: " << min_b << endl;
    cout << "info: finding minimum bandwidth." << endl;
  }
  else if (o_node=="d") {
    cout << "info: finding minimum degree..." << endl;
    unsigned min_d = Nnodes;
    for (unsigned i=0; i<Nnodes; ++i) {
      const unsigned deg = boost::degree(boost::vertex(i,G),G);
      if (deg<min_d) {
        min_d = deg;
        start = (int) i;
      }
    }
    cout << "info: node/degree: " << start << '/' << min_d << endl;
    cout << "info: finding minimum degree." << endl;
  }
  else {
    // use given index instead
    istringstream is(o_node);
    is >> start;
  }


  if (start>=(int) Nnodes) {
    cout << "warn: can't find start node, not renumbering." << endl;
    return;
  }


  cout << "info: renumbering..."  << endl;
  {
    // (re)calculate renumber mappings and bandwidth
    renumber->calc(new2old,old2new,start,G);
    boost::property_map< Graph,boost::vertex_index_t >::type index_map =
      boost::get(boost::vertex_index,G);
    const unsigned calc_b = boost::bandwidth( G,
      boost::make_iterator_property_map(&old2new[0],index_map,old2new[0]) );

    // apply mapping if improvement is possible
    ostringstream os;
    os << "renumbering  s:" << start << "  bw:" << orig_b << ">" << calc_b;
    cout << "info: " << os.str() << endl;
    if (calc_b<orig_b) {
      apply(old2new,m);
      dump(f,os.str(),G,old2new);
    }
    else {
      cout << "warn: bandwidth not improved, not renumbering." << endl;
    }
  }
  cout << "info: renumbering."  << endl;

}


void t_ren::dump(ostream& o, const string& t, const Graph& G, const vector< unsigned >& A2B)
{
  o << "ZONE T=\"" << t << "\" I="
    << boost::num_edges(G)*2 + boost::num_vertices(G)
    << endl;

  // upper/lower diagonal then diagonal entries
  boost::graph_traits< Graph >::edge_iterator e, ef;
  for (tie(e,ef) = boost::edges(G); e!=ef; ++e) {
    const int u = (int) A2B[ boost::source(*e,G) ];
    const int v = (int) A2B[ boost::target(*e,G) ];
    o << u << ' ' << -v << endl
      << v << ' ' << -u << endl;
  }
  boost::graph_traits< Graph >::vertex_iterator v, vf;
  for (tie(v,vf) = boost::vertices(G); v!=vf; ++v) {
    const int u = (int) A2B[ *v ];
    o << u << ' ' << -u << endl;
  }
}


void t_ren::apply(const vector< unsigned >& A2B, mmesh& m)
{
  const unsigned Ndim   = m.d();
  const unsigned Nvars  = m.v();
  const unsigned Nnodes = m.n();
  if (!Ndim || !Nvars || !Nnodes)
    return;

  // renumber node list (swapping)
  vector< vector< double > > vv(Nvars,vector< double >(Nnodes,0.));
  for (unsigned v=0; v<Nvars; ++v)
    for (unsigned i=0; i<Nnodes; ++i)
      vv[v][A2B[i]] = m.vv[v][i];
  m.vv.swap(vv);

  // renumber element lists (in-place)
  unsigned Nelem   = 0;
  unsigned Nerrors = 0;
  unsigned Ncorr = 0;
  unsigned Nswap = 0;
  for (unsigned i=0; i<m.z(); ++i) {
    const unsigned Nenodes = (m.e(i)? m.vz[i].e2n[0].n.size():0);
    vector< vector< double > > c(Nenodes,vector< double >(Ndim));
    for (unsigned j=0; j<m.e(i); ++j) {
      vector< unsigned >& enodes = m.vz[i].e2n[j].n;

      // renumber
      for (vector< unsigned >::iterator n=enodes.begin(); n!=enodes.end(); ++n)
        *n = A2B[*n];

      // check and correct
      ++Nelem;
      for (unsigned i=0; i<Nenodes; ++i)
        for (unsigned j=0; j<Ndim; ++j)
          c[i][j] = m.vv[j][ enodes[i] ];
      if (t_geo::cellgeom_ok(c))
        continue;
      ++Nerrors;
      for (unsigned l=0; l<Nenodes; ++l) {
        ++Nswap;
        swap( enodes[l], enodes[(l+1)%Nenodes] );
        for (unsigned i=0; i<Nenodes; ++i)
          for (unsigned j=0; j<Ndim; ++j)
            c[i][j] = m.vv[j][ enodes[i] ];
        if (t_geo::cellgeom_ok(c)) {
          ++Ncorr;
          break;
        }
      }

    }
  }
  if (Nerrors || Nswap || Ncorr) {
    cout << "info: elements/errors: " << Nelem << '/' << Nerrors << endl;
    cout << "info: swaps/corrections: " << Nswap << '/' << Ncorr << endl;
  }
}

