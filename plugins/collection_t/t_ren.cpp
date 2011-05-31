
#include <sstream>
#include "mfactory.h"
#include "t_ren.h"
#include "t_geo.h"

#include "boost/progress.hpp"
#include "boost/graph/bandwidth.hpp"
#include "boost/graph/connected_components.hpp"

#include "boost/graph/cuthill_mckee_ordering.hpp"
#include "boost/graph/king_ordering.hpp"


using namespace std;
using namespace m;


Register< mtransform,t_ren > mt_ren(11,"-tren","[str] [int] renumber nodes using boost graph library,",
                                       "",     "with [str] method:",
                                       "",     "  \"rcm\"/\"cm\", for reversed/normal Cuthill-McKee, or",
                                       "",     "  \"rking\"/\"king\", for reversed/normal King, or",
                                       "",     "  \"dump\", for just dump sparsity (\"dump.sparsity.plt\"), and",
                                       "",     "with [int] node:",
                                       "",     "  [i] for start node i, or",
                                       "",     "  [-i] for automatic start node, or",
                                       "",     "  \"b\" for minimum bandwidth nodes, or",
                                       "",     "  \"d\" for minimum degree node, or",
                                       "",     "  \"zone_name\" for starting nodes within zone");


// auxiliary definitions
namespace aux {


// boost graph library shortcuts
using boost::get;
using boost::vertex_color;
using boost::vertex_degree;
using boost::vertex_index;
using boost::vertex_index_t;


// renumbering method skeleton
struct renumber {
  virtual ~renumber() {}
  virtual unsigned calc(
    // outputs: bandwidth, renumbering maps new-to-old and old-to-new,
    // inputs: start node (negative for automatic) and mesh graph
    vector< unsigned >& B2A, vector< unsigned >& A2B,
    const int start, Graph& G ) = 0;
  virtual void adjustperm(vector< unsigned >& perm) {}
};


// renumbering method: Cuthill-McKee
struct ren_rcm : renumber {
  unsigned calc( vector< unsigned >& B2A, vector< unsigned >& A2B,
    const int start, Graph& G )
  {
    // calculate permutation (B2A)
    for (unsigned i=0; i<B2A.size(); ++i)
      B2A[i] = i;
    if (start<0) {
      boost::cuthill_mckee_ordering( G,
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G) );
    }
    else {
      boost::cuthill_mckee_ordering( G,boost::vertex((unsigned) start,G),
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G) );
    }
    adjustperm(B2A);

    // calculate inverse permutation (A2B)
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;

    // return bandwidth
    return boost::bandwidth( G,
        boost::make_iterator_property_map(&A2B[0],get(vertex_index,G),A2B[0]) );
  }
};


// renumbering method: King
struct ren_rking : renumber {
  unsigned calc( vector< unsigned >& B2A, vector< unsigned >& A2B,
    const int start, Graph& G )
  {
    // calculate permutation (B2A)
    for (unsigned i=0; i<B2A.size(); ++i)
      B2A[i] = i;
    boost::property_map< Graph,vertex_index_t >::type index_map =
      get(vertex_index,G);
    if (start<0) {
      boost::king_ordering( G,
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }
    else {
      boost::king_ordering( G,vertex((unsigned) start,G),
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }
    adjustperm(B2A);

    // calculate inverse permutation (A2B)
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;

    // return bandwidth
    return boost::bandwidth( G,
        boost::make_iterator_property_map(&A2B[0],index_map,A2B[0]) );
  }
};


// renumbering methods (reversed): Cuthill-McKee, King
struct ren_cm    : ren_rcm    { void adjustperm(vector< unsigned >& perm) { std::reverse(perm.begin(),perm.end()); } };
struct ren_king  : ren_rking  { void adjustperm(vector< unsigned >& perm) { std::reverse(perm.begin(),perm.end()); } };


}  // namespace aux
// auxiliary definitions


void t_ren::transform(GetPot& o, mmesh& m)
{
  const string o_method = o.get(o.inc_cursor(),"");
        string o_node   = o.get(o.inc_cursor(),"");


  cout << "info: building graph..." << endl;
#if 1
  const unsigned Nnodes = m.n();
  Graph G(Nnodes);

  {
    vector< vector< unsigned > > M(Nnodes);

    // use all available zones to get all nodes neighbors
    for (unsigned iz=0; iz<m.z(); ++iz) {
      const vector< melem >& e2n = m.vz[iz].e2n;
      for (unsigned j=0; j<e2n.size(); ++j) {
        const vector< unsigned >& nodes = e2n[j].n;
        const unsigned _Nnodes = nodes.size();
        for (unsigned k=0; k<_Nnodes; ++k)
          for (unsigned l=1; l<_Nnodes; ++l) {  // guarantee "except itself"
            const unsigned I = nodes[k];
            const unsigned J = nodes[(k+l)%_Nnodes];
            M[I].push_back(J);
            M[J].push_back(I);
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

  if (!boost::num_edges(G)) {
    cerr << "error: couldn't build graph!" << endl;
    throw 42;
  }

  {
    vector< unsigned > vcomp(Nnodes);
    const int Ncomp = boost::connected_components(G,&vcomp[0]);
    o_node = Ncomp>1? "-1" : o_node;
    cout << "info: number of components: " << Ncomp << (Ncomp>1? " (automatic start node forced)":"") << endl;
  }
  cout << "info: building graph." << endl;


  vector< unsigned > new2old(Nnodes,0),
                     old2new(Nnodes,0);
  for (unsigned i=0; i<Nnodes; ++i)
    new2old[i] = old2new[i] = i;

  ofstream f("dump.sparsity.plt",ios_base::trunc);
  f << "VARIABLES = i j" << endl;
  dump(f,"original",G,old2new);
  if (o_method=="dump")
    return;


  cout << "info: set method..." << endl;
  aux::renumber* renumber = NULL;
  {
    mfactory< aux::renumber >* fac = mfactory< aux::renumber >::instance();
    fac->Register(new ProductionLine< aux::renumber,aux::ren_cm     >("cm",    "Cuthill-McKee"));
    fac->Register(new ProductionLine< aux::renumber,aux::ren_king   >("king",  "King"));
    fac->Register(new ProductionLine< aux::renumber,aux::ren_rcm    >("rcm",   "reversed Cuthill-McKee"));
    fac->Register(new ProductionLine< aux::renumber,aux::ren_rking  >("rking", "reversed King"));
    renumber = fac->Create(o_method);
  }
  if (renumber==NULL) {
    cerr << "error: incorrect method" << endl;
    throw 42;
  }
  cout << "info: set method." << endl;


  int start    = (int) Nnodes;
  int fromzone = -1;
  if (o_node=="b") {
    unsigned min_b = Nnodes;
    cout << "info: finding minimum bandwidth..." << endl;
    boost::progress_display pbar(Nnodes);
    for (int s=0; s<(int) Nnodes; ++s, ++pbar) {
      const unsigned calc_b = renumber->calc(new2old,old2new,s,G);
      if (calc_b<min_b) {
        min_b = calc_b;
        start = s;
      }
    }
    cout << "info: found s:" << start << " bw: " << min_b << endl;
    cout << "info: finding minimum bandwidth." << endl;
  }
  else if (o_node=="d") {
    unsigned min_d = Nnodes;
    cout << "info: finding minimum degree..." << endl;
    boost::progress_display pbar(Nnodes);
    for (unsigned s=0; s<Nnodes; ++s, ++pbar) {
      const unsigned deg = boost::degree(boost::vertex(s,G),G);
      if (deg<min_d) {
        min_d = deg;
        start = (int) s;
      }
    }
    cout << "info: found s: " << start << " degree:" << min_d << endl;
    cout << "info: finding minimum degree." << endl;
  }
  else {
    // try zone names to see if one matches
    for (unsigned iz=0; iz<m.z(); ++iz)
      if (o_node==m.vz[iz].n)
        fromzone = (int) iz;
    if (fromzone>=0) {

      cout << "info: finding minimum bandwidth (in zone)..." << endl;
      unsigned min_b = Nnodes;

      cout << "info: building zone node set..." << endl;
      std::set< unsigned > nodeset;
      for (unsigned ie=0; ie<m.e(fromzone); ++ie)
        for (std::vector< unsigned >::const_iterator n=m.vz[fromzone].e2n[ie].n.begin(); n!=m.vz[fromzone].e2n[ie].n.end(); ++n)
          nodeset.insert(*n);
      cout << "info: building zone node set." << endl;

      boost::progress_display pbar((unsigned long) nodeset.size());
      for (std::set< unsigned >::const_iterator n=nodeset.begin(); n!=nodeset.end(); ++n, ++pbar) {
        const unsigned calc_b = renumber->calc(new2old,old2new,*n,G);
        if (calc_b<min_b) {
          min_b = calc_b;
          start = *n;
        }
      }
      cout << "info: found s:" << start << " bw: " << min_b << endl;
      cout << "info: finding minimum bandwidth (in zone)." << endl;

    }
    else {
      // if no zones match, use given index instead
      istringstream is(o_node);
      is >> start;
    }
  }


  if (fromzone<0 && start>=(int) Nnodes) {
    cout << "warn: can't find start node, not renumbering." << endl;
    return;
  }


  cout << "info: renumbering..."  << endl;
  {
    // calculate bandwidth and (re)calculate mappings
    const unsigned orig_b = boost::bandwidth(G),
                   calc_b = renumber->calc(new2old,old2new,start,G);

    // apply mapping if improvement is possible
    ostringstream os;
    os << "renumbering  s:" << start << "  bw:" << orig_b << ">" << calc_b;
    cout << "info: " << os.str() << endl;
    if (calc_b>=orig_b)
      cout << "warn: bandwidth not improved, but renumbering anyway." << endl;
    apply(old2new,m);
    dump(f,os.str(),G,old2new);
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
      for (unsigned k=0; k<Nenodes; ++k)
        for (unsigned l=0; l<Ndim; ++l)
          c[k][l] = m.vv[l][ enodes[k] ];
      if (t_geo::cellgeom_ok(c))
        continue;
      ++Nerrors;
      for (unsigned l=0; l<Nenodes; ++l) {
        ++Nswap;
        swap( enodes[l], enodes[(l+1)%Nenodes] );
        for (unsigned k=0; k<Nenodes; ++k)
          for (unsigned d=0; d<Ndim; ++d)
            c[k][d] = m.vv[d][ enodes[k] ];
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

