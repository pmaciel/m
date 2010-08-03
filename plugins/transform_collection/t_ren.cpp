
#include <list>
#include <sstream>
#include "mfactory.h"
#include "t_ren.h"
#include "t_geo.h"

#include "boost/graph/cuthill_mckee_ordering.hpp"
#include "boost/graph/king_ordering.hpp"
#include "boost/graph/minimum_degree_ordering.hpp"
#include "boost/graph/bandwidth.hpp"
#include "boost/graph/connected_components.hpp"


using namespace std;
using namespace m;


//reverse Cuthill-McKee algorithm
Register< mtransform,t_ren > mt_ren(11,"-tren","[str] [int] renumber nodes using boost graph library,",
                                       "",     "with [str] method:",
                                       "",     "  \"rcm\"/\"cm\", for (reversed/normal) Cuthill-McKee algorithm, or",
                                       "",     "  \"rk\"/\"k\", for (reversed/normal) King algorithm, or",
                                       "",     "  \"mindeg\", for minimum degree ordering algorithm, or",
                                       "",     "  \"dump\", for just dump sparsity (\"dump.sparsity.plt\"), and",
                                       "",     "with [int] node:",
                                       "",     "  [i] for starting at node i, or",
                                       "",     "  [-i] for automatic starting node, or",
                                       "",     "  \"b\" for minimum bandwidth nodes, or",
                                       "",     "  \"d\" for minimum degree node");


// auxiliary definitions
namespace aux {


// renumbering method skeleton
struct renumber {
  virtual ~renumber() {}
  virtual void calc(
    // outputs: (half) bandwidth, zone title (for dumping),
    // and renumbering maps, from new to old and old to new
    unsigned& bw, string& tag, vector< unsigned >& B2A, vector< unsigned >& A2B,
    // inputs: start node (negative for automatic) and mesh graph
    int start, Graph& G ) = 0;
};


// renumbering method: reverse Cuthill-McKee
struct ren_rcm : renumber {
  void calc(
    unsigned& bw, string& tag, vector< unsigned >& B2A, vector< unsigned >& A2B,
    int start, Graph& G )
  {
    using namespace boost;

    // calculate permutation
    if (start<0) {
      cuthill_mckee_ordering( G,
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G) );
    }
    else {
      cuthill_mckee_ordering( G,vertex((unsigned) start,G),
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G) );
    }

    // calculate inverse permutation
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;

    // calculate bandwidth and set tag
    property_map< Graph,vertex_index_t >::type index_map = get(vertex_index,G);
    bw = bandwidth(G,make_iterator_property_map(&A2B[0],index_map,A2B[0]));
    ostringstream s;
    s << "RCM[" << start << "] bandwidth: " << bw;
    tag = s.str();
  }
};


// renumbering method: non-reversed Cuthill-McKee
struct ren_cm : renumber {
  void calc(
    unsigned& bw, string& tag, vector< unsigned >& B2A, vector< unsigned >& A2B,
    int start, Graph& G )
  {
    using namespace boost;

    // calculate permutation
    if (start<0) {
      cuthill_mckee_ordering( G,
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G) );
    }
    else {
      cuthill_mckee_ordering( G,vertex((unsigned) start,G),
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G) );
    }

    // calculate inverse permutation
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;

    // calculate bandwidth and set tag
    property_map< Graph,vertex_index_t >::type index_map = get(vertex_index,G);
    bw = bandwidth(G,make_iterator_property_map(&A2B[0],index_map,A2B[0]));
    ostringstream s;
    s << "CM[" << start << "] bandwidth: " << bw;
    tag = s.str();
  }
};


// renumbering method: reversed King
struct ren_rk : renumber {
  void calc(
    unsigned& bw, string& tag, vector< unsigned >& B2A, vector< unsigned >& A2B,
    int start, Graph& G )
  {
    using namespace boost;

    // calculate permutation
    property_map< Graph,vertex_index_t >::type index_map = get(vertex_index,G);
    if (start<0) {
      king_ordering( G,
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }
    else {
      king_ordering( G,vertex((unsigned) start,G),
        B2A.rbegin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }

    // calculate inverse permutation
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;

    // calculate bandwidth and set tag
    bw = bandwidth(G,make_iterator_property_map(&A2B[0],index_map,A2B[0]));
    ostringstream s;
    s << "RK[" << start << "] bandwidth: " << bw;
    tag = s.str();
  }
};

// renumbering method: non-reversed King
struct ren_k : renumber {
  void calc(
    unsigned& bw, string& tag, vector< unsigned >& B2A, vector< unsigned >& A2B,
    int start, Graph& G )
  {
    using namespace boost;

    // calculate permutation
    property_map< Graph,vertex_index_t >::type index_map = get(vertex_index,G);
    if (start<0) {
      king_ordering( G,
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }
    else {
      king_ordering( G,vertex((unsigned) start,G),
        B2A.begin(),get(vertex_color,G),get(vertex_degree,G),index_map );
    }

    // calculate inverse permutation
    for (unsigned i=0; i<B2A.size(); ++i)
      A2B[B2A[i]] = i;

    // calculate bandwidth and set tag
    bw = bandwidth(G,make_iterator_property_map(&A2B[0],index_map,A2B[0]));
    ostringstream s;
    s << "K[" << start << "] bandwidth: " << bw;
    tag = s.str();
  }
};


// renumbering method: minimum degree ordering
struct ren_mindeg : renumber {
  void calc(
    unsigned& bw, string& tag, vector< unsigned >& B2A, vector< unsigned >& A2B,
    int start, Graph& G )
  {
    /*
     * start represents delta, meaning: "Multiple elimination control variable.
     * If it is larger than or equal to zero then multiple elimination is
     * enabled. The value of delta specifies the difference between the minimum
     * degree and the degree of vertices that are to be eliminated."
     */
    using namespace boost;

    static vector< unsigned > supernode_sizes(B2A.size());
    supernode_sizes.assign(B2A.size(),1);

    // calculate permutation (and inverse permutation)
    property_map< Graph,vertex_index_t >::type index_map = get(vertex_index,G);
    minimum_degree_ordering( G,
      get(vertex_degree,G),
      B2A.begin(),
      A2B.begin(),
      make_iterator_property_map(&supernode_sizes[0],index_map,supernode_sizes[0]),
      start,index_map );

    // calculate bandwidth and set tag
    bw = bandwidth(G,make_iterator_property_map(&A2B[0],index_map,A2B[0]));
    ostringstream s;
    s << "MinDeg[" << start << "] bandwidth: " << bw;
    tag = s.str();
  }
};


}  // namespace aux
// auxiliary definitions


void t_ren::transform(GetPot& o, mmesh& m)
{
  using namespace boost;


#if 1
  const unsigned Nnodes = m.n();
  Graph G(Nnodes);


  cout << "info: building graph..." << endl;
  // find a zone with "dimensionality" same as number of dimensions,
  // hopefully connecting all nodes in the mesh -- if not, this won't work well
  for (unsigned j=0; j<m.z(); ++j) {
    if (m.d(j)==m.d()) {
      vector< vector< unsigned > > M(Nnodes);
      const vector< melem >& e2n = m.vz[j].e2n;

      // get all nodes neighbors then remove duplicate entries and construct graph
      // 1* l=1 guarantees "except itself"
      // 2* graph is built with edges with higher indices than "itself node", to
      //    guarantee no parallel edges are introduced
      for (unsigned j=0; j<e2n.size(); ++j) {
        const vector< unsigned >& nodes = e2n[j].n;
        const unsigned Nnodes = nodes.size();
        for (unsigned k=0; k<Nnodes; ++k)
          for (unsigned l=1; l<Nnodes; ++l) {  //*1
            const unsigned i = nodes[k];
            const unsigned j = nodes[(k+l)%Nnodes];
            M[i].push_back(j);
            M[j].push_back(i);
          }
      }

      unsigned i=0;
      for (vector< vector< unsigned > >::iterator n=M.begin(); n!=M.end(); ++n, ++i) {
        sort(n->begin(),n->end());
        n->erase(unique(n->begin(),n->end()),n->end());
        for (vector< unsigned >::const_iterator j=n->begin(); j<n->end(); ++j)
          if (*j>i)  //*2
            add_edge(i,*j,G);
      }
      break;
    }
  }
  if (!num_edges(G)) {
    cerr << "error: didn't find appropriate zone with d=" << m.d() << endl;
    throw 42;
  }
  cout << "info: building graph." << endl;
#else

#if 1 // test graph 1: small graph
  enum NODES {N1,N2,N3,N4,N5,N6,N7,N8,Nnodes};
  Graph G(Nnodes);
  add_edge(N1,N5,G);
  add_edge(N5,N3,G);
  add_edge(N3,N2,G);
  add_edge(N2,N6,G);
  add_edge(N2,N8,G);
  add_edge(N6,N8,G);
  add_edge(N4,N7,G);
  add_edge(N4,N5,G); // connects its two components
#else // test graph 2: small graph, three disconnected components
  enum NODES {N1,N2,N3,N4,N5,N6,Nnodes};
  Graph G(Nnodes);
  add_edge(N1,N2,G);
  add_edge(N2,N5,G);
  add_edge(N5,N1,G);
  add_edge(N3,N6,G);
#endif

#endif


#if 0 // test just graphs connected components (and exit)
  vector< int > component(num_vertices(G));
  int num = connected_components(G,&component[0]);

  cout << "Total number of components: " << num << endl;
  for (size_t i=0; i!=component.size(); ++i)
    cout << "Vertex " << i <<" is in component " << component[i] << endl;
  cout << endl;
  exit(0);
#endif


  vector< unsigned > new2old(Nnodes);
  vector< unsigned > old2new(Nnodes);
  for (unsigned i=0; i<Nnodes; ++i)
    new2old[i] = old2new[i] = i;
  unsigned calc_b = Nnodes;
  unsigned orig_b = bandwidth(G);
  cout << "info: original bandwidth: " << orig_b << endl;

  ofstream f("dump.sparsity.plt",ios_base::trunc);
  f << "VARIABLES = i j" << endl;
  dump(f,"original",G,old2new);

  const string o_method = o.get(o.inc_cursor(),"");
  const string o_node   = o.get(o.inc_cursor(),"");
  if (o_method=="dump")
    return;


  cout << "info: set method..." << endl;
  aux::renumber* renumber = (o_method=="rcm"?    new aux::ren_rcm :
                            (o_method=="cm"?     new aux::ren_cm :
                            (o_method=="rk"?     new aux::ren_rk :
                            (o_method=="k"?      new aux::ren_k :
                            (o_method=="mindeg"? new aux::ren_mindeg :
                                                 (aux::renumber*) NULL )))));
  if (renumber==NULL) {
    cerr << "error: incorrect method" << endl;
    throw 42;
  }
  cout << "info: set method." << endl;


  int start = (int) Nnodes;
  string info;
  if (o_node=="b") {
    cout << "info: finding minumum bandwidth..." << endl;
    unsigned min_b = orig_b;
    string info;
    for (int s=0; s<(int) Nnodes; ++s) {
      renumber->calc(
        calc_b,info,new2old,old2new,
        s,G );
      cout << "info: " << info << (calc_b<min_b? "*":"") << endl;
      if (calc_b<min_b) {
        min_b = calc_b;
        start = s;
      }
    }
    cout << "info: node/bandwidth: " << start << '/' << min_b << endl;
    cout << "info: finding minumum bandwidth." << endl;
  }
  else if (o_node=="d") {
    cout << "info: finding minumum degree..." << endl;
    unsigned min_d = Nnodes;
    for (unsigned i=0; i<Nnodes; ++i) {
      const unsigned deg = degree(vertex(i,G),G);
      if (deg<min_d) {
        min_d = deg;
        start = (int) i;
      }
    }
    cout << "info: node/degree: " << start << '/' << min_d << endl;
    cout << "info: finding minumum degree..." << endl;
  }
  else {
    // use given index instead
    istringstream is(o_node);
    is >> start;
  }


  if (start<(int) Nnodes) {
    cout << "info: build renumbering maps, using start node: " << start << endl;
    renumber->calc(
      calc_b,info,new2old,old2new,
      start,G );
    cout << "info: " << info << endl;
    dump(f,info,G,old2new);
    cout << "info: build renumbering maps." << endl;

    if (calc_b<orig_b) {
      cout << "info: apply renumbering map..." << endl;
      apply(old2new,m);
      cout << "info: apply renumbering map." << endl;
      return;
    }
  }


  cout << "info: didn't improve existing bandwidth, not renumbering." << endl;
}


void t_ren::dump(ostream& o, const string& t, const Graph& G, const vector< unsigned >& A2B)
{
  using namespace boost;
  o << "ZONE T=\"" << t << "\" I=" << num_edges(G)*2 + num_vertices(G) << endl;

  // upper/lower diagonal then diagonal entries
  graph_traits< Graph >::edge_iterator e, ef;
  for (tie(e,ef) = edges(G); e!=ef; ++e) {
    const int u = (int) A2B[ source(*e,G) ];
    const int v = (int) A2B[ target(*e,G) ];
    o << u << ' ' << -v << endl
      << v << ' ' << -u << endl;
  }
  graph_traits< Graph >::vertex_iterator v, vf;
  for (tie(v,vf) = vertices(G); v!=vf; ++v) {
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
  cout << "info: elements/errors: " << Nelem << '/' << Nerrors << endl;
  cout << "info: swaps/corrections: " << Nswap << '/' << Ncorr << endl;
}

