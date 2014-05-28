
#include <algorithm>
#include <sstream>
#include <set>

#include "mfactory.h"

#include "utils.h"
#include "t_zedge.h"


using namespace m;


Register< mtransform,t_zedge > mt_zedge("-tzedge","[str:...] zone edge extraction (2D zones only)");


// auxiliary definitions
namespace aux {

  // edge definition (includes
  typedef std::pair< unsigned, unsigned > eedge;

  struct eedge_lt {
    bool operator()(const eedge& e1, const eedge& e2) const {
      const unsigned
        e1s1(std::min(e1.first,e1.second)),
        e2s1(std::min(e2.first,e2.second));
      return e1s1<e2s1? true  :
             e1s1>e2s1? false :
             std::max(e1.first,e1.second)<std::max(e2.first,e2.second);
    }
  };


}  // namespace aux


void t_zedge::transform(GetPot& o, mmesh& m)
{
  using std::string;
  using std::vector;
  std::cout << "::zone edge..." << std::endl;


  // append a edges zone for each input zone
  // NOTE: done by index, for modifying the zones vector invalidates references/iterators
  vector< string > znames(utils::split(o.get(o.inc_cursor(),""),':'));
  for (vector< string >::const_iterator n=znames.begin(); n!=znames.end(); ++n) {
    unsigned idx;
    try                { idx = utils::getzindex(m,*n); }
    catch (const int&) { continue; }
    mzone& z = m.vz[idx];

    // build set of the zone unique edges
    // (if an element edge exists already, remove it instead)
    std::set< aux::eedge, aux::eedge_lt > edges;
    for (vector< melem >::const_iterator eli=z.e2n.begin(); z.d()==2 && eli!=z.e2n.end(); ++eli) {
      vector< aux::eedge > eedges;
      switch (z.t) {
        case (m::FETRIANGLE):
          eedges.push_back(std::make_pair(eli->n[0],eli->n[1]));
          eedges.push_back(std::make_pair(eli->n[1],eli->n[2]));
          eedges.push_back(std::make_pair(eli->n[2],eli->n[0]));
          break;
        case (m::FEQUADRILATERAL):
          eedges.push_back(std::make_pair(eli->n[0],eli->n[1]));
          eedges.push_back(std::make_pair(eli->n[1],eli->n[2]));
          eedges.push_back(std::make_pair(eli->n[2],eli->n[3]));
          eedges.push_back(std::make_pair(eli->n[3],eli->n[0]));
          break;
        case (m::FEPOLYGON):
        default:
          break;
      }
      for (vector< aux::eedge >::const_iterator edi=eedges.begin(); edi!=eedges.end(); ++edi)
        edges.erase(*edi)? true : edges.insert(*edi).second;
    }

    // add surviving edges to a (new) destination zone connectivity table
    if (edges.empty()) {
      std::cerr << "info: zone \"" << *n << "\" edges: 0 (only non-polygonal 2D zones are supported)" << std::endl;
    }
    else {
      m.vz.push_back(mzone(*n+"_edge",m::FELINESEG));
      vector< melem > &e2n = m.vz.back().e2n;
      melem e;
      e.n.resize(2);
      e2n.reserve(edges.size());
      for (std::set< aux::eedge, aux::eedge_lt >::const_iterator ei=edges.begin(); ei!=edges.end(); ++ei) {
        e.n[0] = ei->first;
        e.n[1] = ei->second;
        e2n.push_back(e);
      }
      std::cout << "info: zone \"" << *n << "\" edges: " << edges.size() << std::endl;
    }

  }


  std::cout << "::zone edge." << std::endl;
}

