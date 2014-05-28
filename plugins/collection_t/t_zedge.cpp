
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


  // find index of named zones
  // NOTE: done by index, for modifying the zones vector invalidates references/iterators
  vector< string > znames(utils::split(o.get(o.inc_cursor(),""),':'));
  for (vector< string >::const_iterator n=znames.begin(); n!=znames.end(); ++n) {
    try {

      // check if zones exist and comply with requirements (2D)
      const unsigned i(utils::getzindex(m,*n));
      if (m.d(i)!=2)
        std::cerr << "warn: zone \"" << *n << "\": not 2D, skipped" << std::endl;
        throw 42;

      // build set of unique edges
      // NOTE: avoid duplicates by attempting to remove them first
      std::set< aux::eedge, aux::eedge_lt > edges;
      for (vector< melem >::const_iterator eli=m.vz[i].e2n.begin(); eli!=m.vz[i].e2n.end(); ++eli) {

        vector< aux::eedge > edges_in;
        switch (m.vz[i].t) {

          case (m::FETRIANGLE):
            edges_in.push_back(std::make_pair(eli->n[0],eli->n[1]));
            edges_in.push_back(std::make_pair(eli->n[1],eli->n[2]));
            edges_in.push_back(std::make_pair(eli->n[2],eli->n[0]));
            break;

          case (m::FEQUADRILATERAL):
            edges_in.push_back(std::make_pair(eli->n[0],eli->n[1]));
            edges_in.push_back(std::make_pair(eli->n[1],eli->n[2]));
            edges_in.push_back(std::make_pair(eli->n[2],eli->n[3]));
            edges_in.push_back(std::make_pair(eli->n[3],eli->n[0]));
            break;

          case (m::FEPOLYGON):
            std::cerr << "warn: zone \"" << *n << "\": polygonal zones not supported" << std::endl;
            throw 42;
            break;

          default:
            break;

        }
        for (vector< aux::eedge >::const_iterator ei=edges_in.begin(); ei!=edges_in.end(); ++ei)
          if (!edges.erase(*ei))
            edges.insert(*ei);

      }

      // add surviving edges to a (new) destination zone connectivity table
      if (!edges.empty()) {
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
        std::cout << "info: zone \"" << *n << "\" edges: " << e2n.size() << std::endl;
      }
      else {
        std::cerr << "warn: zone \"" << *n << "\": no edges found" << std::endl;
      }

    }
    catch (const int&) {
      std::cerr << "warn: skipping zone \"" << *n << "\"" << std::endl;
    }
  }


  std::cout << "::zone edge." << std::endl;
}

