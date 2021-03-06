#ifndef t_ren_h
#define t_ren_h

#include "boost/graph/adjacency_list.hpp"
#include "mkernel.h"


// graph definition
typedef boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
    boost::property< boost::vertex_color_t,  boost::default_color_type,
    boost::property< boost::vertex_degree_t, unsigned > >
  > RenumberingGraph;


// renumber nodes using boost graph library
class t_ren : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
  void dump(std::ostream& o, const std::string& t, const RenumberingGraph& G, const std::vector< unsigned >& A2B);
 private:
  void apply(const std::vector< unsigned >& A2B, m::mmesh& m);
};


#endif

