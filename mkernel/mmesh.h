#ifndef mmesh_h
#define mmesh_h

#include <vector>
#include <string>


namespace m {


// possible element/zone types
// (keeping simple by assuming zone has just one element type,
// it is stored as member of mzone)
enum mtype {
  ORDERED=0,  // equivalent to none
  FELINESEG=1,
  FETRIANGLE=2,
  FEQUADRILATERAL=3,
  FETETRAHEDRON=4,
  FEBRICK=5,
  FEPOLYGON=6,
  FEPOLYHEDRON=7,
  PRISM3=101,   // (introduced) prism, triangular base (or wedge)
  PYRAMID4=102  // (introduced) pyramid, quadrangular base
};


// description of an element
struct melem {
  std::vector< unsigned > n;
};


// description of a zone
struct mzone {
  std::string n;             // name
  std::vector< melem > e2n;  // element-node connectivity
  mtype                t;    // element/zone type

  // dimensionality
  unsigned d() const {
    return (t==FELINESEG?       1 :
           (t==FETRIANGLE?      2 :
           (t==FEQUADRILATERAL? 2 :
           (t==FETETRAHEDRON?   3 :
           (t==FEBRICK?         3 :
           (t==FEPOLYGON?       2 :
           (t==FEPOLYHEDRON?    3 :
           (t==PRISM3?          3 :
           (t==PYRAMID4?        3 :
                                0 )))))))));
  }
};


// description of a mesh
struct mmesh {

  // mesh dimensions and vector fields detection
  unsigned d() const;
  std::vector< bool > vvectors() const;

  // mesh operations:
  // merge: merge another mesh, into single point cloud
  // compress: apply point cloud simplification and renumbering
  // extract: compressed mesh structure with just this zone)
  // edge: detect the edge of a given zone
  void merge(const mmesh& another);
  void compress();
  mmesh extract(const std::string& zn);
  mmesh edge(const std::string& zn);

  // zones;
  // dimensionality, nb. elements and type
  std::vector< mzone > vz;
  inline unsigned d(unsigned i)   const { return i<z()? vz[i].d():0;        }
  inline unsigned e(unsigned i=0) const { return i<z()? vz[i].e2n.size():0; }
  inline unsigned t(unsigned i)   const { return i<z()? vz[i].t:ORDERED;    }

  // point cloud names and values;
  // nb. nodes, variables (inc. dimensions) and zones
  std::vector< std::string >           vn;
  std::vector< std::vector< double > > vv;
  inline unsigned n() const { return vv.size()? vv[0].size():0; }
  inline unsigned v() const { return vn.size(); }
  inline unsigned z() const { return vz.size(); }

};


}  // namespace m


#endif

