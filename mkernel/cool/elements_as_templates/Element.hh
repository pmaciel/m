#ifndef COOLFluiD_Muffin_Element_hh
#define COOLFluiD_Muffin_Element_hh

#include "cool/Framework/Node.hh"
#include "cool/Framework/State.hh"

/*
 * the internal data for the simplicial elements follows face-opposite-node
 * representation. COOLFluiD representation follows Gambit; this representation
 * is accessed through cf_-preprended member functions.
 *
 * these functions also provide a placeholder for interfacing representations.
 */

namespace COOLFluiD {
  namespace Muffin {


// Possible element types
enum elem_type { ORDERED=0,
                 FELINESEG=1,
                 FETRIANGLE=2,
                 //FEQUADRILATERAL=3,
                 FETETRAHEDRON=4//,
                 //FEBRICK=5,
                 //FEPOLYGON=6,
                 //FEPOLYHEDRON=7
                 };


/// Generic element geometric properties, by dimensions and type
template< int D, int T >
struct AElement {
#define NOTIMPLEMENTED cf_always_assert_desc("not implemented",false)

  // constructor, destructor
  AElement() { NOTIMPLEMENTED; }
  ~AElement() {}

  // element and faces properties
  CFreal element(std::vector< Framework::Node* >& n) { NOTIMPLEMENTED; return 0.; }
  template< int F >
  AElement< D,F >& face(const CFuint& i) { NOTIMPLEMENTED; return (*new AElement< D,F >()); }

  // state derivative (on CFreal, given states and internal states)
  CFreal dd(const CoordXYZ& C, const std::vector< CFreal >& v);
  std::vector< CFreal > dd(const CoordXYZ& C, const std::vector< Framework::State* >& v);
  std::vector< CFreal > dd(const CoordXYZ& C) {
    std::vector< Framework::State* > v(S);
    for (CFuint i=0; i<S; ++i)
      v[i] = &states[i];
    return dd(C,v);
  }

  // COOLFluiD representation interface
  CFreal cf_norm2(const CFuint& i);
  RealVector& cf_normal(const CFuint& i);
  template< int F >
  Element< D,F >& cf_face(const CFuint& i) { NOTIMPLEMENTED; return (*new Element< D,F >()); }

  // member variables
  const CFuint N;  // number of nodes
  const CFuint S;  // ... states
  const CFuint F;  // ... faces
  CFreal s;                                // element size
  std::vector< Framework::Node > nodes;    // ... nodes
  std::vector< Framework::State > states;  // ... states
  std::vector< RealVector > normal;  // inwards-facing normals
  std::vector< CFreal > norm2;       // normals length squared

#undef NOTIMPLEMENTED
};


/// Generic P1-P1 simplex state derivative (N CFreal)
template< int D, int T >
CFreal AElement< D,T >::dd(const CoordXYZ& C, const std::vector< CFreal >& v)
{
  cf_assert_desc("unexpected number of states",v.size()==N);
  cf_assert_desc("unexpected derivative component",C<D);
  const CFreal dNface(N-1);

  CFreal r = 0;
  for (size_t i=0; i<N; ++i)
    r += v[i]*normal[i][C];
  return r/(dNface*s);
}


/// Generic P1-P1 simplex state derivative (N given states)
template< int D, int T >
std::vector< CFreal > AElement< D,T >::dd(const CoordXYZ& C, const std::vector< Framework::State >& v)
{
  cf_assert_desc("unexpected number of states",v.size()==S);
  cf_assert_desc("unexpected derivative component",C<D);
  const CFuint nbVars = v[0].size();
  cf_assert_desc("unexpected state size: 0",nbVars>0);
  const CFreal dNface(N-1);

  std::vector< CFreal > r(S,0.);
  for (size_t j=0; j<nbVars; ++j) {
    for (size_t i=0; i<N; ++i) {
      cf_assert_desc("unexpected state size",v[i].size()==nbVars);
      r[j] += v[i][j]*normal[i][C];
    }
    r[j] /= (dNface*s);
  }
  return r;
}


/// Generic P1-P1 simplex COOLFluiD/Gambit representation interface
template< int D, int T >
CFreal AElement< D,T >::cf_norm2(const CFuint& i)
{
  return norm2[ N==1? i
              :(N==2? i
              :(N==3? (i==0? 2:(i==1? 0:(i==2? 1:-1)))
              :(N==4? (i==1? 2:(i==2? 0:(i==3? 1:(i==0? 3:-1))))
                   : -1 ))) ];
}


/// Generic P1-P1 simplex COOLFluiD/Gambit representation interface
template< int D, int T >
RealVector& AElement< D,T >::cf_normal(const CFuint& i)
{
  return normal[ N==1? i
               :(N==2? i
               :(N==3? (i==0? 2:(i==1? 0:(i==2? 1:-1)))
               :(N==4? (i==1? 2:(i==2? 0:(i==3? 1:(i==0? 3:-1))))
                    : -1 ))) ];
}


  }  // namespace Muffin
}  // namespace COOLFluiD


// specialization of 2 and 3-dimensional elements
#include "Element2D.ci"
#include "Element3D.ci"


#endif // COOLFluiD_Muffin_Element_hh

