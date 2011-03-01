#ifndef t_plas_h
#define t_plas_h

#include "mkernel.h"
#include "plas.h"


// mesh parameters
struct s_driver_mesh {
  std::vector< std::vector< int > > bndFaces;    // <---- remove!
  std::vector< std::vector< int > > bndDomElms;  // <---- remove!

  std::vector< std::vector< int > > nodElms;
  std::vector< std::vector< int > > elmNeighbs;
  std::vector< double > nodVolumes;
  std::vector< std::vector< std::vector< double > > > elmNorms;
};


// parameters and flow properties
struct s_driver_param {
  double dt;
  int
    iter,
    numIter;
};


// additional zone properties (for convenience)
struct s_zoneprops {
  s_zoneprops() : iswall(false), nelems(0), e_type(0), e_nnodes(0), e_nfaces(0) {}
  bool iswall;
  int
    nelems,
    e_type,
    e_nnodes,
    e_nfaces;
};


// module to use the Particle Lagrangian Solver (PLaS)
class t_plas : public m::mtransform,
               public plas {
 public:
  void transform(GetPot& o, m::mmesh& m);

 private:  // plas interface functions
   void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);
   void setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);

   int    getBndDomElm        (int bnd, int bface)                                      { return dmesh.bndDomElms[bnd][bface]; }
   double getBndFaceNormComp(int bnd, int face, int dim) {
     return dmesh.elmNorms
         [ dmesh.bndDomElms[bnd][face] ]
         [ dmesh.bndFaces[bnd][face] ]
         [ dim ];
   }
   int    getElmNeighbour     (int elm, int eface)                                      { return dmesh.elmNeighbs[elm][eface]; }
   void   getElmNodes         (const int iz, const int ie, std::vector< unsigned >& en) { en = M->vz[iz].e2n[ie].n; }
   double getElmNormComp      (int elm, int eface, int dim)                             { return dmesh.elmNorms[elm][eface][dim]; }
   int    getElmType          (const int iz, const int ie)                              { return (m_zinner_props[iz].nelems? m_zinner_props[iz].e_type : m_zbound_props[iz].e_type); }
   double getEulerianTimeScale(int nod)                                                 { return 0.; }
   double getNodVol           (int nod)                                                 { return dmesh.nodVolumes[nod]; }
   double getOldQuantity      (const PLAS_QUANTITY& Q, int nod)                         { return (m_quantityold_idx[Q]<0? 0. : ( M->vv[ m_quantityold_idx[Q] ][nod] )); }
   double getQuantity         (const PLAS_QUANTITY& Q, int nod)                         { return (m_quantity_idx   [Q]<0? 0. : ( M->vv[ m_quantity_idx   [Q] ][nod] )); }
   int    getWallBndFlag      (int bnd)                                                 { return m_zbound_props[bnd].iswall? 1:0; }
   void   screenOutput        (const std::string& text)                                 { std::cout << "t_plas: info: " << text << std::endl; }
   void   screenWarning       (const std::string& text)                                 { std::cout << "t_plas: warn: " << text << std::endl; }

 private:  // member variables
  m::mmesh *M;
  std::vector< s_zoneprops >
    m_zinner_props,
    m_zbound_props;
  std::vector< int >
    m_quantity_idx,
    m_quantityold_idx;

  // FIXME refactor!
  s_driver_mesh dmesh;
  s_driver_param dparam;

};


#endif
