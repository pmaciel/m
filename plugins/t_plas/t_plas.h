#ifndef t_plas_h
#define t_plas_h

#include "mkernel.h"
#include "plas.h"


// Preprocessor constants
#define AIR      1
#define WATER    2
#define NITROGEN 3


// mesh parameters
struct s_driver_mesh {
  std::vector< std::vector< int > > bndFaces;    // <---- remove!
  std::vector< std::vector< int > > bndDomElms;  // <---- remove!

  std::vector< std::vector< int > > nodElms;
  std::vector< std::vector< int > > elmNeighbs;
  std::vector< double > nodVolumes;
  std::vector< double > elmVolumes;
  std::vector< std::vector< std::vector< double > > > elmNorms;
};


// parameters and flow properties
struct s_driver_param {
  double
    rho,
    mu,
    cp,
    k,
    dt;
  int
    iter,
    numIter,
    material;
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

 private:  // member functions
  double plasdriver_CalcAreaTriangle(double c[3][2]);
  double plasdriver_CalcVolumeTetra(double c[4][3]);
  void   plasdriver_InitFlowField(int material);
  void   plasdriver_GetFaceNodes(int iz, int ie, int face, int *nodes);

 private:  // plas interface functions
   void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);
   void setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);

   double getBndFaceRefCoord(int bnd, int bface, int dim) {
     return M->vv[dim][ M->vz[bnd].e2n[bface].n[0] ];
   }

   double getElmFaceMiddlePoint(int elm, int eface, int dim) {
     double coord = 0.;
     for (int iz=0, nelems=0; iz<(int) m_zinner_props.size(); ++iz) {
       nelems += m_zinner_props[iz].nelems;
       if (nelems > elm) {

         const int ie = elm - nelems + m_zinner_props[iz].nelems;
         int
           faceNodes[4],
           ctr = 0;
         plasdriver_GetFaceNodes(iz,ie,eface,faceNodes);

         for (int i=0; i<4; ++i)
           if (faceNodes[i]!=-1) {
             coord += M->vv[dim][faceNodes[i]];
             ++ctr;
           }
         coord /= ctr;

         break;
       }
     }
     return coord;
   }

   double getBndFaceNormComp(int bnd, int face, int dim) {
     return dmesh.elmNorms
         [ dmesh.bndDomElms[bnd][face] ]
         [ dmesh.bndFaces[bnd][face] ]
         [ dim ];
   }
   double getElmNormComp               (int elm, int eface, int dim) { return dmesh.elmNorms[elm][eface][dim]; }
   double getElmVol                    (int elm)                     { return dmesh.elmVolumes[elm]; }
   double getEulerianTimeScale         (int nod)                     { return 0.; }
   double getNodVol                    (int nod)                     { return dmesh.nodVolumes[nod]; }
   double getQuantity   (const PLAS_QUANTITY& Q, int nod) { return (m_quantity_idx   [Q]<0? 0. : ( M->vv[ m_quantity_idx   [Q] ][nod] )); }
   double getOldQuantity(const PLAS_QUANTITY& Q, int nod) { return (m_quantityold_idx[Q]<0? 0. : ( M->vv[ m_quantityold_idx[Q] ][nod] )); }
   int    getBndDomElm                 (int bnd, int bface)          { return dmesh.bndDomElms[bnd][bface]; }
   int    getElementType               (int elm) {
     int nelems = 0;
     for (std::vector< s_zoneprops >::const_iterator z=m_zinner_props.begin(); z!=m_zinner_props.end(); ++z) {
       nelems += z->nelems;
       if (elm<nelems)
         return z->e_type;
     }
     return 0;
   }
   int    getElmNeighbour              (int elm, int eface)          { return dmesh.elmNeighbs[elm][eface]; }
   int getElmNode(int elm, int enod) {
     int nelems = 0;
     for (size_t iz=0; iz<m_zinner_props.size(); ++iz) {
       if (elm < nelems + m_zinner_props[iz].nelems)
         return M->vz[iz].e2n[ elm-nelems ].n[enod];
       else
         nelems += m_zinner_props[iz].nelems;
     }
     return -1;
   }
   int    getNumBndFaces               (int bnd)                     { return m_zbound_props[bnd].nelems; }
   int    getWallBndFlag               (int bnd)                     { return m_zbound_props[bnd].iswall? 1:0; }
   void   screenOutput                 (const std::string& text)     { std::cout << "t_plas: info: " << text << std::endl; }
   void   screenWarning                (const std::string& text)     { std::cout << "t_plas: warn: " << text << std::endl; }

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
