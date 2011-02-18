#ifndef t_plas_h
#define t_plas_h

#include "mkernel.h"
#include "plas.h"


// Preprocessor constants
#define AIR      1
#define WATER    2
#define NITROGEN 3


// module to use the Particle Lagrangian Solver (PLaS)
class t_plas : public m::mtransform,
               public plas {
 public:
  void transform(GetPot& o, m::mmesh& m);

 private:  // member functions
  double plasdriver_CalcAreaTriangle(double c[3][2]);
  double plasdriver_CalcVolumeTetra(double c[4][3]);
  void   plasdriver_InitFlowField(int material);
  void   plasdriver_GetFaceNodes(int elm, int face, int *nodes);

  void   plasdriver_ReadGambitNeutralFile();
  void   plasdriver_CalcElmsAroundNode();
  void   plasdriver_CalcElementNeighbours();
  void   plasdriver_CalcElementNormals();
  void   plasdriver_CalcElementVolumes();
  void   plasdriver_CalcNodalVolumes();

 private:  // plas interface functions
   void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);
   void setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);
   double getBndFaceRefCoord           (int bnd, int bface, int dim);
   double getElmFaceMiddlePoint        (int elm, int eface, int dim);
   double getBndFaceNormComp           (int bnd, int face, int dim)    { return dmesh.elmNorms[dmesh.bndDomElms[bnd][face]][dmesh.bndFaces[bnd][face]][dim]; }
   double getElmNormComp               (int elm, int eface, int dim)   { return dmesh.elmNorms[elm][eface][dim]; }
   double getElmVol                    (int elm)                       { return dmesh.elmVolumes[elm]; }
   double getEulerianTimeScale         (int nod)                       { return 0.; }
   double getNodCoord                  (int nod, int dim)              { return M->vv[dim][nod]; }
   double getNodVol                    (int nod)                       { return dmesh.nodVolumes[nod]; }
   double getPerBndOffset              (int bnd, int dim)              { return 0.; }  // TODO: implement
   double getPressure                  (int nod)                       { return dflow.p[nod]; }
   double getPressureOld               (int nod)                       { return dflow.p[nod]; }
   double getTemperature               (int nod)                       { return dflow.T[nod]; }
   double getTemperatureOld            (int nod)                       { return dflow.T[nod]; }
   double getVelocityComp              (int nod, int dim)              { return dflow.u[nod][dim]; }
   double getVelocityCompOld           (int nod, int dim)              { return dflow.u[nod][dim]; }
   double getVelocityDerivativeComp    (int nod, int idim, int jdim)   { return 0.; }
   double getVelocityDerivativeCompOld (int nod, int idim, int jdim)   { return 0.; }
   int    getBndDomElm                 (int bnd, int bface, int dummy) { return dmesh.bndDomElms[bnd][bface]; }
   int    getElementType               (int elm)                       { return dmesh.elmTypes[elm]; }
   int    getElmNeighbour              (int elm, int eface)            { return dmesh.elmNeighbs[elm][eface]; }
   int    getElmNode                   (int elm, int enod)             { return dmesh.elmNodes[elm][enod]; }
   int    getNumBndFaces               (int bnd)                       { return dmesh.numBndFaces[bnd]; }
   int    getPerBndFlag                (int bnd)                       { return dmesh.bndTypes[bnd]<0?  1:0; }
   int    getWallBndFlag               (int bnd)                       { return dmesh.bndTypes[bnd]==1? 1:0; }
   int    EndElementSearch             (double *pos)                   { return dmesh.numElm-1; }
   int    StartElementSearch           (double *pos)                   { return 0; }
   void   screenOutput                 (const std::string& text)                    { std::cout << "t_plas: info: " << text << std::endl; }
   void   screenWarning                (const std::string& text)                    { std::cout << "t_plas: warn: " << text << std::endl; }

 private:  // member variables
  m::mmesh *M;

  std::vector< bool > m_ziswall;

#if 0
  struct s_zoneprops {
    s_zoneprops() : nelem(0), etype(m::ORDERED) {}
    unsigned nelem;
    m::mtype etype;
  };
  std::vector< s_zoneprops >
    m_zinner,
    m_zbound;
#endif


  // mesh paramaters data structure
  struct {
    int
      numElm,
      *elmTypes,
      *numElmNodes,
      **elmNodes,
      *numNodElms,
      **nodElms,
      *numElmFaces,
      **elmNeighbs,
      *bndTypes,
      *numBndFaces,
      **bndFaces,
      **bndDomElms;
    double
      ***elmNorms,
      *nodVolumes,
      *elmVolumes,
      minElmVolume,
      maxElmVolume;
  } dmesh;

  // flow field data structure
  struct {
    double
      rho,
      mu,
      cp,
      k,
      *p,
      **u,
      *T;
  } dflow;

  // input parameters from drivers.conf file
  struct {
    int
      iter,
      numIter,
      numUnk,
      *bnd,
      material;
    double dt;
  } dparam;

};


#endif
