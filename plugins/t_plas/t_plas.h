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

  void   plasdriver_ReadGambitNeutralFile(const XMLNode& x);
  void   plasdriver_CalcElmsAroundNode();
  void   plasdriver_CalcElementNeighbours();
  void   plasdriver_CalcElementNormals();
  void   plasdriver_CalcElementVolumes();
  void   plasdriver_CalcNodalVolumes();

 private:  // plas interface functions
   void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);
   void setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);

   double getBndFaceRefCoord(int bnd, int bface, int dim) {
     return M->vv[dim][ M->vz[bnd].e2n[bface].n[0] ];
   }

   double getElmFaceMiddlePoint(int elm, int eface, int dim) {
     int faceNodes[4];
     plasdriver_GetFaceNodes(elm,eface,faceNodes);

     int ctr = 0;
     double coord = 0.;
     for (int ifac=0; ifac<4; ++ifac)
       if (faceNodes[ifac]!=-1) {
         coord += M->vv[dim][faceNodes[ifac]];
         ++ctr;
       }
     coord /= ctr;

     return coord;
   }

   double getBndFaceNormComp           (int bnd, int face, int dim)  { return dmesh.elmNorms[dmesh.bndDomElms[bnd][face]][dmesh.bndFaces[bnd][face]][dim]; }
   double getElmNormComp               (int elm, int eface, int dim) { return dmesh.elmNorms[elm][eface][dim]; }
   double getElmVol                    (int elm)                     { return dmesh.elmVolumes[elm]; }
   double getEulerianTimeScale         (int nod)                     { return 0.; }
   double getNodCoord                  (int nod, int dim)            { return M->vv[dim][nod]; }
   double getNodVol                    (int nod)                     { return dmesh.nodVolumes[nod]; }
   double getPressure                  (int nod)                     { return dflow.p[nod]; }
   double getPressureOld               (int nod)                     { return dflow.p[nod]; }
   double getTemperature               (int nod)                     { return dflow.T[nod]; }
   double getTemperatureOld            (int nod)                     { return dflow.T[nod]; }
   double getVelocityComp              (int nod, int dim)            { return dflow.u[nod][dim]; }
   double getVelocityCompOld           (int nod, int dim)            { return dflow.u[nod][dim]; }
   double getVelocityDerivativeComp    (int nod, int idim, int jdim) { return 0.; }
   double getVelocityDerivativeCompOld (int nod, int idim, int jdim) { return 0.; }
   int    getBndDomElm                 (int bnd, int bface)          { return dmesh.bndDomElms[bnd][bface]; }
   int    getElementType               (int elm)                     { return dmesh.elmTypes[elm]; }
   int    getElmNeighbour              (int elm, int eface)          { return dmesh.elmNeighbs[elm][eface]; }
   int    getElmNode                   (int elm, int enod)           { return (int) M->vz[0].e2n[elm].n[enod]; }
   int    getNumBndFaces               (int bnd)                     { return dmesh.numBndFaces[bnd]; }
   int    getWallBndFlag               (int bnd)                     { return m_zprops[bnd].iswall? 1:0; }
   void   screenOutput                 (const std::string& text)     { std::cout << "t_plas: info: " << text << std::endl; }
   void   screenWarning                (const std::string& text)     { std::cout << "t_plas: warn: " << text << std::endl; }

 private:  // member variables
  m::mmesh *M;

  struct s_zoneprops {
    s_zoneprops() : isinner(false), iswall(false) {}
    bool isinner;
    bool iswall;
  };
  std::vector< s_zoneprops > m_zprops;

  // mesh paramaters data structure
  struct {
    int numElm;
    double minElmVolume, maxElmVolume;
    std::vector< int >
      elmTypes,
      numElmNodes,
      numNodElms,
      numElmFaces,
      numBndFaces;
    std::vector< std::vector< int > >
      nodElms,
      elmNeighbs,
      bndFaces,
      bndDomElms;
    std::vector< double >
      nodVolumes,
      elmVolumes;
    std::vector< std::vector< std::vector< double > > > elmNorms;
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
      material;
    double dt;
  } dparam;

};


#endif
