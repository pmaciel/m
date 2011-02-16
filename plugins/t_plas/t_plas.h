#ifndef t_plas_h
#define t_plas_h

#include "mkernel.h"
#include "plas.h"


// Preprocessor constants
#define AIR      1
#define WATER    2
#define NITROGEN 3


// Driver input parameters for drivers.conf file
typedef struct driver_parameters{
  int
    iter,
    numIter,
    numUnk,
    numBnd,
    *bnd,
    material;
  double dtEul;
  char gridString[100];
} DRIVER_PARAMETERS;


// Mesh paramaters data structure
typedef struct driver_gambit_mesh{
  int
    numNod,
    numElm,
    numBnd,
    numDim,
    numGrps,
    *elmTypes,
    *bndTypes,
    *numElmNodes,
    **elmNodes,
    *numBndFaces,
    **bndFaces,
    **bndDomElms,
    *numNodElms,
    **nodElms,
    *numElmFaces,
    **elmNeighbs;
  double
    **coords,
    ***elmNorms,
    *nodVolumes,
    *elmVolumes,
    domainVolume,
    minElmVolume,
    maxElmVolume;
  char
    **bndNames,
    text[100][100];
} DRIVER_GAMBIT_MESH;


// Flow field data structure
typedef struct driver_flow_field{
  double
    rho,
    mu,
    cp,
    k,
    *p,
    **u,
    *T;
} DRIVER_FLOW_FIELD;


// module to use the Particle Lagrangian Solver (PLaS)
class t_plas : public m::mtransform, plas {
 public:
  void transform(GetPot& o, m::mmesh& m);

  // member functions
 private:
  double plasdriver_CalcAreaTriangle(double c[3][2]);
  void   plasdriver_CalcElementNeighbours(DRIVER_GAMBIT_MESH *dmesh);
  void   plasdriver_CalcElementNormals(DRIVER_GAMBIT_MESH *dmesh);
  void   plasdriver_CalcElementVolumes(DRIVER_GAMBIT_MESH *dmesh);
  void   plasdriver_CalcElmsAroundNode(DRIVER_GAMBIT_MESH *dmesh);
  void   plasdriver_CalcNodalVolumes(DRIVER_GAMBIT_MESH *dmesh);
  double plasdriver_CalcVolumeTetra(double c[4][3]);
  void   plasdriver_FreeGambitMemory(DRIVER_GAMBIT_MESH *dmesh);
  void   plasdriver_InitFlowField(DRIVER_GAMBIT_MESH *dmesh, int material, DRIVER_FLOW_FIELD *dflow);
  void   plasdriver_ReadDriverDataFile(DRIVER_PARAMETERS *dparam);
  void   plasdriver_ReadGambitNeutralFile(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam);
  void   plasdriver_WriteTecplot(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam, DRIVER_FLOW_FIELD *dflow);
  void   plasdriver_GetFaceNodes(DRIVER_GAMBIT_MESH *dmesh, int elm, int face, int *nodes);

 private:  // plas interface functions
   void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);
   void setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);
   double getBndFaceRefCoord(int bnd, int bface, int dim);
   double getElmFaceMiddlePoint(int elm, int eface, int dim);
   double getBndFaceNormComp           (int bnd, int face, int dim)    { return dmesh.elmNorms[dmesh.bndDomElms[bnd][face]][dmesh.bndFaces[bnd][face]][dim]; }
   double getElmNormComp               (int elm, int eface, int dim)   { return dmesh.elmNorms[elm][eface][dim]; }
   double getElmVol                    (int elm)                       { return dmesh.elmVolumes[elm]; }
   double getEulerianTimeScale         (int nod)                       { return 0.; }
   double getNodCoord                  (int nod, int dim)              { return dmesh.coords[nod][dim]; }
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
   void   setPartitioningData          (PLAS_PART_DATA *part)          { part->numNodPairs = 0; }
   void   screenOutput                 (char *text)                    { printf("INFO: %s",text); }
   void   screenWarning                (char *text)                    { printf("WARN: %s",text); }

 private:  // member variables
  DRIVER_PARAMETERS  dparam;
  DRIVER_GAMBIT_MESH dmesh;
  DRIVER_FLOW_FIELD  dflow;
};


#endif
