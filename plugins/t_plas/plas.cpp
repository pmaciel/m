
#include "plas_driver.h"

/*
 * This file includes the interface between the PLaS driver
 * program and PLaS itself. All functions below are called
 * during PLaS runtime. The driver program has to provide data
 * about the mesh and the flow field to PLaS.
 *
 * Documentation of input/output parameters in plas.h
 */

extern DRIVER_PARAMETERS dparam;
extern DRIVER_GAMBIT_MESH dmesh;
extern DRIVER_FLOW_FIELD dflow;

// Set flow solver parameters on initialization.
void plasinterface_setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->flowSolver = FLOWSOLVER_DRIVER;
  fp->numDim = dmesh.numDim;
  fp->numUnk = dparam.numUnk;
  fp->numNod = dmesh.numNod;
  fp->numElm = dmesh.numElm;
  fp->numBnd = dmesh.numBnd;
  fp->rhoCont = dflow.rho;
  fp->muCont = dflow.mu;
  fp->nuCont = dflow.mu/dflow.rho;
  fp->cpCont = dflow.cp;
  fp->kCont = dflow.k;
  fp->dtEul = dparam.dtEul;
  fp->domainVolume = dmesh.domainVolume;
  fp->minElmVolume = dmesh.minElmVolume;
  fp->maxElmVolume = dmesh.maxElmVolume;

  fp->time = 0.;
  fp->writeOutput = 1;
}

// Set flow solver parameters on time step.
void plasinterface_setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time += dparam.dtEul;
  fp->iter = dparam.iter;
}

// Set partitioning data.
void plasinterface_setPartitioningData(PLAS_PART_DATA *part)
{
  part->numNodPairs = 0;
}

// Provide node of an element to PLaS.
int plasinterface_getElmNode(int elm, int enod)
{
  return dmesh.elmNodes[elm][enod];
}

// Provide neighbour of an element to PLaS.
int plasinterface_getElmNeighbour(int elm, int eface)
{
  return dmesh.elmNeighbs[elm][eface];
}

// Provide arbitrary reference point of boundary face to PLaS.
double plasinterface_getBndFaceRefCoord(int bnd, int bface, int dim)
{
  int faceNodes[4];

  plasdriver_GetFaceNodes(&dmesh,dmesh.bndDomElms[bnd][bface],dmesh.bndFaces[bnd][bface],faceNodes);

  return dmesh.coords[faceNodes[0]][dim];
}

// Provide domain element associated to boundary face to PLaS.
int plasinterface_getBndDomElm(int bnd, int bface, int dummy)
{
  return dmesh.bndDomElms[bnd][bface];
}

// Provide node coordinate to PLaS.
double plasinterface_getNodCoord(int nod, int dim)
{
  return dmesh.coords[nod][dim];
}

// Provide component of element normal to PLaS.
double plasinterface_getElmNormComp(int elm, int eface, int dim)
{
  return dmesh.elmNorms[elm][eface][dim];
}

// Provide component of boundary face normal to PLaS.
double plasinterface_getBndFaceNormComp(int bnd, int face, int dim)
{
  return dmesh.elmNorms[dmesh.bndDomElms[bnd][face]][dmesh.bndFaces[bnd][face]][dim];
}

// Provide component of element face middle-point vector.
double plasinterface_getElmFaceMiddlePoint(int elm, int eface, int dim)
{
  int ifac,faceNodes[4];

  int ctr = 0;
  double coord = 0.0;

  plasdriver_GetFaceNodes(&dmesh,elm,eface,faceNodes);

  for(ifac=0; ifac<4; ifac++){
    if(faceNodes[ifac]!=-1){
      coord += dmesh.coords[faceNodes[ifac]][dim];
      ctr++;
    }
  }
  coord /= ctr;

  return coord;
}

// Provide nodal area/volume to PLaS.
double plasinterface_getNodVol(int nod)
{
  return dmesh.nodVolumes[nod];
}

// Provide element area/volume to PLaS.
double plasinterface_getElmVol(int elm)
{
  return dmesh.elmVolumes[elm];
}

// Provide number of faces of a boundary to PLaS.
int plasinterface_getNumBndFaces(int bnd)
{
  return dmesh.numBndFaces[bnd];
}

// Provide information about which boundary is a wall to PLaS.
int plasinterface_getWallBndFlag(int bnd)
{
  int flag;

  if(dmesh.bndTypes[bnd]==1){
    flag = 1;
  } else{
    flag = 0;
  }

  return flag;
}

// Provide information about which boundary is periodic.
int plasinterface_getPerBndFlag(int bnd)
{
  int flag;

  if(dmesh.bndTypes[bnd]<0){
    flag = 1;
  } else{
    flag = 0;
  }

  return flag;
}

// Provide information about periodic boundary offset to PLaS.
double plasinterface_getPerBndOffset(int bnd, int dim)
{
  //------------------
  // To be implemented
  //------------------

  return 0.0;
}

// Provide type of an element to PLaS.
int plasinterface_getElementType(int elm)
{
  return dmesh.elmTypes[elm];
}

// Provide nodal velocity component at time step n to PLaS.
double plasinterface_getVelocityComp(int nod, int dim)
{
  return dflow.u[nod][dim];
}

// Provide nodal velocity component at time step n-1 to PLaS.
double plasinterface_getVelocityCompOld(int nod, int dim)
{
  return dflow.u[nod][dim];
}

// Provide velocity component derivative at time step n.
double plasinterface_getVelocityDerivativeComp(int nod, int idim, int jdim)
{
  return 0.0;
}

// Provide velocity component derivative at time step n-1.
double plasinterface_getVelocityDerivativeCompOld(int nod, int idim, int jdim)
{
  return 0.0;
}

// Provide nodal temperature at time step n to PLaS.
double plasinterface_getTemperature(int nod)
{
  return dflow.T[nod];
}

// Provide nodal temperature at time step n-1 to PLaS.
double plasinterface_getTemperatureOld(int nod)
{
  return dflow.T[nod];
}

// Provide nodal pressure at time step n to PLaS.
double plasinterface_getPressure(int nod)
{
  return dflow.p[nod];
}

// Provide nodal pressure at time step n-1 to PLaS.
double plasinterface_getPressureOld(int nod)
{
  return dflow.p[nod];
}

// Provide starting element for local brute force search.
int plasinterface_StartElementSearch(double *pos)
{
  return 0;
}

// Provide ending element for local brute force search.
int plasinterface_EndElementSearch(double *pos)
{
  return dmesh.numElm-1;
}

// Provide Eulerian time scale (tke/disspation) for a node.
double plasinterface_getEulerianTimeScale(int nod)
{
  return 0.0;
}

// Pass data for screen output.
void plasinterface_screenOutput(char *text)
{
  printf("   %s",text);
}

// Pass data for a screen warning.
void plasinterface_screenWarning(char *text)
{
  printf("WARNING: %s",text);
}

