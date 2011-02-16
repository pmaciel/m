
#include "t_plas.h"

/*
 * This file includes the interface between the PLaS driver
 * program and PLaS itself. All functions below are called
 * during PLaS runtime. The driver program has to provide data
 * about the mesh and the flow field to PLaS.
 *
 * Documentation of input/output parameters in plas.h
 */


// Set flow solver parameters on initialization.
void plasinterface_setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->flowSolver = FLOWSOLVER_DRIVER;
  fp->numDim = p_dmesh->numDim;
  fp->numUnk = p_dparam->numUnk;
  fp->numNod = p_dmesh->numNod;
  fp->numElm = p_dmesh->numElm;
  fp->numBnd = p_dmesh->numBnd;
  fp->rhoCont = p_dflow->rho;
  fp->muCont = p_dflow->mu;
  fp->nuCont = p_dflow->mu/p_dflow->rho;
  fp->cpCont = p_dflow->cp;
  fp->kCont = p_dflow->k;
  fp->dtEul = p_dparam->dtEul;
  fp->domainVolume = p_dmesh->domainVolume;
  fp->minElmVolume = p_dmesh->minElmVolume;
  fp->maxElmVolume = p_dmesh->maxElmVolume;

  fp->time = 0.;
  fp->writeOutput = 1;
}

// Set flow solver parameters on time step.
void plasinterface_setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time += p_dparam->dtEul;
  fp->iter = p_dparam->iter;
}

// Set partitioning data.
void plasinterface_setPartitioningData(PLAS_PART_DATA *part)
{
  part->numNodPairs = 0;
}

// Provide node of an element to PLaS.
int plasinterface_getElmNode(int elm, int enod)
{
  return p_dmesh->elmNodes[elm][enod];
}

// Provide neighbour of an element to PLaS.
int plasinterface_getElmNeighbour(int elm, int eface)
{
  return p_dmesh->elmNeighbs[elm][eface];
}

// Provide arbitrary reference point of boundary face to PLaS.
double plasinterface_getBndFaceRefCoord(int bnd, int bface, int dim)
{
  int faceNodes[4];

  ::pt_plas->plasdriver_GetFaceNodes(p_dmesh,p_dmesh->bndDomElms[bnd][bface],p_dmesh->bndFaces[bnd][bface],faceNodes);

  return p_dmesh->coords[faceNodes[0]][dim];
}

// Provide domain element associated to boundary face to PLaS.
int plasinterface_getBndDomElm(int bnd, int bface, int dummy)
{
  return p_dmesh->bndDomElms[bnd][bface];
}

// Provide node coordinate to PLaS.
double plasinterface_getNodCoord(int nod, int dim)
{
  return p_dmesh->coords[nod][dim];
}

// Provide component of element normal to PLaS.
double plasinterface_getElmNormComp(int elm, int eface, int dim)
{
  return p_dmesh->elmNorms[elm][eface][dim];
}

// Provide component of boundary face normal to PLaS.
double plasinterface_getBndFaceNormComp(int bnd, int face, int dim)
{
  return p_dmesh->elmNorms[p_dmesh->bndDomElms[bnd][face]][p_dmesh->bndFaces[bnd][face]][dim];
}

// Provide component of element face middle-point vector.
double plasinterface_getElmFaceMiddlePoint(int elm, int eface, int dim)
{
  int ifac,faceNodes[4];

  int ctr = 0;
  double coord = 0.0;

  ::pt_plas->plasdriver_GetFaceNodes(p_dmesh,elm,eface,faceNodes);

  for(ifac=0; ifac<4; ifac++){
    if(faceNodes[ifac]!=-1){
      coord += p_dmesh->coords[faceNodes[ifac]][dim];
      ctr++;
    }
  }
  coord /= ctr;

  return coord;
}

// Provide nodal area/volume to PLaS.
double plasinterface_getNodVol(int nod)
{
  return p_dmesh->nodVolumes[nod];
}

// Provide element area/volume to PLaS.
double plasinterface_getElmVol(int elm)
{
  return p_dmesh->elmVolumes[elm];
}

// Provide number of faces of a boundary to PLaS.
int plasinterface_getNumBndFaces(int bnd)
{
  return p_dmesh->numBndFaces[bnd];
}

// Provide information about which boundary is a wall to PLaS.
int plasinterface_getWallBndFlag(int bnd)
{
  int flag;

  if(p_dmesh->bndTypes[bnd]==1){
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

  if(p_dmesh->bndTypes[bnd]<0){
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
  return p_dmesh->elmTypes[elm];
}

// Provide nodal velocity component at time step n to PLaS.
double plasinterface_getVelocityComp(int nod, int dim)
{
  return p_dflow->u[nod][dim];
}

// Provide nodal velocity component at time step n-1 to PLaS.
double plasinterface_getVelocityCompOld(int nod, int dim)
{
  return p_dflow->u[nod][dim];
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
  return p_dflow->T[nod];
}

// Provide nodal temperature at time step n-1 to PLaS.
double plasinterface_getTemperatureOld(int nod)
{
  return p_dflow->T[nod];
}

// Provide nodal pressure at time step n to PLaS.
double plasinterface_getPressure(int nod)
{
  return p_dflow->p[nod];
}

// Provide nodal pressure at time step n-1 to PLaS.
double plasinterface_getPressureOld(int nod)
{
  return p_dflow->p[nod];
}

// Provide starting element for local brute force search.
int plasinterface_StartElementSearch(double *pos)
{
  return 0;
}

// Provide ending element for local brute force search.
int plasinterface_EndElementSearch(double *pos)
{
  return p_dmesh->numElm-1;
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

