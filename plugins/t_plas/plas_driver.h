
/*
 *  This file contains definitions of internal PLaS driver data
 *  structures and declares all internal functions.
 */


#ifndef PLAS_DRIVER_H
#define PLAS_DRIVER_H


// Included headers
#include <cstdlib>
#include <cstdio>
#ifdef MPI
#include <mpi.h>
#endif

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


// Function declarations
double plasdriver_CalcAreaTriangle(double c[3][2]);
void   plasdriver_CalcElementNeighbours(DRIVER_GAMBIT_MESH *dmesh);
void   plasdriver_CalcElementNormals(DRIVER_GAMBIT_MESH *dmesh);
void   plasdriver_CalcElementVolumes(DRIVER_GAMBIT_MESH *dmesh);
void   plasdriver_CalcElmsAroundNode(DRIVER_GAMBIT_MESH *dmesh);
void   plasdriver_CalcNodalVolumes(DRIVER_GAMBIT_MESH *dmesh);
double plasdriver_CalcVolumeTetra(double c[4][3]);
void   plasdriver_FreeGambitMemory(DRIVER_GAMBIT_MESH *dmesh);
void   plasdriver_GetFaceNodes(DRIVER_GAMBIT_MESH *dmesh, int elm, int face, int *faceNodes);
void   plasdriver_InitFlowField(DRIVER_GAMBIT_MESH *dmesh, int material, DRIVER_FLOW_FIELD *dflow);
void   plasdriver_ReadDriverDataFile(DRIVER_PARAMETERS *dparam);
void   plasdriver_ReadGambitNeutralFile(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam);
void   plasdriver_WriteTecplot(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam, DRIVER_FLOW_FIELD *dflow);


#endif  // PLAS_DRIVER_H

