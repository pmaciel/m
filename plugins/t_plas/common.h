
/*
 * This file contains definitions of internal PLaS data
 * structures and declares all internal functions.
 */

#ifndef PLAS_COMMON_H
#define PLAS_COMMON_H


// Included headers
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#ifdef MPI
#include <mpi.h>
#endif

#include "plas.h"


// Entity element data structures
typedef struct _entity_element_data{
  int numElmNodes;          // Number of element nodes
  int *elmNodes;            // Element nodes
  int numElmFaces;          // Number of element faces
  double **elmFaceVectors;  // Element face middle points
  double **elmNorms;        // Element face normals
} ENTITY_ELEMENT_DATA;


// Local entity variables data structure
typedef struct _local_entity_variables{
  int flag;                   // Flag to determine if the entity is active
  int node;                   // Number of corresponding grid node
  int elm;                    // Number of corresponding grid element
  double diam;                // Diameter
  double *pos;                // Entity position at time step n
  double *posOld;             // Entity position at time step n-1
  double *vel;                // Entity velocity at time step n
  double *velOld;             // Entity velocity at time step n-1
  double *relVel;             // Relative velocity
  double normVel;             // Norm of the relative velocity
  double temp;                // Temperature
  double relTemp;             // Relative temperature
  double reynolds;            // Dispersed Reynold number
  double nusselt;             // Nusselt number
  double sherwood;            // Sherwood number
  double schmidt;             // Schmidt number
  double spalding;            // Spalding number
  double prandtl;             // Prandtl number
  double dragCoeff;           // Drag coefficient
  double liftCoeff;           // Lift coefficient
  double kinRespTime;         // Kinematic response time
  double thermRespTime;       // Thermal response time
  double massTrCoeff;         // Mass transfer coefficient beta
  double pressBubble;         // Bubble interal pressure
  double rhoBubble;           // Bubble internal density
  double concInterf;          // Bubble surface concentration
  double concDelta;           // Bubble concentration boundary layer
  ENTITY_ELEMENT_DATA edata;  // Corresponding grid element data
} LOCAL_ENTITY_VARIABLES;


// Local flow variables data structure
typedef struct _local_flow_variables{
  double *vel;           // Local instantaneous flow velocity
  double *velDt;         // Flow velocity time derivative
  double **velDx;        // Flow velocity space derivative
  double *vort;          // Local instantaneous flow vorticity
  double pressure;       // Local instantaneous flow pressure
  //double concentration;  // Local instanteneous flow (solvent) concentration
  double temp;           // Local instantaneous flow temperature
} LOCAL_FLOW_VARIABLES;


// Function declarations
void   plas_AllocateLocalEntityVar(int numDim, LOCAL_ENTITY_VARIABLES *ent);       
void   plas_AllocateLocalFlowVar(int numDim, LOCAL_FLOW_VARIABLES *flow);       
void   plas_BroadcastParameters(PLAS_DATA *data);       
void   plas_CalcBackCoupling(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double *force, double tFactor);       
void   plas_CalcBoundaryUnitNormal(int numDim, int ibnd, int ifac, double *normVec);       
void   plas_CalcCellwiseData(PLAS_DATA *data);       
double plas_CalcConcInterf(PLAS_DATA *data, double pressBubble);       
void   plas_CalcCouplingForcesBubble(PLAS_DATA *data,LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor);       
void   plas_CalcCouplingForcesParticle(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor);       
double plas_CalcDispReynolds(double viscosity, double diameter, double normVel);       
double plas_CalcDragCoeff(int flowType, double reynolds);       
void   plas_CalcEntityCoefficients(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow);       
double plas_CalcKinematicResponseTime(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent);       
double plas_CalcLiftCoeff(int flowType);       
double plas_CalcMassTransferCoeff(PLAS_DATA *data, double sherwood, double spalding);       
void   plas_CalcMatVectScalarProduct_2D(double *value, double **m, double *a);       
void   plas_CalcMatVectScalarProduct_3D(double *value, double **m, double *a);       
void   plas_CalcMaterialData(PLAS_DATA *data, double T, double p);       
void   plas_CalcNodeImpactFactors(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, double *dist);       
double plas_CalcNusseltNumber(int evapModel, double reynolds, double spalding, double prandtl);       
double plas_CalcPrandtlNumber(PLAS_DATA *data);       
double plas_CalcPressBubble(PLAS_DATA *data, double diameter, double pressure);       
double plas_CalcRhoBubble(PLAS_DATA *data, double temperature, double pressBubble);       
void   plas_CalcRotationMatrix_2D(double phi, double **m);       
void   plas_CalcRotationMatrix_3D(double phi, double **m, int axis);       
double plas_CalcSchmidtNumber(PLAS_DATA *data);       
double plas_CalcSherwoodNumber(int evapModel, double reynolds, double schmidt, double spalding);       
double plas_CalcSpaldingNumber(PLAS_DATA *data, double pressure);       
double plas_CalcThermalResponseTime(PLAS_DATA *data, double diameter);       
void   plas_CalcTrajectory(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double dtLagr);       
double plas_CalcVectorAngle(int numDim, double *a, double *b);       
double plas_CalcVectorLength(int numDim, double *a);       
void   plas_CalcVectorRotation_3D(double psi, double **m, double *a);       
void   plas_CalcVorticity(int numDim, LOCAL_FLOW_VARIABLES *flow);       
double plas_CalcWallFaceDistance(int numDim, double *pos, int ibnd, int ifac);       
void   plas_CheckNaN(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent);       
int    plas_Coalescence(PLAS_DATA *data, double dj, double di, double *uRelPrPr, double *x, double *y, double *z, double Mi, double Mj, double *uiPrPr, double *uiPrPrNew, double *uiNew, double *ui);       
void   plas_CollisionModel(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, int numDens, double dtLagr);       
void   plas_CreateStatsFile(PLAS_DATA *data, char *outpString);       
void   plas_CreateTecplotFile(PLAS_DATA *data, char *outpString);       
void   plas_DeallocateLocalEntityVar(LOCAL_ENTITY_VARIABLES *ent);       
void   plas_DeallocateLocalFlowVar(int numDim, LOCAL_FLOW_VARIABLES *flow);       
void   plas_FindExitFace(int numBnd, int numDim, LOCAL_ENTITY_VARIABLES *ent, int *f, int *i, int *j);       
void   plas_FindMinimumElementFaceDistance(int numDim, LOCAL_ENTITY_VARIABLES *ent, int *idx, double *dmin);       
int    plas_FindNearestElementNode(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent);       
void   plas_ImposeExternal(PLAS_DATA *data);       
void   plas_ImposeProductionDomains(PLAS_DATA *data);       
void   plas_Interpolate(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step);       
void   plas_InterpolatePressure(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step);       
void   plas_InterpolateTemperature(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step);       
void   plas_InterpolateVelocity(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step);       
void   plas_LoadInitialDistribution(PLAS_DATA *data, char *inpString);       
int    plas_MpiAllMaxInt(int val);       
void   plas_MpiAllMaxIntArray(int *val, int size);       
int    plas_MpiAllMinInt(int val);       
double plas_MpiAllSumDouble(double val);       
void   plas_MpiAllSumDoubleArray(double *val, int size);       
int    plas_MpiAllSumInt(int val);       
void   plas_MpiAllSumIntArray(int *val, int size);       
void   plas_MpiBarrier();       
void   plas_MpiBroadcastDouble(double *variable, int size, int root);       
void   plas_MpiBroadcastInt(int *variable, int size, int root);       
int    plas_MpiGetNumProc();       
int    plas_MpiGetRank();       
void   plas_NormalizeVector(int numDim, double *a);       
void   plas_PassEntities(PLAS_DATA *data);       
double plas_RandomDouble();       
double plas_RandomGaussian(float m, float s);       
void   plas_RandomInitialDistribution(PLAS_DATA *data);       
int    plas_RandomInteger(int min, int max);       
void   plas_ReadParameters(PLAS_DATA *data);       
void   plas_SearchDomainParallel(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent);       
void   plas_SearchSuccessive(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent);       
double plas_SetDiameter(PLAS_DATA *data);       
void   plas_SetElementFaces(int numDim, LOCAL_ENTITY_VARIABLES *ent);       
void   plas_SetElementGeometry(int numDim, LOCAL_ENTITY_VARIABLES *ent);       
void   plas_SetElementNodes(int numDim, LOCAL_ENTITY_VARIABLES *ent);       
void   plas_SolveGaussSeidel(int numDim, double **mat, double *s, double *rhs);       
void   plas_WallBounce(int numDim, double elasticity, LOCAL_ENTITY_VARIABLES *ent, int ibnd, int ifac);       
void   plas_WriteStatsFile(PLAS_DATA *data, char *outpString, int iter, double time);       
void   plas_WriteTecplotFile(PLAS_DATA *data, char *outpString, int iter, double time);       


#endif  // PLAS_COMMON_H       

