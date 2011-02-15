
/*
 * This file contains definitions of the PLaS data structure
 * and interface functions to the flow solver. These functions
 * have to be implemented on the flow solver side and provide
 * data to PLaS.
 */

#ifndef PLAS_PLAS_H
#define PLAS_PLAS_H


/*
 * Preprocessor constants
 */

#define DFLAG_DISABLED       0
#define DFLAG_ENABLED        1
#define DFLAG_CREATED        2
#define DFLAG_LEFT           3
#define DFLAG_PASS           4
#define FLOW_PARTIC          1
#define FLOW_DROPLET         2
#define FLOW_BUBBLY          3
#define MAT_COPPER           1
#define MAT_POLY             2
#define MAT_WATER            3
#define MAT_NHEPTANE         4
#define MAT_HYDROGEN         5
#define MAT_OXYGEN           6
#define MAT_AIR              7
#define MAXELMFACES          6
#define MAXELMNODES          8
#define MAXELMNORMS          6
#define FORCE_PIC            1
#define FORCE_PROJ           2
#define COLL_UNCORRELATED    1
#define COLL_CORRELATED      2
#define ELM_SIMPLEX          1
#define ELM_PRISM            2
#define ELM_QUAD             3
#define ELM_HEX              4
#define ELM_PYRAMID          5
#define FLOWSOLVER_DRIVER    1
#define FLOWSOLVER_MORPHEUS  2
#define FLOWSOLVER_SFELES    3
#define FLOWSOLVER_COOLFLUID 4
#define SEARCH_ONEPROC       1
#define SEARCH_PARALLEL      2
#ifndef PI
#define PI 3.14159265
#endif
#ifndef Ru
#define Ru 8314.472
#endif


/*
 * Data of a dispersed entity (particle, droplet or bubble)
 */

typedef struct _plas_entity_data{
  int flag;             // Entity enabled/disabled flag
  int element;          // Grid element containing the entity
  int node;             // Nearest grid node to the entity
  double diameter;      // Entity diameter
  double *position;     // Entity position coordiantes
  double *velocity;     // Entity velocity components
  double *velocityOld;  // Entity velocity components of time step n-1
  double temperature;   // Entity temperature;
} PLAS_ENTITY_DATA;


/*
 * Data of the dispersed phase (per node)
 */

typedef struct _plas_phase_data{
  int numDens;         // Number density of entities (per node)
  double volFrac;      // Volume fraction (per node)
  double volFracDt;    // Volume fraction change rate (per node)
  double *volFracDx;   // Volume fraction gradient (per node)
  double avgDiam;      // Mean diameter (per node)
  double stdDiam;      // Diameter standard deviation (per node)
  double *avgVel;      // Mean entity velocity (per node)
  double *stdVel;      // Entity velocity standard deviation (per node)
  double *dispForce;   // Momentum force of the entities (per node)
  double avgRespTime;  // Average response time (per node)
} PLAS_PHASE_DATA;


/*
 * Data of the dispersed phase material
 */

typedef struct _plas_material_data{
  double rhoDisp;           // Density of dispersed entities
  double muDisp;            // Viscosity of dispersed entities
  double cpDisp;            // Specific heat capacity of the dispersed entities
  double kDisp;             // Thermal conductivity of the dispersed entities
  double sigDisp;           // Surface tension of the dispersed entities
  double epsDisp;           // Emissivity of the dispersed entities
  double satPresDisp;       // Saturation pressure of the dispersed entities evaluated at surface of entities
  double vapPres;           // Vapour pressure of the dispersed entities
  double latHeatDisp;       // Specific latent heat of the dispersed entities
  double molarMassDisp;     // Molar mass of the dispersed entities components
  double molarMassDispVap;  // Molar mass of the dispersed entities components at vapour phase
  double binaryDiffCoeff;   // Binary diffusion coefficient
  double massDiffCoeff;     // Mass diffusivity for a gas in a liquid
  double HeDisp;            // Henry Law constant for the dispersed entities
} PLAS_MATERIAL_DATA;


/*
 * Statistics data of PLaS
 */

typedef struct _plas_stats{
  int enabled;         // Number of active entities
  int in;              // Number of entities added
  int out;             // Number of entities removed
  int bounce;          // Number of wall bounces
  int periodic;        // Number of entities passing periodic boundaries
  int passed;          // Number of entities passing process borders
  int lost;            // Number of lost entities (should be zero)
  int coll;            // Number of particle or bubble collisions
  int coalesc;         // Number of coalescences of bubbles
  int leftproc;        // Internal counter for particles leaving a process
  double dtLagrAvg;    // Average lagrangian time scale
  double reynoldsAvg;  // Average entity Reynolds number
  double nusseltAvg;   // Average entity Nusselt number
  double subIterAvg;   // Average number of Lagrangian sub-iterations
} PLAS_STATS;


/*
 * PLaS input parameters (read from PLaS input file)
 */

typedef struct _plas_input_param{
  int numMaxEnt;               // Maximum number of entities per process
  int numIniEnt;               // Number of initially distributed entities:
  int numProdDom;              // Number of entity production zones
  int *prodDom;                // Geometrical shape of entity inlets
  double **prodParam;          // Coordinates of entity inlets
  double *massFluxes;          // Secondary phase mass flux of entity inlets
  int iniDiamType;             // Type of diameter distribution
  double iniDiam;              // Diameter of dispersed entities
  double iniDiamStd;           // Diameter standard deviation
  double iniVel[3];            // Initial velocity of dispersed entities
  double iniTempDisp;          // Temperature of dispersed entities
  int material;                // Flag: Entity material (defines flow type)
  int momentumCoupl;           // Flag: Momentum coupling
  int volfracCoupl;            // Flag: Volume fraction coupling
  int energyCoupl;             // Flag: Energy coupling
  int collisionModel;          // Flag: Collision model
  int liftForce;               // Flag: Lift force (only for bubbly flow)
  int evapModel;               // Flag: Evaporation model (only for droplet flow)
  int saturModel;              // Flag: Saturation model (only for bubbly flow)
  int perBnd;                  // Flag: Periodic boundaries for dispersed entities
  double gravVec[3];           // Gravity vector
  int restart;                 // Restart flag
  char* writeStatsFilename;    // Write output statistics filename
  char* writeTecplotFilename;  // Write output tecplot filename
  char* confFilename;          // Configuration filename
} PLAS_INPUT_PARAM;


/*
 * Partitioning data
 */

typedef struct _plas_part_data{
  int numNodPairs;   // Number of updatable/ghost node pairs
  int *sendRank;     // Sending ranks
  int *sendNodeIdx;  // Sending rank local node numbers
  int *recvRank;     // Receiving ranks
  int *recvNodeIdx;  // Receiving rank local node numbers
} PLAS_PART_DATA;


/*
 * Data to be set by the flow solver
 */

typedef struct _plas_flowsolver_param{

  //***To be set on initialization (see interface function below)***//

  int flowSolver;       // Flag for the flow solver used
  int numDim;           // Number of space dimensions
  int numUnk;           // Number of unknown variables for the primary phase flow
  int numNod;           // Number of nodes
  int numElm;           // Number of elements
  int numBnd;           // Number of boundaries
  double rhoCont;       // Primary phase density
  double muCont;        // Primary phase viscosity
  double nuCont;        // Primary phase kinematic viscosity
  double cpCont;        // Specific heat coefficient of the flow medium
  double kCont;         // Thermal conductivity of the flow medium
  double dtEul;         // Eulerian time scale
  double domainVolume;  // Total volume (in 2D: area) of the domain
  double minElmVolume;  // Smallest element volume (in 2D: area) in the domain
  double maxElmVolume;  // Largest element volume (in 2D: area) in the domain
  PLAS_PART_DATA part;  // Partitioning data structure

  //***To be set on every iteration (see interface function below)***//

  int iter;             // Current iteration
  double time;          // Current time
  int writeOutput;      // Flag whether the flow solver writes outputs

} PLAS_FLOWSOLVER_PARAM;


/*
 * Internal PLaS runtime variables
 */

typedef struct _plas_runtime_param{
  double lagrTimeFactor;  // Lagrangian time factor (constant)
  double errTol;          // Error tolerance
  int flowType;           // Flag: Type of flow (particle, droplet, bubble)
  double *massResid;      // Mass flux residual
  double wallElasticity;  // Elasticity factor for wall bounces
} PLAS_RUNTIME_PARAM;


/*
 * PLaS data structure (instantiate in the flow solver)
 */

typedef struct _plas_data{
  PLAS_ENTITY_DATA *ed;      // Entity data structure (per entity)
  PLAS_PHASE_DATA *pd;       // Phase data structure (per node)
  PLAS_MATERIAL_DATA md;     // Material data structure
  PLAS_STATS sd;             // Statistics data structure
  PLAS_INPUT_PARAM ip;       // Input file parameter structure
  PLAS_FLOWSOLVER_PARAM fp;  // Flowsolver parameter structure
  PLAS_RUNTIME_PARAM rp;     // PLaS internal runtime parameter structure
  int numExtEnt;             // Number of bubbles coming from external code
  double *extEntPos;         // Positions of bubbles coming from external code
  double *extEntVel;         // Velocities of bubbles coming from external code
  double *extEntTemp;        // Temperatures of bubbles coming from external code
  double *extEntDiam;        // Diameters of bubbles coming from external code
} PLAS_DATA;


/*
 * Declaration of PLaS interface functions
 */

void initPLaS(PLAS_DATA *data);
void runPLaS(PLAS_DATA *data);
void terminatePLaS(PLAS_DATA *data);

void   plas_CalcCrossProduct_3D(double *value, double *a, double *b);       
double plas_CalcVectScalarProduct(int numDim, double *a, double *b);       
double plas_ReadDoubleParam(FILE *inpFile);       
int    plas_ReadIntParam(FILE *inpFile);       
void   plas_TerminateOnError(char *errMessage);       


////////////////////////////////////////////////////////////////////////////////////////
//  Declaration of PLaS interface functions to be implemented in calling flow solver  //
////////////////////////////////////////////////////////////////////////////////////////

// --------------------------------------------
// Set flow solver parameters on initialization
// --------------------------------------------
// fp = pointer to PLaS data structure
//
void plasinterface_setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);

// ---------------------------------------
// Set flow solver parameters on time step
// ---------------------------------------
// fp = pointer to PLaS data structure
//
void plasinterface_setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp);

//----------------------
// Set partitioning data
//----------------------
// part = partitioning data structure
//
void plasinterface_setPartitioningData(PLAS_PART_DATA *part);

// ----------------------------------
// Provide node of an element to PLaS
// ----------------------------------
// elm = element index
// enod = node of the element
// return value = node index
//
int plasinterface_getElmNode(int elm, int enod);

// ---------------------------------------
// Provide neighbour of an element to PLaS
// ---------------------------------------
// elm = element index
// eface = face of the element
// return value = neighbour element index
//
int plasinterface_getElmNeighbour(int elm, int eface);

// ----------------------------------------------------------
// Provide arbitrary reference point of boundary face to PLaS
// ----------------------------------------------------------
// bnd = boundary index
// bface = face of the boundary
// dim = coordinate index
// return value = reference point of the face
//
double plasinterface_getBndFaceRefCoord(int bnd, int bface, int dim);

// ----------------------------------------------------------
// Provide domain element associated to boundary face to PLaS
// ----------------------------------------------------------
// bnd = boundary index
// bface = face of the boundary
// dummy = not needed
// return value = element index
//
int plasinterface_getBndDomElm(int bnd, int bface, int dummy);

// -------------------------------
// Provide node coordinate to PLaS
// -------------------------------
// nod = node index
// dim = coordinate index
// return value = coordinate of the node
//
double plasinterface_getNodCoord(int nod, int dim);

// -------------------------------------------
// Provide component of element normal to PLaS
// -------------------------------------------
// elm = element index
// eface = face of the element
// dim = coordinate index
// return value = coordinate of the normal
//
double plasinterface_getElmNormComp(int elm, int eface, int dim);

// -------------------------------------------------
// Provide component of boundary face normal to PLaS
// -------------------------------------------------
// elm = element index
// bface = face of the boundary
// dim = coordinate index
// return value = coordinate of the normal
//
double plasinterface_getBndFaceNormComp(int bnd, int bface, int dim);

// -----------------------------------------------------
// Provide component of element face middle-point vector
// -----------------------------------------------------
// elm = element index
// eface = face of the element
// dim = coordinate index
// return value = coordinate of the face middle-point
//
double plasinterface_getElmFaceMiddlePoint(int elm, int eface, int dim);

// ---------------------------------
// Provide nodal area/volume to PLaS
// ---------------------------------
// nod = node index
// return value = nodal cell volume
//
double plasinterface_getNodVol(int nod);

// -----------------------------------
// Provide element area/volume to PLaS
// -----------------------------------
// elm = element index
// return value = element volume
//
double plasinterface_getElmVol(int elm);

// ---------------------------------------------
// Provide number of faces of a boundary to PLaS
// ---------------------------------------------
// bnd = boundary index
// return value = number of boundary faces
//
int plasinterface_getNumBndFaces(int bnd);

// ----------------------------------------------------------
// Provide information about which boundary is a wall to PLaS
// ----------------------------------------------------------
// bnd = boundary index
// return value = value if boundary is a wall, else zero
//
int plasinterface_getWallBndFlag(int bnd);

// ----------------------------------------------------
// Provide information about which boundary is periodic
// ----------------------------------------------------
// bnd = boundary index
// return value = value if boundary is periodic, else zero
//
int plasinterface_getPerBndFlag(int bnd);

// ----------------------------------------------------------
// Provide information about periodic boundary offset to PLaS
// ----------------------------------------------------------
// bnd = boundary index
// dim = coordinate index
// return value = periodic boundary offset
//
double plasinterface_getPerBndOffset(int bnd, int dim);

// ----------------------------------
// Provide type of an element to PLaS
// ----------------------------------
// elm = element index
// return value = element type (codes in #defines)
//
int plasinterface_getElementType(int elm);

// -------------------------------------------------------
// Provide nodal velocity component at time step n to PLaS
// -------------------------------------------------------
// nod = node index
// dim = coordinate index
// return value = velocity at time step n
//
double plasinterface_getVelocityComp(int nod, int dim);

// ---------------------------------------------------------
// Provide nodal velocity component at time step n-1 to PLaS
// ---------------------------------------------------------
// nod = node index
// dim = coordinate index
// return value = velocity at time step n-1
//
double plasinterface_getVelocityCompOld(int nod, int dim);

// ----------------------------------------------------
// Provide velocity component derivative at time step n
// ----------------------------------------------------
// nod = node index
// idim = variable to take derivatrive of
// jdim = coordinate direction of derivative
// return value = velocity derivative at time step n-1
//
double plasinterface_getVelocityDerivativeComp(int nod, int idim, int jdim);

// ------------------------------------------------------
// Provide velocity component derivative at time step n-1
// ------------------------------------------------------
// nod = node index
// idim = variable to take derivatrive of
// jdim = coordinate direction of derivative
// return value = velocity derivative at time step n-1
//
double plasinterface_getVelocityDerivativeCompOld(int nod, int idim, int jdim);

//-------------------------------------------------
// Provide nodal temperature at time step n to PLaS
//-------------------------------------------------
// nod = node idx
// return value = temperature at time step n
//
double plasinterface_getTemperature(int nod);

//---------------------------------------------------
// Provide nodal temperature at time step n-1 to PLaS
//---------------------------------------------------
// nod = node idx
// return value = temperature at time step n-1
//
double plasinterface_getTemperatureOld(int nod);

//----------------------------------------------
// Provide nodal pressure at time step n to PLaS
//----------------------------------------------
// nod = node idx
// return value = pressure at time step n
//
double plasinterface_getPressure(int nod);

//------------------------------------------------
// Provide nodal pressure at time step n-1 to PLaS
//------------------------------------------------
// nod = node idx
// return value = pressure at time step n-1
//
double plasinterface_getPressureOld(int nod);

// -----------------------------------------------------
// Provide starting element for local brute force search
// -----------------------------------------------------
// pos = position vector (xyz)
// return value = lowest element for BF search (0)
//
int plasinterface_StartElementSearch(double *pos);

// ---------------------------------------------------
// Provide ending element for local brute force search
// ---------------------------------------------------
// pos = position vector (xyz)
// return value = highest element for BF search (numElm-1)
//
int plasinterface_EndElementSearch(double *pos);

// -------------------------------------------------------
// Provide Eulerian time scale (tke/disspation) for a node
// -------------------------------------------------------
// nod = node index
// return value = Eulerian time scale (k/eps)
//
double plasinterface_getEulerianTimeScale(int nod);

// ---------------------------
// Pass data for screen output
// ---------------------------
// text = screen output
//
void plasinterface_screenOutput(char *text);

// ------------------------------
// Pass data for a screen warning
// ------------------------------
// text = screen warning
//
void plasinterface_screenWarning(char *text);

#endif  // PLAS_PLAS_H

