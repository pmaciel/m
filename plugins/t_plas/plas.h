#ifndef PLAS_PLAS_H
#define PLAS_PLAS_H


#include <string>
#include <vector>
#include "ext/xmlParser.h"


/// Preprocessor constants
#define DFLAG_DISABLED    0
#define DFLAG_ENABLED     1
#define DFLAG_CREATED     2
#define DFLAG_LEFT        3

#define MAXELMFACES       6
#define MAXELMNODES       8
#define FORCE_PIC         1
#define FORCE_PROJ        2
#define COLL_UNCORRELATED 1
#define COLL_CORRELATED   2

#ifndef PI
#define PI 3.14159265
#endif
#ifndef Ru
#define Ru 8314.472
#endif


/// type of flow associated to the dispersed phase material
enum plas_flowtype_t { FLOW_PARTIC, FLOW_DROPLET, FLOW_BUBBLY };


/// supported element types (following Gambit convention numbering)
enum plas_elmtype_t { ELM_UNDEFINED, ELM_TRIANGLE, ELM_TETRAHEDRON, ELM_WEDGE, ELM_QUAD, ELM_BRICK, ELM_PYRAMID, ELM_EDGE, ALL_ELEMENTS };


/// Available quantities (coordinates, pressure, temperature, velocity and its spacial derivatives)
enum plas_quantity_t {
  QUANTITY_UNDEFINED,
  COORD,                            // vector
  PRESSURE       = COORD+3,         // scalar
  TEMPERATURE,                      // scalar
  VELOCITY,                         // vector
  VELOCITY_X_D   = VELOCITY+3,      // vector
  VELOCITY_Y_D   = VELOCITY_X_D+3,  // vector
  VELOCITY_Z_D   = VELOCITY_Y_D+3,  // vector
  ALL_QUANTITIES = VELOCITY_Z_D+3
};


/// Data of a dispersed entity (particle, droplet or bubble)
struct PLAS_ENTITY_DATA
{
  int flag;             // Entity enabled/disabled flag
  int element;          // Grid element containing the entity
  int zone;             // Element zone index
  int node;             // Nearest grid node to the entity
  double diameter;      // Entity diameter
  double *position;     // Entity position coordiantes
  double *velocity;     // Entity velocity components
  double *velocityOld;  // Entity velocity components of time step n-1
  double temperature;   // Entity temperature;
};


/// Data of the dispersed phase (per node)
struct PLAS_PHASE_DATA
{
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
};


/// Statistics data of PLaS
struct PLAS_STATS
{
  int enabled;         // Number of active entities
  int in;              // Number of entities added
  int out;             // Number of entities removed
  int bounce;          // Number of wall bounces
  int lost;            // Number of lost entities (should be zero)
  int coll;            // Number of particle or bubble collisions
  int coalesc;         // Number of coalescences of bubbles
  double dtLagrAvg;    // Average lagrangian time scale
  double reynoldsAvg;  // Average entity Reynolds number
  double nusseltAvg;   // Average entity Nusselt number
  double subIterAvg;   // Average number of Lagrangian sub-iterations
};


/// PLaS input parameters (read from PLaS input file)
struct PLAS_INPUT_PARAM
{
  int numMaxEnt;           // Maximum number of entities per process
  int numIniEnt;           // Number of initially distributed entities:
  int numProdDom;          // Number of entity production zones
  int *prodDom;            // Geometrical shape of entity inlets
  double **prodParam;      // Coordinates of entity inlets
  double *massFluxes;      // Secondary phase mass flux of entity inlets
  int iniDiamType;         // Type of diameter distribution
  double iniDiam;          // Diameter of dispersed entities
  double iniDiamStd;       // Diameter standard deviation
  double iniVel[3];        // Initial velocity of dispersed entities
  double iniTempDisp;      // Temperature of dispersed entities
  int momentumCoupl;       // Flag: Momentum coupling
  int volfracCoupl;        // Flag: Volume fraction coupling
  int energyCoupl;         // Flag: Energy coupling
  int collisionModel;      // Flag: Collision model
  int liftForce;           // Flag: Lift force (only for bubbly flow)
  int evapModel;           // Flag: Evaporation model (only for droplet flow)
  int saturModel;          // Flag: Saturation model (only for bubbly flow)
  double gravVec[3];       // Gravity vector

  double lagrTimeFactor;   // Lagrangian time factor (constant)
  double errTol;           // Error tolerance
  double wallElasticity;   // Elasticity factor for wall bounces

  std::string
    writeStatsFilename,    // Write output statistics filename
    writeTecplotFilename;  // Write output tecplot filename
};


/// Data to be set by the flow solver
struct PLAS_FLOWSOLVER_PARAM
{
  // to be set on initialization
  int
    numDim,          // Number of space dimensions
    numUnk,          // Number of unknown variables for the primary phase flow
    numNod;          // Number of nodes
  double dtEul;      // Eulerian time scale
  std::vector< unsigned >
    nInnerElements,  // number of inner elements, per zone
    nBoundElements;  // number of boundary elements, per zone

  // updated on calls to plas::run()
  int iter;          // Current iteration
  double time;       // Current time
};


/// Entity element data structures
struct ENTITY_ELEMENT_DATA
{
  ENTITY_ELEMENT_DATA(const int dim) :
    elmNodes(MAXELMNODES,-1),
    elmFaceVectors(MAXELMFACES,std::vector< double >(dim,0.)),
    elmNorms      (MAXELMFACES,std::vector< double >(dim,0.)) {}
  ~ENTITY_ELEMENT_DATA() {}
  int
    numElmNodes,                // Number of element nodes
    numElmFaces;                // Number of element faces
  std::vector< int > elmNodes;  // Element nodes
  std::vector< std::vector< double > >
    elmFaceVectors,             // Element face middle points
    elmNorms;                   // Element face normals
};


/// Local entity variables data structure
struct LOCAL_ENTITY_VARIABLES
{
  LOCAL_ENTITY_VARIABLES(const int dim) :
    elm(-1),
    zone(-1),
    node(-1),
    pos   (dim,0.),
    posOld(dim,0.),
    vel   (dim,0.),
    velOld(dim,0.),
    relVel(dim,0.),
    edata(dim) {}
  ~LOCAL_ENTITY_VARIABLES() {}
  int
    flag,           // Flag to determine if the entity is active
    elm,            // Number of corresponding grid element
    zone,           // Element zone index
    node;           // Number of corresponding grid node

  double
    diam,           // Diameter
    normVel,        // Norm of the relative velocity
    temp,           // Temperature
    relTemp,        // Relative temperature
    reynolds,       // Dispersed Reynold number
    nusselt,        // Nusselt number
    sherwood,       // Sherwood number
    schmidt,        // Schmidt number
    spalding,       // Spalding number
    prandtl,        // Prandtl number
    dragCoeff,      // Drag coefficient
    liftCoeff,      // Lift coefficient
    kinRespTime,    // Kinematic response time
    thermRespTime,  // Thermal response time
    massTrCoeff,    // Mass transfer coefficient beta
    pressBubble,    // Bubble interal pressure
    rhoBubble,      // Bubble internal density
    concInterf;     // Bubble surface concentration

  std::vector< double >
    pos,                      // Entity position at time step n
    posOld,                   // Entity position at time step n-1
    vel,                      // Entity velocity at time step n
    velOld,                   // Entity velocity at time step n-1
    relVel;                   // Relative velocity

  ENTITY_ELEMENT_DATA edata;  // Corresponding grid element data
};


/// Local flow variables data structure
struct LOCAL_FLOW_VARIABLES
{
  LOCAL_FLOW_VARIABLES(const int dim) :
    vel  (dim,0.),
    velDt(dim,0.),
    vort (dim,0.),
    velDx(dim,std::vector< double >(dim,0.)) {}
  ~LOCAL_FLOW_VARIABLES() {}
  std::vector< double >
    vel,       // Local instantaneous flow velocity
    velDt,     // Flow velocity time derivative
    vort;      // Local instantaneous flow vorticity
  std::vector< std::vector< double > >
    velDx;     // Flow velocity space derivative
  double
    pressure,  // Local instantaneous flow pressure
    temp;      // Local instantaneous flow temperature
};



/**
 * structure to hold material physical properties
 */
struct PLAS_MATERIAL_DATA
{
  // methods
  PLAS_MATERIAL_DATA() : T(273.15/*K*/), p(101325./*Pa*/) {}
  virtual void update(double _T, double _p) { T=_T; p=_p; }

  // independent physical properties
  double
    T,
    p;
  plas_flowtype_t flowtype;  // type of flow (particle, droplet, bubble)

  // generic (dependent) physical properties
  double
    rho,  // density
    mu,   // dynamic viscosity
    cp,   // specific heat capacity
    k;    // thermal conductivity

  // dispersed phase-specific properties
  double
    sig,              // surface tension
    eps,              // emissivity
    satPres,          // saturation pressure evaluated at surface of entities
    vapPres,          // vapour pressure
    latHeat,          // specific latent heat
    molarMass,        // molar mass
    molarMassVap,     // molar mass at vapour phase
    binaryDiffCoeff,  // binary diffusion coefficient
    massDiffCoeff,    // mass diffusivity for a gas in a liquid
    He;               // Henry law constant
};


/// PLaS interface class
class plas {


 public:  // interface public methods

   /**
    * initialize simulation-independent solver data
    */
   plas();


   /**
   * This function initializes the PLaS solver. It has to be
   * called before doing any run of PLaS from the driving flow
   * solver.
   */
  void initialize(const XMLNode& x);

  /**
   * This is the main routine of the PLaS solver. It has to be
   * called at each time step of the driving flow solver.
   */
  void run();

  /**
   * This routine terminates PLaS. It has to be called after
   * the last run of PLaS from the driving flow solver.
   */
  virtual ~plas();


 private:  // interface methods to implement in derived class

  /**
   * Set flow solver parameters on initialization
   * @param[in] fp pointer to PLaS data structure
   */
  virtual void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp) = 0;

  /**
   * Provide element node indices to PLaS
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @param[out] enodes element node indices array (preallocated)
   */
  virtual void getElmNodes(const int iz, const int ie, int *enodes) = 0;

  /**
   * Provide element type to PLaS
   * @param[in] element zone
   * @param[in] element index in the zone
   * @return element type
   */
  virtual plas_elmtype_t getElmType(const int iz, const int ie) = 0;

  /**
   * Provide information about which boundary is a wall to PLaS
   * @param[in] bnd boundary index
   * @return non-zero if boundary is a wall, else zero
   */
  virtual int getWallBndFlag(int bnd) = 0;

  /**
   * Provide nodal quantity at time step n-1 (previous) to PLaS
   * @param[in] Q quantity to provide (as quantity enumerated index)
   * @param[in] nod node index
   * @param[out] v quantity value
   * @param[in] d quantity dimension, for tensor quantities
   * @return first entry in quantity value pointer (for scalar quantities)
   */
  virtual double getOldQuantity(const plas_quantity_t& Q, const int nod, double *v=NULL, const int d=1) = 0;

  /**
   * Provide nodal quantity at time step n (current) to PLaS
   * @param[in] Q quantity to provide (as quantity enumerated index)
   * @param[in] nod node index
   * @param[out] v quantity value
   * @param[in] d quantity dimension, for tensor quantities
   * @return first entry in quantity value pointer (for scalar quantities)
   */
  virtual double getQuantity(const plas_quantity_t& Q, const int nod, double *v=NULL, const int d=1) = 0;

  /**
   * Provide Eulerian time scale (turbulence kinetic energy over dissipation)
   * for a node
   * @param[in] nod node index
   * @return Eulerian time scale (k/eps)
   */
  virtual double getEulerianTimeScale(int nod) = 0;

  /**
   * Pass data for screen output
   * @param[in] text screen output
   */
  virtual void screenOutput(const std::string& text) = 0;

  /**
   * Pass data for a screen warning
   * @param[in] text screen warning
   */
  virtual void screenWarning(const std::string& text) = 0;


 protected:  // internal methods

  /**
   * This file includes the computation of back-coupling terms
   * from the dispersed phase to the continuous phase.
   *
   * This function computes the mass and momentum back coupling
   * terms for a dispersed entity.
   */
  void plas_CalcBackCoupling(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double *force, double tFactor);


  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function calculates the unit normal vector of a
   * boundary face.
   */
  void plas_CalcBoundaryUnitNormal(int numDim, int ibnd, int ifac, double *unitVec);


  /**
   * This file includes a function to manage cellwise data of
   * the secondary phase.
   *
   * This function computes the cellwise data (e.g. the volume
   * fraction) of the secondary phase and corrects the data
   * on the boundaries of a multiprocessor job.
   */
  void plas_CalcCellwiseData();


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the concentation in bubble's surface.
   */
  double plas_CalcConcInterf(double pressBubble);


  /**
   * This file includes the computation of back-coupling terms
   * from the dispersed phase to the continuous phase.
   *
   * This function computes the momentum back coupling forces
   * for a bubble.
   */
  void plas_CalcCouplingForcesBubble(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor);


  /**
   * This file includes the computation of back-coupling terms
   * from the dispersed phase to the continuous phase.
   *
   * This function calculates the momentum back coupling forces
   * for a particle or droplet.
   */
  void plas_CalcCouplingForcesParticle(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor);


  /**
   * This function computes the 3D cross product of two vectors.
   */
  void plas_CalcCrossProduct_3D(double *value, double *a, double *b);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Reynolds number.
   */
  double plas_CalcDispReynolds(double viscosity, double diameter, double normVel);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity drag coefficient.
   */
  double plas_CalcDragCoeff(int flowType, double reynolds);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes all flow coefficients.
   */
  void plas_CalcEntityCoefficients(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the kinematic response time of the
   * dispersed entity.
   */
  double plas_CalcKinematicResponseTime(LOCAL_ENTITY_VARIABLES *ent);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity lift coefficient.
   */
  double plas_CalcLiftCoeff(int flowType);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity mass transfer
   * coefficient.
   */
  double plas_CalcMassTransferCoeff(double sherwood, double spalding);


  /**
   * This function computes the 2D scalar product of a matrix
   * and a vector.
   */
  void plas_CalcMatVectScalarProduct_2D(double *value, double **m, double *a);


  /**
   * This function computes the 3D scalar product of a matrix
   * and a vector.
   */
  void plas_CalcMatVectScalarProduct_3D(double *value, double **m, double *a);


  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function computes the node impact factors of an
   * element according to the entity position.
   */
  void plas_CalcNodeImpactFactors(const LOCAL_ENTITY_VARIABLES *ent, double *imp);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Nusselt number.
   */
  double plas_CalcNusseltNumber(int evapModel, double reynolds, double spalding, double prandtl);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Prandtl number.
   */
  double plas_CalcPrandtlNumber();


  /**
   * This file includes functions to compute flow coefficients.
   */
  double plas_CalcPressBubble(double diameter, double pressure);


  /**
   * This file includes functions to compute flow coefficients.
   */
  double plas_CalcRhoBubble(double temperature, double pressBubble);


  /**
   * This function computes a rotation matrix in 2D.
   */
  void plas_CalcRotationMatrix_2D(double phi, double **m);


  /**
   * This function computes a rotation matrix in 3D.
   */
  void plas_CalcRotationMatrix_3D(double phi, double **m, int axis);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Schmidt number.
   */
  double plas_CalcSchmidtNumber();


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Sherwood number.
   */
  double plas_CalcSherwoodNumber(int evapModel, double reynolds, double schmidt, double spalding);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Spalding number.
   */
  double plas_CalcSpaldingNumber(double pressure);


  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the thermal response time of the
   * dispersed entity.
   */
  double plas_CalcThermalResponseTime(double diameter);


  /**
   * This file includes the functionality to perform  trajectory
   * integrations of dispersed entities.
   *
   * This function composes the trajectory equation.
   */
  void plas_CalcTrajectory(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double dtLagr);


  /**
   * This function computes the angle between two vectors.
   */
  double plas_CalcVectorAngle(int numDim, double *a, double *b);


  /**
   * This function computes the length of a vector.
   */
  double plas_CalcVectorLength(int numDim, double *a);


  /**
   * This function rotates a 3D vector.
   */
  void plas_CalcVectorRotation_3D(double phi, double **m, double *a);


  /**
   * This function computes the scalar product of two vectors.
   */
  double plas_CalcVectorScalarProduct(int numDim, double *a, double *b);


  /**
   * Calculates the size of a triangle (area)
   * @param[in] n1 element node index 1 (Gambit convention)
   * @param[in] n2 element node index 2
   * @param[in] n3 element node index 3
   * @return area of the triangle
   */
  double plas_CalcSizeTriangle(const unsigned n1, const unsigned n2, const unsigned n3);


  /**
   * Calculates the size of a tetrahedron (volume)
   * @param[in] n1 element node index 1 (Gambit convention)
   * @param[in] n2 element node index 2
   * @param[in] n3 element node index 3
   * @param[in] n4 element node index 4
   * @return volume of the tetrahedron
   */
  double plas_CalcSizeTetrahedron(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4);


  /**
   * Calculates the size (volume/area) of an element
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @return volume/area of the given element
   */
  double plas_CalcElmSize(const unsigned iz, const unsigned ie);


  /**
   * This function computes the vorticity of the fluid flow.
   */
  void plas_CalcVorticity(int numDim, LOCAL_FLOW_VARIABLES *flow);


  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function calculates the distance of an entity to a
   * wall face.
   */
  double plas_CalcWallFaceDistance(int numDim, double *pos, int ibnd, int ifac);


  /**
   * This file includes the functionality to perform  trajectory
   * integrations of dispersed entities.
   *
   * This function performs a not-a-number check for the entity
   * position and velocity.
   */
  void plas_CheckNaN(LOCAL_ENTITY_VARIABLES *ent);


  /**
   * This function computes coalescense of bubbles according
   * to the thin film theory.
   *
   * The routine is not used for the time being. Before using
   * it, please check carefully the implementation.
   */
  int plas_Coalescence(double dj, double di, double *uijRelPrPr, double *x, double *y, double *z, double Mi, double Mj, double *uiPrPr, double *uiPrPrNew, double *uiNew, double *ui);


  /**
   * This function computes inter-entity collisions based on the
   * stochastic model of Sommerfeld.
   */
  void plas_CollisionModel(LOCAL_ENTITY_VARIABLES *ent, int numDens, double dtLagr);


  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function writes the PLaS statistics to a file.
   */
  void plas_CreateStatsFile(const std::string &outpString);


  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function creates a Tecplot file of the dispersed
   * phase entities.
   */
  void plas_CreateTecplotFile(const std::string &outpString);


  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function finds the index of the boundary face through
   * which an entity left the computational domain.
   */
  bool plas_FindExitFace(LOCAL_ENTITY_VARIABLES *ent, int *i, int *j);


  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function finds the smallest element face distance to
   * an entity.
   */
  void plas_FindMinimumElementFaceDistance(int numDim, LOCAL_ENTITY_VARIABLES *ent, int *idx, double *dmin);


  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function finds the closest element node to an entity.
   */
  int plas_FindNearestElementNode(LOCAL_ENTITY_VARIABLES *ent);

  /**
   * Provide domain element associated to boundary face to PLaS
   * @param[in] bnd boundary index
   * @param[in] bface face of the boundary
   * @return element index
   */
  int plas_getBndDomElm(int bnd, int bface);

  /**
   * Provide component of boundary face normal to PLaS
   * @param[in] elm element index
   * @param[in] bface face of the boundary
   * @param[in] dim coordinate index
   * @return coordinate of the normal
   */
  double plas_getBndFaceNormComp(int bnd, int bface, int dim);

  /**
   * Provide component of element face middle-point vector
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @param[in] iface element face index
   * @param[out] coordinates of the face middle-point (preallocated)
   */
  void plas_getElmFaceMiddlePoint(const int iz, const int ie, const int iface, double *fmp);

  /**
   * This function gets the nodes of a boundary face of a given element
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @param[in] iface element face index
   * @param[out] fnodes element face node indices array (preallocated)
   * @return returned face element type
   */
  void plas_getElmFaceNodes(const int iz, const int ie, const int iface, int *fnodes);

  /**
   * This function gets the face type of a given element
   * @param[in] element type
   * @param[in] iface element face index
   * @return returned face element type
   */
  plas_elmtype_t plas_getElmFaceType(const plas_elmtype_t et, const int iface);

  /**
   * Provide neighbour of an element to PLaS
   * @param[in] elm element index
   * @param[in] eface face of the element
   * @return neighbour element index
   */
  int plas_getElmNeighbour(int elm, int eface);

  /**
   * Provide element number of nodes
   * @param[in] element type
   * @return element number of faces
   */
  int plas_getElmNFaces(const plas_elmtype_t et);


  /**
   * Provide element number of nodes
   * @param[in] element type
   * @return element number of nodes
   */
  int plas_getElmNNodes(const plas_elmtype_t et);


  /**
   * Provide component of element normal to PLaS
   * @param[in] elm element index
   * @param[in] eface face of the element
   * @param[in] dim coordinate index
   * @return coordinate of the normal
   */
  double plas_getElmNormComp(int elm, int eface, int dim);

  /**
   * Provide nodal quantity at time step n-1 (previous) to PLaS
   * @param[in] Q quantity to provide (as simple index)
   * @param[in] nod node index
   * @param[out] v quantity value
   * @param[in] d quantity dimension, for tensor quantities
   * @return first entry in quantity value pointer (for scalar quantities)
   */
  double plas_getOldQuantity(const int q, const int nod, double *v=NULL, const int d=1) {
    return getOldQuantity(plas_quantity_t(q),nod,v,d);
  }


  /**
   * Provide nodal quantity at time step n (current) to PLaS
   * @param[in] Q quantity to provide (as simple index)
   * @param[in] nod node index
   * @param[out] v quantity value
   * @param[in] d quantity dimension, for tensor quantities
   * @return first entry in quantity value pointer (for scalar quantities)
   */
  double plas_getQuantity(const int q, const int nod, double *v=NULL, const int d=1) {
    return getQuantity(plas_quantity_t(q),nod,v,d);
  }


  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function imposes entities that were generated by an
   * external code module (e.g. electrochemistry).
   */
  void plas_ImposeExternal();


  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function generates entites in specified spatial
   * domains of the computational space.
   */
  void plas_ImposeProductionDomains();


  /**
   * This function manages the interpolation of variables from
   * the fluid flow solver.
   */
  void plas_Interpolate(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step);


  /**
   * This function calculates all the quantities interpolation
   * @param[in] Q quantity to interpolate (as quantity enumerated index)
   * @param[in] ent local entity description
   * @param[in] fold factor to weigh the quantity old value with
   * @param[in] new factor to weigh the quantity new value with
   * @return quantity interpolated value
   */
  double plas_InterpolateQuantity(const plas_quantity_t &Q, const LOCAL_ENTITY_VARIABLES &ent, const std::vector< double >& en_impactfactor, double fold, double fnew);


  /**
   * This function calculates all the quantities interpolation
   * @param[in] q quantity to interpolate (as simple index)
   * @param[in] ent local entity description
   * @param[in] fold factor to weigh the quantity old value with
   * @param[in] new factor to weigh the quantity new value with
   * @return quantity interpolated value
   */
  double plas_InterpolateQuantity(const int &q, const LOCAL_ENTITY_VARIABLES &ent, const std::vector< double >& en_impactfactor, double fold, double fnew) {
    return plas_InterpolateQuantity(plas_quantity_t(q),ent,en_impactfactor,fold,fnew);
  }


  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function loads a distribution of entities from a file.
   */
  void plas_LoadInitialDistribution(const std::string &inpString);


  /**
   * This function normalizes a vector to length one
   */
  void plas_NormalizeVector(int numDim, double *a);


  /**
   * This function generates random doubles in a given range
   */
  double plas_RandomDouble(double _min=0., double _max=1.);


  /**
   * Generate a random position inside an element in a given zone
   * @param[in] zone element zone index
   * @param[in] elm element index in the zone
   * @param[out] p random position (preallocated to number of dimensions)
   */
  void plas_RandomElmPosition(const int zone, const int elm, double *p);


  /**
   * This file contains routines to generate random numbers.
   *
   * This functions generates random numbers according to a
   * Gaussian distribution with mean m and standard deviation s.
   *
   * (c) Copyright 1994, Everett F. Carter Jr.
   * Permission is granted by the author to use this software
   * for any application provided this copyright notice
   * is preserved.
   */
  double plas_RandomGaussian(float m, float s);


  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function performs a random initial distribution of
   * dispersed entites.
   */
  void plas_RandomInitialDistribution();


  /**
   * This file contains routines to generate random numbers.
   *
   * This functions generates random integers.
   */
  int plas_RandomInteger(int min, int max);


  /**
   * This file contains all functionality to read in data from
   * the PLaS.conf data file.
   *
   * This function reads the PLaS.conf data file.
   */
  void plas_ReadParameters(const XMLNode& x);


  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function performs a successive neigbour search routine
   * following the algorithm of Loehner et al.
   */
  void plas_SearchSuccessive(LOCAL_ENTITY_VARIABLES *ent);


  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function sets an entity diameter according to a
   * distribution function.
   */
  double plas_SetDiameter();


  /**
   * This function sets the geometry of an element assigned to
   * a dispersed entity.
   */
  void plas_SetElementGeometry(int numDim, LOCAL_ENTITY_VARIABLES *ent);


  /**
   * This function sets the volumes of elements and additionally distributes
   * them to the nodes.
   */
  void plas_SetupElementsAndNodalVolumes();


  /**
   * This file includes the functionality to perform  trajectory
   * integrations of dispersed entities.
   *
   * This function solves the trajectory equation.
   */
  void plas_SolveGaussSeidel(int numDim, double **mat, double *s, double *rhs);


  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This routine terminates PLaS due to a fatal error. It
   * writes out an error message.
   */
  void plas_TerminateOnError(const std::string& errMessage);


  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function performs a wall bounce for an entity that
   * crossed a wall face of the domain boundary in the last
   * trajectory updata. Position and velocity are corrected.
   */
  void plas_WallBounce(int numDim, double elasticity, LOCAL_ENTITY_VARIABLES *ent, int ibnd, int ifac);


  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function writes the PLaS statistics to a file.
   */
  void plas_WriteStatsFile(const std::string &outpString, int iter, double time);


  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function appends Tecplot output of the dispersed
   * phase entities into the previously created file
   */
  void plas_WriteTecplotFile(const std::string &outpString, int iter, double time);


 public: // internal data (PLAS_DATA)

  PLAS_ENTITY_DATA      *ed;          // Entity data structure (per entity)
  PLAS_PHASE_DATA       *pd;          // Phase data structure (per node)
  PLAS_STATS            sd;           // Statistics data structure
  PLAS_INPUT_PARAM      ip;           // Input file parameter structure
  PLAS_FLOWSOLVER_PARAM fp;           // Flowsolver parameter structure
  int                   numExtEnt;    // Number of bubbles coming from external code
  double                *extEntPos;   // Positions of bubbles coming from external code
  double                *extEntVel;   // Velocities of bubbles coming from external code
  double                *extEntTemp;  // Temperatures of bubbles coming from external code
  double                *extEntDiam;  // Diameters of bubbles coming from external code
  double *massResid;                  // Mass flux residual
  PLAS_MATERIAL_DATA
    *mdd,                             // Material data, for the dispersed phase
    *mdc;                             // Material data, for the continuous phase

  std::vector< double > volumeNod;                 // nodal volume (dual mesh)
  std::vector< std::vector< double > > volumeElm;  // elements volume

};

#endif  // PLAS_PLAS_H

