#ifndef PILLAZ_PILLAZ_H
#define PILLAZ_PILLAZ_H


#include <string>
#include <vector>
#include "ext/xmlParser.h"


/// Preprocessor constants
#define DFLAG_DISABLED    0
#define DFLAG_ENABLED     1
#define DFLAG_CREATED     2
#define DFLAG_LEFT        3
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
enum pillaz_flowtype_t { FLOW_PARTIC, FLOW_DROPLET, FLOW_BUBBLY };


/// supported element types (following Gambit convention numbering)
enum pillaz_elmtype_t { ELM_UNDEFINED, ELM_TRIANGLE, ELM_TETRAHEDRON, ELM_WEDGE, ELM_QUAD, ELM_BRICK, ELM_PYRAMID, ELM_EDGE, ALL_ELEMENTS };


/// Available quantities (coordinates, pressure, temperature, velocity and its spacial derivatives)
enum pillaz_quantityscalar_t { QSCALAR_UNDEFINED, PRESSURE, TEMPERATURE, ALL_QSCALAR };
enum pillaz_quantityvector_t { QVECTOR_UNDEFINED, COORD, VELOCITY, VELOCITY_X_D, VELOCITY_Y_D, VELOCITY_Z_D, ALL_QVECTOR };


/// Description of an element address
struct pillaz_elementaddress {
  pillaz_elementaddress(const int _izone, const int _ielem, const int _iface=-1) : izone(_izone), ielem(_ielem), iface(_iface) {}
  pillaz_elementaddress() : izone(-1), ielem(-1), iface(-1) { reset(); }
  int
    izone,
    ielem,
    iface;
  void reset() { izone = ielem = iface = -1; }
  bool valid() const { return (izone>=0 && ielem>=0); }
  pillaz_elementaddress& operator=(const pillaz_elementaddress& other) {
    this->izone = other.izone;
    this->ielem = other.ielem;
    this->iface = other.iface;
    return *this;
  }
  bool operator==(const pillaz_elementaddress& other) {
    // does not compare faces
    return (this->izone==other.izone && this->ielem==other.ielem);
  }
};


/// Data of a dispersed entity (particle, droplet or bubble)
struct PILLAZ_ENTITY_DATA
{
  pillaz_elementaddress eaddress;
  int flag;             // Entity enabled/disabled flag
  int node;             // Nearest grid node to the entity
  double diameter;      // Entity diameter
  double *position;     // Entity position coordiantes
  double *velocity;     // Entity velocity components
  double *velocityOld;  // Entity velocity components of time step n-1
  double temperature;   // Entity temperature;
};


/// Data of the dispersed phase (per node)
struct PILLAZ_PHASE_DATA
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


/// Statistics data of Pillaz
struct PILLAZ_STATS
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


/// Pillaz input parameters (read from Pillaz input file)
struct PILLAZ_INPUT_PARAM
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
struct PILLAZ_FLOWSOLVER_PARAM
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

  // updated on calls to pillaz::run()
  int iter;          // Current iteration
  double time;       // Current time
};


/// Entity element data structures
struct ENTITY_ELEMENT_DATA
{
  std::vector< int > elmNodes;                          // element nodes
  std::vector< std::vector< double > > elmFaceVectors;  // element face middle points
  std::vector< std::vector< double > > elmNorms;        // element face normals
};


/// Local entity variables data structure
struct LOCAL_ENTITY_VARIABLES
{
  LOCAL_ENTITY_VARIABLES(const int dim) :
    node(-1),
    pos   (dim,0.),
    posOld(dim,0.),
    vel   (dim,0.),
    velOld(dim,0.),
    relVel(dim,0.) {}
  ~LOCAL_ENTITY_VARIABLES() {}
  int
    flag,           // Flag to determine if the entity is active
    node;           // Number of corresponding grid node
  pillaz_elementaddress eaddress;

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
struct PILLAZ_MATERIAL_DATA
{
  // methods
  PILLAZ_MATERIAL_DATA() : T(273.15/*K*/), p(101325./*Pa*/) {}
  virtual void update(double _T, double _p) { T=_T; p=_p; }

  // independent physical properties
  double
    T,
    p;
  pillaz_flowtype_t flowtype;  // type of flow (particle, droplet, bubble)

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


/// Pillaz interface class
class pillaz {


 public:  // interface public methods

   /**
   * This function initializes the Pillaz solver. It has to be
   * called before doing any run of Pillaz from the driving flow
   * solver.
   */
  void initialize(const XMLNode& x);

  /**
   * This is the main routine of the Pillaz solver. It has to be
   * called at each time step of the driving flow solver.
   */
  void run();

  /**
   * This routine terminates Pillaz. It has to be called after
   * the last run of Pillaz from the driving flow solver.
   */
  virtual ~pillaz();


 private:  // interface methods to implement in derived class

  /**
   * Set flow solver parameters on initialization
   * @param[in] fp pointer to Pillaz data structure
   */
  virtual void setFlowSolverParamOnInit(PILLAZ_FLOWSOLVER_PARAM *_fp) = 0;

  /**
   * Provide boundary element node indices to Pillaz
   * @param[in] iz boundary element zone
   * @param[in] ie element index in the boundary zone
   * @param[out] enodes element node indices array (preallocated)
   */
  virtual void getBndElmNodes(const int iz, const int ie, int *enodes) = 0;

  /**
   * Provide boundary element type to Pillaz
   * @param[in] element boundary zone
   * @param[in] element index in the boundary zone
   * @return element type
   */
  virtual pillaz_elmtype_t getBndElmType(const int iz, const int ie) = 0;

  /**
   * Provide inner element node indices to Pillaz
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @param[out] enodes element node indices array (preallocated)
   */
  virtual void getElmNodes(const int iz, const int ie, int *enodes) = 0;

  /**
   * Provide inner element type to Pillaz
   * @param[in] element inner zone
   * @param[in] element index in the inner zone
   * @return element type
   */
  virtual pillaz_elmtype_t getElmType(const int iz, const int ie) = 0;

  /**
   * Provide information about which boundary is a wall to Pillaz
   * @param[in] bnd boundary index
   * @return non-zero if boundary is a wall, else zero
   */
  virtual int getWallBndFlag(int bnd) = 0;

  /**
   * Provide nodal scalar quantity at time step n (current) to Pillaz
   * @param[in] Q scalar quantity to provide
   * @param[in] nod node index
   * @return first entry in quantity value pointer (for scalar quantities)
   */
  virtual double getQuantityScalar(const pillaz_quantityscalar_t& Q, const int nod) = 0;

  /**
   * Provide nodal scalar quantity at time step n-1 (previous) to Pillaz
   * (default implementation is not time-accurate)
   * @param[in] Q scalar quantity to provide
   * @param[in] nod node index
   * @return first entry in quantity value pointer (for scalar quantities)
   */
  virtual double getQuantityScalarOld(const pillaz_quantityscalar_t& Q, const int nod) {
    return getQuantityScalar(Q,nod);
  }

  /**
   * Provide nodal vector quantity at time step n (current) to Pillaz
   * @param[in] Q vector quantity to provide
   * @param[in] nod node index
   * @param[in] dim quantity dimension
   * @param[out] v quantity value
   */
  virtual void getQuantityVector(const pillaz_quantityvector_t& Q, const int nod, const int dim, double *v) = 0;

  /**
   * Provide nodal vector quantity at time step n-1 (previous) to Pillaz
   * (default implementation is not time-accurate)
   * @param[in] Q vector quantity to provide
   * @param[in] nod node index
   * @param[in] dim quantity dimension
   * @param[out] v quantity value
   */
  virtual void  getQuantityVectorOld(const pillaz_quantityvector_t& Q, const int nod, const int dim, double *v) {
    getQuantityVector(Q,nod,dim,v);
  }

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
  void pillaz_CalcBackCoupling(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double *force, double tFactor);

  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function calculates the unit normal vector of a
   * boundary face.
   */
  void pillaz_CalcBoundaryUnitNormal(int numDim, int ibnd, int ifac, double *unitVec);

  /**
   * This file includes a function to manage cellwise data of
   * the secondary phase.
   *
   * This function computes the cellwise data (e.g. the volume
   * fraction) of the secondary phase and corrects the data
   * on the boundaries of a multiprocessor job.
   */
  void pillaz_CalcCellwiseData();

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the concentation in bubble's surface.
   */
  double pillaz_CalcConcInterf(double pressBubble);

  /**
   * This file includes the computation of back-coupling terms
   * from the dispersed phase to the continuous phase.
   *
   * This function computes the momentum back coupling forces
   * for a bubble.
   */
  void pillaz_CalcCouplingForcesBubble(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor);

  /**
   * This file includes the computation of back-coupling terms
   * from the dispersed phase to the continuous phase.
   *
   * This function calculates the momentum back coupling forces
   * for a particle or droplet.
   */
  void pillaz_CalcCouplingForcesParticle(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double tFactor);

  /**
   * This function computes the 3D cross product of two vectors.
   */
  void pillaz_CalcCrossProduct_3D(double *value, double *a, double *b);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Reynolds number.
   */
  double pillaz_CalcDispReynolds(double viscosity, double diameter, double normVel);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity drag coefficient.
   */
  double pillaz_CalcDragCoeff(int flowType, double reynolds);

  /**
   * Caculate area/length-scaled element face normal
   * @param[in] element zone index
   * @param[in] elm element index
   * @param[in] eface face of the element
   * @param[out] normal area/length-scaled face normal vector
   */
  void pillaz_CalcElmFaceNormal(const int zone, const int elm, const int eface, double *normal);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes all flow coefficients.
   */
  void pillaz_CalcEntityCoefficients(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the kinematic response time of the
   * dispersed entity.
   */
  double pillaz_CalcKinematicResponseTime(LOCAL_ENTITY_VARIABLES *ent);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity lift coefficient.
   */
  double pillaz_CalcLiftCoeff(int flowType);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity mass transfer
   * coefficient.
   */
  double pillaz_CalcMassTransferCoeff(double sherwood, double spalding);

  /**
   * This function computes the 2D scalar product of a matrix
   * and a vector.
   */
  void pillaz_CalcMatVectScalarProduct_2D(double *value, double **m, double *a);

  /**
   * This function computes the 3D scalar product of a matrix
   * and a vector.
   */
  void pillaz_CalcMatVectScalarProduct_3D(double *value, double **m, double *a);

  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function computes the node impact factors of an
   * element according to the entity position.
   */
  void pillaz_CalcNodeImpactFactors(const LOCAL_ENTITY_VARIABLES *ent, double *imp);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Nusselt number.
   */
  double pillaz_CalcNusseltNumber(int evapModel, double reynolds, double spalding, double prandtl);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Prandtl number.
   */
  double pillaz_CalcPrandtlNumber();

  /**
   * This file includes functions to compute flow coefficients.
   */
  double pillaz_CalcPressBubble(double diameter, double pressure);

  /**
   * This file includes functions to compute flow coefficients.
   */
  double pillaz_CalcRhoBubble(double temperature, double pressBubble);

  /**
   * This function computes a rotation matrix in 2D.
   */
  void pillaz_CalcRotationMatrix_2D(double phi, double **m);

  /**
   * This function computes a rotation matrix in 3D.
   */
  void pillaz_CalcRotationMatrix_3D(double phi, double **m, int axis);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Schmidt number.
   */
  double pillaz_CalcSchmidtNumber();

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Sherwood number.
   */
  double pillaz_CalcSherwoodNumber(int evapModel, double reynolds, double schmidt, double spalding);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the entity Spalding number.
   */
  double pillaz_CalcSpaldingNumber(double pressure);

  /**
   * This file includes functions to compute flow coefficients.
   *
   * This function computes the thermal response time of the
   * dispersed entity.
   */
  double pillaz_CalcThermalResponseTime(double diameter);

  /**
   * This file includes the functionality to perform  trajectory
   * integrations of dispersed entities.
   *
   * This function composes the trajectory equation.
   */
  void pillaz_CalcTrajectory(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double dtLagr);

  /**
   * This function computes the angle between two vectors.
   */
  double pillaz_CalcVectorAngle(int numDim, double *a, double *b);

  /**
   * This function computes the length of a vector.
   */
  double pillaz_CalcVectorLength(int numDim, double *a);

  /**
   * This function rotates a 3D vector.
   */
  void pillaz_CalcVectorRotation_3D(double phi, double **m, double *a);

  /**
   * This function computes the scalar product of two vectors.
   */
  double pillaz_CalcVectorScalarProduct(int numDim, double *a, double *b);

  /**
   * Calculates the size of a triangle (area)
   * @param[in] n1 element node index 1 (Gambit convention)
   * @param[in] n2 element node index 2
   * @param[in] n3 element node index 3
   * @return area of the triangle
   */
  double pillaz_CalcSizeTriangle(const unsigned n1, const unsigned n2, const unsigned n3);

  /**
   * Calculates the size of a tetrahedron (volume)
   * @param[in] n1 element node index 1 (Gambit convention)
   * @param[in] n2 element node index 2
   * @param[in] n3 element node index 3
   * @param[in] n4 element node index 4
   * @return volume of the tetrahedron
   */
  double pillaz_CalcSizeTetrahedron(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4);

  /**
   * Calculates the size (volume/area) of an element
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @return volume/area of the given element
   */
  double pillaz_CalcElmSize(const unsigned iz, const unsigned ie);

  /**
   * This function computes the vorticity of the fluid flow.
   */
  void pillaz_CalcVorticity(int numDim, LOCAL_FLOW_VARIABLES *flow);

  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function calculates the distance of an entity to a
   * wall face.
   */
  double pillaz_CalcWallFaceDistance(int numDim, double *pos, int ibnd, int ifac);

  /**
   * This file includes the functionality to perform  trajectory
   * integrations of dispersed entities.
   *
   * This function performs a not-a-number check for the entity
   * position and velocity.
   */
  void pillaz_CheckNaN(LOCAL_ENTITY_VARIABLES *ent);

  /**
   * This function computes coalescense of bubbles according
   * to the thin film theory.
   *
   * The routine is not used for the time being. Before using
   * it, please check carefully the implementation.
   */
  int pillaz_Coalescence(double dj, double di, double *uijRelPrPr, double *x, double *y, double *z, double Mi, double Mj, double *uiPrPr, double *uiPrPrNew, double *uiNew, double *ui);

  /**
   * This function computes inter-entity collisions based on the
   * stochastic model of Sommerfeld.
   */
  void pillaz_CollisionModel(LOCAL_ENTITY_VARIABLES *ent, int numDens, double dtLagr);

  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function writes the Pillaz statistics to a file.
   */
  void pillaz_CreateStatsFile(const std::string &outpString);

  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function creates a Tecplot file of the dispersed
   * phase entities.
   */
  void pillaz_CreateTecplotFile(const std::string &outpString);

  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function finds the index of the boundary face through
   * which an entity left the computational domain.
   */
  bool pillaz_FindExitFace(LOCAL_ENTITY_VARIABLES *ent, int *i, int *j);

  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function finds the smallest element face distance to
   * an entity.
   */
  void pillaz_FindMinimumElementFaceDistance(int numDim, LOCAL_ENTITY_VARIABLES *ent, int *idx, double *dmin);

  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function finds the closest element node to an entity.
   */
  int pillaz_FindNearestElementNode(LOCAL_ENTITY_VARIABLES *ent);

  /**
   * Provide component of element face middle-point vector
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @param[in] iface element face index
   * @param[out] coordinates of the face middle-point (preallocated)
   */
  void pillaz_getElmFaceMiddlePoint(const int iz, const int ie, const int iface, double *fmp);

  /**
   * This function gets the nodes of a boundary face of a given element
   * @param[in] iz element zone
   * @param[in] ie element index in the zone
   * @param[in] iface element face index
   * @param[out] fnodes element face node indices array (preallocated)
   * @return returned face element type
   */
  void pillaz_getElmFaceNodes(const int iz, const int ie, const int iface, int *fnodes);

  /**
   * This function gets the face type of a given element
   * @param[in] element type
   * @param[in] iface element face index
   * @return returned face element type
   */
  pillaz_elmtype_t pillaz_getElmFaceType(const pillaz_elmtype_t et, const int iface);

  /**
   * Provide element number of nodes
   * @param[in] element type
   * @return element number of faces
   */
  int pillaz_getElmNFaces(const pillaz_elmtype_t et);

  /**
   * Provide element number of nodes
   * @param[in] element type
   * @return element number of nodes
   */
  int pillaz_getElmNNodes(const pillaz_elmtype_t et);

  /**
   * Provide inner element type to Pillaz
   * @param[in] ea element address (structure)
   * @return element type
   */
  pillaz_elmtype_t pillaz_getElmType(const pillaz_elementaddress& ea) {
    return getElmType(ea.izone,ea.iface);
  }

  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function imposes entities that were generated by an
   * external code module (e.g. electrochemistry).
   */
  void pillaz_ImposeExternal();

  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function generates entites in specified spatial
   * domains of the computational space.
   */
  void pillaz_ImposeProductionDomains();

  /**
   * This function manages the interpolation of variables from
   * the fluid flow solver.
   */
  void pillaz_Interpolate(LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step);

  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function loads a distribution of entities from a file.
   */
  void pillaz_LoadInitialDistribution(const std::string &inpString);

  /**
   * This function normalizes a vector to length one
   */
  void pillaz_NormalizeVector(int numDim, double *a);

  /**
   * This function generates random doubles in a given range
   */
  double pillaz_RandomDouble(double _min=0., double _max=1.);

  /**
   * Generate a random position inside an element in a given zone
   * @param[in] zone element zone index
   * @param[in] elm element index in the zone
   * @param[out] p random position (preallocated to number of dimensions)
   */
  void pillaz_RandomElmPosition(const int zone, const int elm, double *p);

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
  double pillaz_RandomGaussian(float m, float s);

  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function performs a random initial distribution of
   * dispersed entites.
   */
  void pillaz_RandomInitialDistribution();

  /**
   * This file contains routines to generate random numbers.
   *
   * This functions generates random integers.
   */
  int pillaz_RandomInteger(int min, int max);

  /**
   * This file contains all functionality to read in data from
   * the Pillaz.conf data file.
   *
   * This function reads the Pillaz.conf data file.
   */
  void pillaz_ReadParameters(const XMLNode& x);

  /**
   * This file includes all functinality to perform an element
   * search for a dispersed entity.
   *
   * This function performs a successive neigbour search routine
   * following the algorithm of Loehner et al.
   */
  void pillaz_SearchSuccessive(LOCAL_ENTITY_VARIABLES *ent);

  /**
   * This file contains routines for the generation of entities,
   * be it to form an initial set of entities or to create
   * entites due to secondary phase boundary conditions.
   *
   * This function sets an entity diameter according to a
   * distribution function.
   */
  double pillaz_SetDiameter();

  /**
   * This function sets the geometry of an element assigned to
   * a dispersed entity.
   */
  void pillaz_SetElementGeometry(int numDim, LOCAL_ENTITY_VARIABLES *ent);

  /**
   * This file includes the functionality to perform  trajectory
   * integrations of dispersed entities.
   *
   * This function solves the trajectory equation.
   */
  void pillaz_SolveGaussSeidel(int numDim, double **mat, double *s, double *rhs);

  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This routine terminates Pillaz due to a fatal error. It
   * writes out an error message.
   */
  void pillaz_TerminateOnError(const std::string& errMessage);

  /**
   * This file includes all functionality concerning entities
   * interacting with boundaries of the computational domain,
   * such as wall bounces and entities leaving through outlets.
   *
   * This function performs a wall bounce for an entity that
   * crossed a wall face of the domain boundary in the last
   * trajectory updata. Position and velocity are corrected.
   */
  void pillaz_WallBounce(int numDim, double elasticity, LOCAL_ENTITY_VARIABLES *ent, int ibnd, int ifac);

  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function writes the Pillaz statistics to a file.
   */
  void pillaz_WriteStatsFile(const std::string &outpString, int iter, double time);

  /**
   * This file contains all write functionality, to the screen,
   * to output files and to Tecplot.
   *
   * This function appends Tecplot output of the dispersed
   * phase entities into the previously created file
   */
  void pillaz_WriteTecplotFile(const std::string &outpString, int iter, double time);


 public: // internal data (PILLAZ_DATA)

  PILLAZ_ENTITY_DATA      *ed;          // Entity data structure (per entity)
  PILLAZ_PHASE_DATA       *pd;          // Phase data structure (per node)
  PILLAZ_STATS            sd;           // Statistics data structure
  PILLAZ_INPUT_PARAM      ip;           // Input file parameter structure
  PILLAZ_FLOWSOLVER_PARAM fp;           // Flowsolver parameter structure
  int                   numExtEnt;    // Number of bubbles coming from external code
  double                *extEntPos;   // Positions of bubbles coming from external code
  double                *extEntVel;   // Velocities of bubbles coming from external code
  double                *extEntTemp;  // Temperatures of bubbles coming from external code
  double                *extEntDiam;  // Diameters of bubbles coming from external code
  double *massResid;                  // Mass flux residual
  PILLAZ_MATERIAL_DATA
    *mdd,                             // Material data, for the dispersed phase
    *mdc;                             // Material data, for the continuous phase

  std::vector< double > volumeNod;                 // nodal volume (dual mesh)
  std::vector< std::vector< double > > volumeElm;  // elements volume

  // inner element faces normal vectors
  std::vector<        // per (inner) zone
    std::vector<      // per element
      std::vector<    // per face
        std::vector<  // per coordinate
          double      // coordinate value
  > > > > m_innerelem_normals;

  // map of boundary elements to (inner) elements
  std::vector<               // per (boundary) zone
    std::vector<             // per element
      pillaz_elementaddress  // (inner) element address
  > > m_boundtoinner;

  // map of nodes to (inner) elements
  std::vector<               // per node
    std::vector<             // list of neighbors
      pillaz_elementaddress  // (inner) element address
  > > m_nodetoelem;

  // map of (inner) elements to elements (sharing a face)
  std::vector<                 // per (inner) zone
    std::vector<               // per element
      std::vector<             // per face
        pillaz_elementaddress  // (inner) element address
  > > > m_elemtoelem;

};

#endif  // PILLAZ_PILLAZ_H

