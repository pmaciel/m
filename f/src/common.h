#ifndef __COMMON__
#define __COMMON__

#include <vector>
#include <string>
#include <cmath>

#include "ext/xmlParser.h"
#include "mfactory.h"
#include "mkernel.h"
#include "mlinearsystem.h"
#include "utils.h"

// for mitremassembler library
#include "MITReM.h"
#include "ElementMatrixAssembler.h"

// boundary conditions, scalar convection scheme and turbulence model identifiers
enum IBid { IBNONE=0,
  IBFIXV, IBFIXP, IBWALL, IBWALQ, IBSYMM, IBPERI, IBPERE, IBWLFN };
enum ISid { ISNONE=0,
  ISSFOU, ISSNSC, ISSGAL, ISSLDA, ISSLWS, ISSPSI };
enum ITid { ITNONE=0, ITMGKE=1, ITMGKW=2,
  ITKEHR=10, ITKEHG, ITKELB, ITKELS, ITKENA, ITKE2L,
  ITKWWF=20, ITKWHR, ITKWLR, ITKWBS, ITKWSS, ITKWPD };

typedef m::mlinearsystem< double > LS;

struct bcgroup_struct
{
  std::string name;  // boundary group name
  std::string zone;  // ... zone
  int type;          // ... type
  int option;        // ... type option
  int nnode;         // number of nodes
  double qflux;      // volume flux
  std::vector< double > invals;  // input values
  std::vector< double > n;       // unit normal to boundary
  std::vector< double > flux;    // integrated fluxes
  bcgroup_struct() {
    n.assign(3,0.);
    flux.assign(2,0.);
  }
};

struct local_node_struct
{
  int node;        /* node index                   */
  double W[10];    /* solution vector at nodes     */
  double Res[10];  /* residual vector at nodes     */
  double norm[3];  /* inner scaled normal          */
  double C[4];     /* scalar coefficients for cell */
  double norm2;    /* square of normal             */
  double Dt;       /* time-step contribution       */
};

struct bound_node_struct
{
  int node;     /* node index                   */
  int twin;     /* twin node for wall/periodic  */
  int deg;      /* degree (no. boundary faces)  */
  double n[3];  /* inner scaled nodal normal    */
  double modn;  /* modulus of nodal normal      */
  double dist;  /* store e.g. for wall distance */
};

struct mitremassembler_struct
{
  XMLNode x;       // xml for setup
  bool ok;         // if electrochemistry is to be run
  int iterinit;    // start calculations after this number of iterations
  int iv;          // index in W vector
  int Nions;       // number of ions
  double surfacegasfraction_min;  // surface gas fraction minimum
  double surfacegasfraction_max;  // ... maximum
  double linrelx;  // linear relaxation
  std::vector< double > bulk;  // bulk concentrations
  LS *ls;          // linear system solver
  m::mmesh Mj;     // mmesh structure for electrode reactions current density
  bool forcev;                      // if metal potentials are to be forced (and adjusted)
  double Vwelectrode, Vcelectrode,  // working/counter electrodes metal potentials (corrected)
         Awelectrode, Acelectrode;  // working/counter electrodes areas
  MITReM* m_mitrem;                     // MITReM object
  ElementMatrixAssembler* m_assembler;  // ElementMatrixAssembler object
};

/* main functions forward declarations */
void cellgeom(int ic, local_node_struct *No_loc, double *vol, int *inc_min);
void convection(local_node_struct *No_loc, double vol, int inc_min, int celltype, double *Beta2, double *LWfactor, std::vector< double >& k, std::vector< double >& sbuoy, double ***A, double ***K, double ***B, std::vector< double >& Rc);
double F1_function(double k ,double w, double nu, double y, double dkdw);
double Gfunction(const int Ndim, const double gradv[4][3]);
void Iboundary();
void IJacobian_A();
void IJacobian_N();
void IJacobian_P();
void IJacobian_scalar(int iv_scalar);
void IJacobian_turb_coupled();
void IJacobian_turb_uncoupled();
void massflux();
void read_inlet_2D(const std::string& finlet,int ig);
void read_inlet_3D(const std::string& finlet,int ig);
void readsoltp(const std::string& infile, int read_soln);
void readstart(const std::string& ccase);
void rescalc(int iv, int iv_local, double *res, int Nsize);
void scacde(int ic, int iv, local_node_struct *No_loc, double vol, int inc_min, int scheme, double diffco, double source, double coeff, int coeff_calc);
void turb_init(ITid model);
void turb_source_D(int model, double k, double turb2, double nu_l, double wd, double len, double gradkw, double *source_k, double *source_ew, double *deriv_k, double *deriv_ew, double *deriv_kew, double *deriv_ewk);
void turb_source_node();
void turb_source_P(int model, double k, double turb2, double nu_t, double nu_l, double wd, double G, double gradkw, double *source_k, double *source_ew, double v2);
double turb_viscosity(local_node_struct *No_loc, int turmod, int cell_type, double vol);
void turb_wallbc(int iv1, int iv2, LS *ls1, LS *ls2);
void turb_wfuncs(double *res);
void Update_mitrem();
void Update_temperature();
void Update_turbulence_coupled();
void Update_turbulence_sequential();
void Update_turbulence_uncoupled();
void viscous(local_node_struct *No_loc, double vol, int celltype);
void writeres();
void writesoltp(const std::string& outfile);


/* linear system solvers */
extern LS *ls_coupled;  // for "coupled" system
extern LS *ls_scalar;   // for scalar system
extern LS *ls_turb;     // for turbulence system
extern LS *ls_turb1;    // for uncoupled turbulence system (k)
extern LS *ls_turb2;    // for uncoupled turbulence system (epsilon/omega)


/* data structures */
extern m::mmesh M;                            // mmesh structure
extern std::vector< m::melem > e2n;           // mmesh "inner" connectivity
extern std::vector< m::melem > e2n_periodic;  // mmesh periodic connectivity
extern bound_node_struct **Nobg;                  /* boundary-node structure */
extern std::vector< int                   > Ce_type;
extern std::vector< std::vector< double > > No_W;
extern std::vector< bound_node_struct > WFnodes;  // wall functions nodes structure
extern std::vector< bcgroup_struct    > BCgroup;  // b.c. group structure
extern std::vector< int    > Fab_cell;
extern std::vector< int    > Fab_node;
extern std::vector< int    > Fab_inc;
extern std::vector< int    > Fab_group;
extern std::vector< int    > No_group;
extern std::vector< double > No_vol;
extern std::vector< double > No_dt;
extern std::vector< double > No_nuturb;
extern std::vector< double > No_dissipation;
extern std::vector< double > No_wd;
extern std::vector< double > No_lenturb;
extern mitremassembler_struct mitremassembler;

extern std::vector< std::string > m_vars_label;
extern std::vector< double      > m_vars_init;
extern std::vector< double > logL1, logL2, logLi;
extern std::vector< double > resL1, resL2, resLi;

extern int iverr;          // variable for error check
extern int Nsys;           // number of system variables (excluding MITReM)
extern int Nmit;           // number of system variables (for MITReM)
extern int Neqns;          // number of equations solved
extern int Ncoupled;       // number of coupled equations
extern int Nvtfce;         // number of vertices per face
extern int Ncell;          // number of cells
extern int Nnode;          // number of nodes
extern int Nbface;         // number of boundary faces
extern int Nbcgroup;       // number of b.c. groups
extern int turb_iterinit;  // start real turbulent coefficients calculation

extern int iv_turb1;
extern int iv_turb2;
extern int iv_temp;
extern double dNvtcell;
extern double dNvtfce;

extern int Jacobian;
extern int dtrelax;
extern int dtrelax_scalar;
extern int dtrelax_turb;
extern int periodic;
extern int periodic_dirn;
extern double linrlx;
extern double linrlx_scalar;
extern double linrlx_turb;
extern double CFL;
extern double CFL_scalar;
extern double CFL_turb;
extern double periodic_pgrad;
extern double newton_eps;
extern double newton_thresh;
extern double conv_thresh;
extern double epsilon;
extern double C2;


extern std::string file_input;   // grid/solution file for input
extern std::string file_output;  // grid/solution file for output
extern std::string file_log;     // log file
extern std::string file_inlet;   // inlet file (some turbulent testcases)

extern int iter;
extern int Ndim;
extern int Nvtcell;
extern int Niter;
extern int restart;
extern int scaconv;
extern int turmod;
extern int wall_functions;
extern int walldist;
extern int temperature;
extern int buoyancy;
extern int scalar_coupling;
extern int turbulence_coupling;
extern int diffusion;
extern int varden;

extern double grav[3];
extern double Umax;
extern double Vtot;
extern double turb_reflen;
extern double nulam;
extern double prandtl;
extern double rho;
extern double rhofac;
extern double To;
extern double Po;
extern double masstot;
extern double Qin;
extern double Qout;


/* turbulence model constants */
extern double Cmu;     // k-epsilon (1)
extern double Ceps1;   // (1)
extern double Ceps2;   // (1)
extern double Sig1;    // (1,2): sigma k
extern double Sig2;    // (1,2): sigma epsilon/omega
extern double Ceta;    // k-epsilon Launder-Sharma (3)
extern double Ceta2;   // (3)
extern double Amu;     // k-epsilon Lam-Bremhorst (4)
extern double Bmu;     // (4)
extern double Dmu;     // (4)
extern double Aeps1;   // (4)
extern double Aeps2;   // (4)
extern double Beps2;   // (4)
extern double Ck;      // k-omega (2)
extern double Cw1;     // (2)
extern double Cw2;     // (2)
extern double Cw;      // (2)
extern double Beta;    // (2)
extern double Betas;   // (2)
extern double Beta1;   // (2)
extern double Beta2;   // (2)
extern double Sigk1;   // (2)
extern double Sigk2;   // (2)
extern double Sigw1;   // (2)
extern double Sigw2;   // (2)
extern double Gamma1;  // (2)
extern double Gamma2;  // (2)
extern double A1;      // (2)

#endif

