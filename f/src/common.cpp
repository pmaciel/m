
#include "common.h"

/* linear systems solvers */
LS *ls_coupled;  // for "coupled" system
LS *ls_scalar;   // for scalar system
LS *ls_turb;     // for turbulence system
LS *ls_turb1;    // for uncoupled turbulence system (k)
LS *ls_turb2;    // for uncoupled turbulence system (epsilon/omega)


/* data structures */
m::mmesh M;                            // mmesh structure
std::vector< m::melem > e2n;           // mmesh "inner" connectivity
std::vector< m::melem > e2n_periodic;  // mmesh periodic connectivity
bound_node_struct **Nobg;                  /* boundary-node structure */
std::vector< int                   > Ce_type;
std::vector< std::vector< double > > No_W;
std::vector< bound_node_struct > WFnodes;  // wall functions nodes structure
std::vector< bcgroup_struct    > BCgroup;  // b.c. group structure
std::vector< int    > Fab_cell;
std::vector< int    > Fab_node;
std::vector< int    > Fab_inc;
std::vector< int    > Fab_group;
std::vector< int    > No_group;
std::vector< double > No_vol;
std::vector< double > No_dt;
std::vector< double > No_nuturb;
std::vector< double > No_dissipation;
std::vector< double > No_wd;
std::vector< double > No_lenturb;
mitremassembler_struct mitremassembler;

std::vector< std::string > m_vars_label;
std::vector< double      > m_vars_init;
std::vector< double > logL1, logL2, logLi;
std::vector< double > resL1, resL2, resLi;

int iverr;          // variable for error check
int Nsys;           // number of system variables (excluding MITReM)
int Nmit;           // number of system variables (for MITReM)
int Neqns;          // number of equations solved
int Ncoupled;       // number of coupled equations
int Nvtfce;         // number of vertices per face
int Ncell;          // number of cells
int Nnode;          // number of nodes
int Nbface;         // number of boundary faces
int Nbcgroup;       // number of b.c. groups
int turb_iterinit;  // start real turbulent coefficients calculation

int iv_turb1;
int iv_turb2;
int iv_temp;
double dNvtcell;
double dNvtfce;

int Jacobian;
int dtrelax;
int dtrelax_scalar;
int dtrelax_turb;
int periodic;
int periodic_dirn;
double linrlx;
double linrlx_scalar;
double linrlx_turb;
double CFL;
double CFL_scalar;
double CFL_turb;
double periodic_pgrad;
double newton_eps;
double newton_thresh;
double conv_thresh;
double epsilon;
double C2;

std::string file_input;   // grid/solution file for input
std::string file_output;  // grid/solution file for output
std::string file_log;     // log file
std::string file_inlet;   // inlet file (some turbulent testcases)

int iter;
int Ndim;
int Nvtcell;
int Nitermax;  // number of non-linear iterations (maximum)
int Nitermin;  // ... (minimum)
int restart;
int scaconv;
int turmod;
int wall_functions;
int walldist;
int temperature;
int buoyancy;
int scalar_coupling;
int turbulence_coupling;
int diffusion;
int varden;

double grav[3];
double Umax;
double Vtot;
double turb_reflen;
double nulam;
double prandtl;
double rho;
double rhofac;
double To;
double Po;
double masstot;
double Qin;
double Qout;


/* turbulence model constants */
double Cmu;     // k-epsilon (1)
double Ceps1;   // (1)
double Ceps2;   // (1)
double Sig1;    // (1,2): sigma k
double Sig2;    // (1,2): sigma epsilon/omega
double Ceta;    // k-epsilon Launder-Sharma (3)
double Ceta2;   // (3)
double Amu;     // k-epsilon Lam-Bremhorst (4)
double Bmu;     // (4)
double Dmu;     // (4)
double Aeps1;   // (4)
double Aeps2;   // (4)
double Beps2;   // (4)
double Ck;      // k-omega (2)
double Cw1;     // (2)
double Cw2;     // (2)
double Cw;      // (2)
double Beta;    // (2)
double Betas;   // (2)
double Beta1;   // (2)
double Beta2;   // (2)
double Sigk1;   // (2)
double Sigk2;   // (2)
double Sigw1;   // (2)
double Sigw2;   // (2)
double Gamma1;  // (2)
double Gamma2;  // (2)
double A1;      // (2)

