//---------------------------------------------------------------------------

#ifndef MathematicsPhysicsH
#define MathematicsPhysicsH

#define _USE_MATH_DEFINES

#include <cmath>

//---------------------------------------------------------------------------

const double kB_CONST = 1.3806503e-23;      // Boltzmann's constant [J/K]
const double e_CONST = 1.60217646e-19;      // elementary charge [C]
const double NA_CONST = 6.02214199e23;      // Avogadro's constant [1/mol]
const double F_CONST = e_CONST*NA_CONST;    // Faraday's constant [C/mol]
const double R_CONST = kB_CONST*NA_CONST;   // ideal gas constant [J/mol K]
const double eulergamma_CONST = 0.57721566;

#ifndef M_PI
#define M_PI 3.14159265
#endif

unsigned fac(unsigned N);
unsigned kroneck(unsigned i, unsigned j);
double   expInt(double x);
void     solveGauss(double** A, double* B, unsigned N);

//---------------------------------------------------------------------------

#endif

