//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------

#include "ElectrolyteModel_Ideal.h"
#include <iostream>
#include <math.h>
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrolyteModel_Ideal::ElectrolyteModel_Ideal(ElectrolyteSolution* electrolyteSolution_)
  : ElectrolyteModel(electrolyteSolution_)
{
  conductivityCorrectionFactor = 1.;
  nIons = electrolyteSolution->getNIons();
  z = new int[nIons];
  c = new double[nIons];
  D = new double[nIons];
}
//---------------------------------------------------------------------------
ElectrolyteModel_Ideal::~ElectrolyteModel_Ideal()
{
  delete[] z;
  delete[] c;
  delete[] D;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_Ideal::init(bool verbose)
{
  T = electrolyteSolution->getSolutionTemperature();
  for (unsigned i=0; i<nIons; i++)
  {
    z[i] = electrolyteSolution->getIonChargeNumber(i);
    c[i] = electrolyteSolution->getIonConcentration(i);
    D[i] = electrolyteSolution->getIonDiffusionConstant(i);
  }
}
//---------------------------------------------------------------------------
// Dij
double ElectrolyteModel_Ideal::calcDiffusionFactor(unsigned i, unsigned j) const
{
  return kroneck(i,j)*D[i]*conductivityCorrectionFactor;
}
//---------------------------------------------------------------------------
// Wi = zi*e*ui
double ElectrolyteModel_Ideal::calcMigrationFactor(unsigned i) const
{
  return F_CONST*z[i]*D[i]*conductivityCorrectionFactor/(R_CONST*T);
}
//---------------------------------------------------------------------------
// Conductivity
double ElectrolyteModel_Ideal::calcConductivity() const
{
  double kappa = 0.;
  for (unsigned i=0; i<nIons; i++) {
    kappa += c[i]*D[i]*z[i]*z[i];
  }
  return kappa*conductivityCorrectionFactor*F_CONST*F_CONST/(R_CONST*T);
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Ideal::calcOnsagerCoefficient_SF(unsigned i, unsigned j) const
{
  return c[i]*NA_CONST*kroneck(i,j)*D[i]*conductivityCorrectionFactor;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Ideal::calcActivityCoefficient_MM(unsigned i) const
{
  return 1.;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Ideal::calcActivityCoefficientDerivative_MM(unsigned i, unsigned j) const
{
  return 0.;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Ideal::calcOsmoticCoefficient_MM(double totalConcentration) const
{
  return 1.;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Ideal::calcMeanActivityCoefficient_MM(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const
{
  return 1.;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_Ideal::calcBinaryOnsagerCoefficients_SF(int* stoichCation, int* stoichAnion)
{
  Ls00 = D[0]*conductivityCorrectionFactor*c[0];
  Ls01 = 0.;
  Ls11 = D[1]*conductivityCorrectionFactor*c[1];
  unsigned nAssociatedIons = nIons-2;
  for (unsigned i=0; i<nAssociatedIons; i++)
  {
    unsigned i2 = i+2;
    Ls00 += D[i2]*conductivityCorrectionFactor*c[i2]*stoichCation[i]*stoichCation[i];
    Ls11 += D[i2]*conductivityCorrectionFactor*c[i2]*stoichAnion[i]*stoichAnion[i];
  }
  Ls00 /= R_CONST*T;
  Ls11 /= R_CONST*T;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Ideal::calcConductivityCorrectionFactor(double experimentalConductivity)
{
  conductivityCorrectionFactor = experimentalConductivity/calcConductivity();
  return conductivityCorrectionFactor;
}
//---------------------------------------------------------------------------

