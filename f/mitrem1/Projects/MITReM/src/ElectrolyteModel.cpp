//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------

#include "ElectrolyteModel.h"
#include <iostream>
#include <math.h>
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrolyteModel::ElectrolyteModel(ElectrolyteSolution* electrolyteSolution_)
  : electrolyteSolution(electrolyteSolution_)
{
}
//---------------------------------------------------------------------------
ElectrolyteModel::~ElectrolyteModel()
{
}
//---------------------------------------------------------------------------
void ElectrolyteModel::init(bool verbose)
{
  this->verbose = verbose;
  T = electrolyteSolution->getSolutionTemperature();
}
//---------------------------------------------------------------------------
double ElectrolyteModel::calcEquivalentBinaryOnsagerCoefficient_CC_SF(double normality) const
{
  return Ls00/normality;
}
//---------------------------------------------------------------------------
double ElectrolyteModel::calcEquivalentBinaryOnsagerCoefficient_CA_SF(double normality) const
{
  return Ls01/normality;
}
//---------------------------------------------------------------------------
double ElectrolyteModel::calcEquivalentBinaryOnsagerCoefficient_AA_SF(double normality) const
{
  return Ls11/normality;
}
//---------------------------------------------------------------------------
double ElectrolyteModel::calcEquivalentConductivity(double normality) const
{
  double Lambda = calcConductivity();
  return Lambda/normality;
}
//---------------------------------------------------------------------------
double ElectrolyteModel::calcMolarThermodynamicDiffusionCoefficient_SF(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const
{
  int zCation = electrolyteSolution->getIonChargeNumber(0);
  int zAnion = electrolyteSolution->getIonChargeNumber(1);
  double L = -zCation*zAnion/(sCation*sAnion) * (Ls00*Ls11 - Ls01*Ls01)/(zCation*zCation*Ls00 + 2.*zCation*zAnion*Ls01 + zAnion*zAnion*Ls11);
  return L*R_CONST*T*(sCation+sAnion)/electrolyteConcentration;
}
//---------------------------------------------------------------------------
double ElectrolyteModel::calcCationTransportNumber_SF() const
{
  int zCation = electrolyteSolution->getIonChargeNumber(0);
  int zAnion = electrolyteSolution->getIonChargeNumber(1);
  double tCation = (zCation*zCation*Ls00 + zCation*zAnion*Ls01)/(zCation*zCation*Ls00 + 2.*zCation*zAnion*Ls01 + zAnion*zAnion*Ls11);
  return tCation;
}
//---------------------------------------------------------------------------

