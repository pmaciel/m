//---------------------------------------------------------------------------

//#define k            kineticParameters[0]
//#define cSaturation  kineticParameters[1]
//#define cDissGas    electrolyteSolution->getIonConcentration(dissolvedGas)

//---------------------------------------------------------------------------

#include "GasReaction.h"
#include <math.h>
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
GasReaction::GasReaction(ElectrolyteSolution* electrolyteSolution_)
: electrolyteSolution(electrolyteSolution_)
{
  kineticParameters = new double[2];
}
//---------------------------------------------------------------------------
GasReaction::~GasReaction()
{
  delete[] kineticParameters;
}
//---------------------------------------------------------------------------
double GasReaction::calcReactionRate() const
{
  return kineticParameters[0]*
    std::max(0., electrolyteSolution->getIonConcentration(dissolvedGas) - kineticParameters[1]);
}
//---------------------------------------------------------------------------
double GasReaction::calcReactionRateDerivativeCDissGas() const
{
  if (electrolyteSolution->getIonConcentration(dissolvedGas) > kineticParameters[1])
  {
    return kineticParameters[0];
  }
  else
  {
    return 0.;
  }
}
//---------------------------------------------------------------------------


