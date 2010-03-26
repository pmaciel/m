//---------------------------------------------------------------------------

#include "Electrostatics.h"
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//*************************************************************************//
// ELECTRONEUTRALITY CONDITION
//*************************************************************************//

//--- METHODS ---------------------------------------------------------------
double Electrostatics_Electroneutrality::calcConcentrationFactor(unsigned j) const
{
  return electrolyteSolution->getIonChargeNumber(j);
}
//---------------------------------------------------------------------------
double Electrostatics_Electroneutrality::calcPotentialFactor() const
{
  return 0.;
}
//---------------------------------------------------------------------------


//*************************************************************************//
// POISSON EQUATION
//*************************************************************************//

//--- METHODS ---------------------------------------------------------------
double Electrostatics_Poisson::calcConcentrationFactor(unsigned j) const
{
//  return electrolyteSolution->getIonChargeNumber(j)*F_CONST;
  return electrolyteSolution->getIonChargeNumber(j)*F_CONST/electrolyteSolution->getSolventDielectricConstant();
}
//---------------------------------------------------------------------------
double Electrostatics_Poisson::calcPotentialFactor() const
{
//  return electrolyteSolution->getSolventDielectricConstant();
  return 1.;
}
//---------------------------------------------------------------------------

