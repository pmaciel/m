//---------------------------------------------------------------------------

#define T        electrolyteSolution->getSolutionTemperature()
#define cRed(j)         electrolyteSolution->getIonConcentration(agentsRed[j])
#define cOxi(j)         electrolyteSolution->getIonConcentration(agentsOxi[j])
#define U        electrolyteSolution->getSolutionPotential()
#define kOxi      kineticParameters[0]
#define kRed      kineticParameters[1]
#define aOxi      kineticParameters[2]
#define aRed      kineticParameters[3]

//---------------------------------------------------------------------------

#include "ElecReaction_BV.h"
#include <math.h>
#include "MathematicsPhysics.h"
#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReaction_BV::ElecReaction_BV(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_)
: ElecReaction(electrolyteSolution_, nAgentsRed_, nAgentsOxi_)
{
  kineticParameters = new double[4];
}
//---------------------------------------------------------------------------
ElecReaction_BV::~ElecReaction_BV()
{
  delete[] kineticParameters;
}
//---------------------------------------------------------------------------
double ElecReaction_BV::calcReactionRate(double V) const
{
  double cRed = 1.;
  double cOxi = 1.;
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    cRed *= pow(cRed(j),orderRed[j]);
  }
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    cOxi *= pow(cOxi(j),orderOxi[j]);
  }
  return kOxi*exp(aOxi*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cRed
       - kRed*exp(-aRed*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cOxi;
}
//---------------------------------------------------------------------------
double ElecReaction_BV::calcReactionRateDerivativeU(double V) const
{
  double cRed = 1.;
  double cOxi = 1.;
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    cRed *= pow(cRed(j),orderRed[j]);
  }
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    cOxi *= pow(cOxi(j),orderOxi[j]);
  }
  return - kOxi*aOxi*nElectrons*F_CONST/(R_CONST*T)*exp( aOxi*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cRed
       - kRed*aRed*nElectrons*F_CONST/(R_CONST*T)*exp(-aRed*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cOxi;
}
//---------------------------------------------------------------------------
double ElecReaction_BV::calcReactionRateDerivativeV(double V) const
{
  double cRed = 1.;
  double cOxi = 1.;
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    cRed *= pow(cRed(j),orderRed[j]);
  }
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    cOxi *= pow(cOxi(j),orderOxi[j]);
  }
  return   kOxi*aOxi*nElectrons*F_CONST/(R_CONST*T)*exp( aOxi*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cRed
       + kRed*aRed*nElectrons*F_CONST/(R_CONST*T)*exp(-aRed*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cOxi;
}
//---------------------------------------------------------------------------
double ElecReaction_BV::calcReactionRateDerivativeCRed(double V, unsigned i) const
{
  double cRed = 1.;
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    if (j == i) cRed *= pow(cRed(j),orderRed[j]-1);
    else cRed *= pow(cRed(j),orderRed[j]);
  }
  return orderRed[i]*kOxi*exp(aOxi*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cRed;
}
//---------------------------------------------------------------------------
double ElecReaction_BV::calcReactionRateDerivativeCOxi(double V, unsigned i) const
{
  double cOxi = 1.;
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    if (j == i) cOxi *= pow(cOxi(j),orderOxi[j]-1);
    else cOxi *= pow(cOxi(j),orderOxi[j]);
  }
  return -orderOxi[i]*kRed*exp(-aRed*nElectrons*F_CONST/(R_CONST*T)*(V-U))*cOxi;
}
//---------------------------------------------------------------------------
double ElecReaction_BV::calcEquilibriumPotential() const
{
  double cRed = 1.;
  double cOxi = 1.;
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    cRed *= pow(cRed(j),orderRed[j]);
  }
  if (cRed == 0) cRed = 1.;
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    cOxi *= pow(cOxi(j),orderOxi[j]);
  }
  if (cOxi == 0) cOxi = 1.;
  return R_CONST*T/(nElectrons*F_CONST)*log(kRed*cOxi/(kOxi*cRed));
}
//---------------------------------------------------------------------------

