//---------------------------------------------------------------------------

#define Tmp      electrolyteSolution->getSolutionTemperature()
#define cRed(j)  electrolyteSolution->getIonConcentration(agentsRed[j])
#define cOxi(j)  electrolyteSolution->getIonConcentration(agentsOxi[j])
#define cAds     electrolyteSolution->getIonConcentration(electrolyteSolution->getNIons()-1)
#define U        electrolyteSolution->getSolutionPotential()
#define kOxi     kineticParameters[0]
#define kRed     kineticParameters[1]
#define aOxi     kineticParameters[2]
#define aRed     kineticParameters[3]
#define KAds     kineticParameters[4]
#define aAds     kineticParameters[5]

//---------------------------------------------------------------------------

#include "ElecReaction_BVads.h"
#include <math.h>
#include "MathematicsPhysics.h"
#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReaction_BVads::ElecReaction_BVads(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_)
: ElecReaction(electrolyteSolution_, nAgentsRed_, nAgentsOxi_)
{
  kineticParameters = new double[6];
}
//---------------------------------------------------------------------------
ElecReaction_BVads::~ElecReaction_BVads()
{
  delete[] kineticParameters;
}
//---------------------------------------------------------------------------
double ElecReaction_BVads::calcReactionRate(double V) const
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
  double thetaFree = 1./(1.+KAds*exp(aAds*F_CONST/(R_CONST*Tmp)*(V-U))*cAds);
  return (kOxi*exp(aOxi*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cRed
       - kRed*exp(-aRed*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cOxi) * thetaFree;
}
//---------------------------------------------------------------------------
double ElecReaction_BVads::calcReactionRateDerivativeU(double V) const
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
  double thetaFree = 1./(1.+KAds*exp(aAds*F_CONST/(R_CONST*Tmp)*(V-U))*cAds);
  return (- kOxi*aOxi*nElectrons*F_CONST/(R_CONST*Tmp)*exp( aOxi*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cRed
          - kRed*aRed*nElectrons*F_CONST/(R_CONST*Tmp)*exp(-aRed*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cOxi) * thetaFree;
}
//---------------------------------------------------------------------------
double ElecReaction_BVads::calcReactionRateDerivativeV(double V) const
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
  double thetaFree = 1./(1.+KAds*exp(aAds*F_CONST/(R_CONST*Tmp)*(V-U))*cAds);
  return (kOxi*aOxi*nElectrons*F_CONST/(R_CONST*Tmp)*exp( aOxi*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cRed
        + kRed*aRed*nElectrons*F_CONST/(R_CONST*Tmp)*exp(-aRed*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cOxi) * thetaFree;
}
//---------------------------------------------------------------------------
double ElecReaction_BVads::calcReactionRateDerivativeCRed(double V, unsigned i) const
{
  double cRed = 1.;
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    if (j == i) cRed *= pow(cRed(j),orderRed[j]-1);
    else cRed *= pow(cRed(j),orderRed[j]);
  }
  double thetaFree = 1./(1.+KAds*exp(aAds*F_CONST/(R_CONST*Tmp)*(V-U))*cAds);
  return orderRed[i]*kOxi*exp(aOxi*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cRed * thetaFree;
}
//---------------------------------------------------------------------------
double ElecReaction_BVads::calcReactionRateDerivativeCOxi(double V, unsigned i) const
{
  double cOxi = 1.;
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    if (j == i) cOxi *= pow(cOxi(j),orderOxi[j]-1);
    else cOxi *= pow(cOxi(j),orderOxi[j]);
  }
  double thetaFree = 1./(1.+KAds*exp(aAds*F_CONST/(R_CONST*Tmp)*(V-U))*cAds);
  return -orderOxi[i]*kRed*exp(-aRed*nElectrons*F_CONST/(R_CONST*Tmp)*(V-U))*cOxi * thetaFree;
}
//---------------------------------------------------------------------------
double ElecReaction_BVads::calcEquilibriumPotential() const
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
  return R_CONST*Tmp/(nElectrons*F_CONST)*log(kRed*cOxi/(kOxi*cRed));
}
//---------------------------------------------------------------------------

