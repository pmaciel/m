//---------------------------------------------------------------------------

//#define...

//---------------------------------------------------------------------------

#include "ElecReaction.h"
#include <math.h>
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReaction::ElecReaction (ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_)
: nAgentsRed(nAgentsRed_), nAgentsOxi(nAgentsOxi_), electrolyteSolution(electrolyteSolution_)
{
  agentsRed = new unsigned[nAgentsRed];
  agentsOxi = new unsigned[nAgentsOxi];
  stoichRed = new int[nAgentsRed];
  stoichOxi = new int[nAgentsOxi];
  orderRed = new double[nAgentsRed];
  orderOxi = new double[nAgentsOxi];
}
//---------------------------------------------------------------------------
ElecReaction::~ElecReaction ()
{
  delete[] agentsRed;
  delete[] agentsOxi;
  delete[] stoichRed;
  delete[] stoichOxi;
  delete[] orderRed;
  delete[] orderOxi;
}
//---------------------------------------------------------------------------
int ElecReaction::getStoichOf(unsigned i) const
{
  for (unsigned j=0; j<nAgentsRed; j++)
  {
    if (i == agentsRed[j])
    {
      return stoichRed[j];
    }
  }
  for (unsigned j=0; j<nAgentsOxi; j++)
  {
    if (i == agentsOxi[j])
    {
      return stoichOxi[j];
    }
  }
  return 0;
}
//---------------------------------------------------------------------------
double ElecReaction::calcReactionCurrentDensity(double V) const
{
  return nElectrons*F_CONST*calcReactionRate(V);
}
//---------------------------------------------------------------------------

