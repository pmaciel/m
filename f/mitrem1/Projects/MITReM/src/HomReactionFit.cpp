//---------------------------------------------------------------------------

#include <cstdlib>

#include "HomReactionFit.h"
#include <math.h>
#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionFit::HomReactionFit(const std::string &label_, Electrolyte** electrolytes, unsigned nElectrolytes)
: label(label_)
{
  bool found = false;
  for (unsigned e=0; e<nElectrolytes; e++)
  {
    for (unsigned r=0; r<electrolytes[e]->getNHomReactions(); r++)
    {
      std::string labelTemp = electrolytes[e]->getHomReactionLabel(r);
      if (label == electrolytes[e]->getHomReactionLabel(r))
      {
        found = true;
        homReactionReference.electrolyte = electrolytes[e];
        homReactionReference.electrolyteIndex = e;
        homReactionReference.reactionIndex = r;
        break;
      }
    }
  }
  if (!found) errorHomReactionNotFound(label);
  epsilon = 1e-3;
}
//---------------------------------------------------------------------------
HomReactionFit::~HomReactionFit()
{
}
//---------------------------------------------------------------------------


//--- METHODS ---------------------------------------------------------------
double HomReactionFit::calcOsmoticCoefficientDerivativeEquilibriumConstant(unsigned t)
{
  double K1 = (1.-epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K1);
  double y1 = homReactionReference.electrolyte->calcOsmoticCoefficient(t);
  double K2 = (1.+epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K2);
  double y2 = homReactionReference.electrolyte->calcOsmoticCoefficient(t);

  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K);
  return (y2-y1)/(2.*epsilon*K);
}
//---------------------------------------------------------------------------
double HomReactionFit::calcActivityCoefficientDerivativeEquilibriumConstant(unsigned t)
{
  double K1 = (1.-epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K1);
  double y1 = homReactionReference.electrolyte->calcActivityCoefficient(t);
  double K2 = (1.+epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K2);
  double y2 = homReactionReference.electrolyte->calcActivityCoefficient(t);

  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K);
  return (y2-y1)/(2.*epsilon*K);
}
//---------------------------------------------------------------------------
double HomReactionFit::calcConductivityDerivativeEquilibriumConstant(unsigned t)
{
  double K1 = (1.-epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K1);
  double y1 = homReactionReference.electrolyte->calcConductivity(t);
  double K2 = (1.+epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K2);
  double y2 = homReactionReference.electrolyte->calcConductivity(t);

  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K);
  return (y2-y1)/(2.*epsilon*K);
}
//---------------------------------------------------------------------------
double HomReactionFit::calcTransportNumberDerivativeEquilibriumConstant(unsigned t)
{
  double K1 = (1.-epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K1);
  double y1 = homReactionReference.electrolyte->calcTransportNumber(t);
  double K2 = (1.+epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K2);
  double y2 = homReactionReference.electrolyte->calcTransportNumber(t);

  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K);
  return (y2-y1)/(2.*epsilon*K);
}
//---------------------------------------------------------------------------
double HomReactionFit::calcDiffusionCoefficientDerivativeEquilibriumConstant(unsigned t)
{
  double K1 = (1.-epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K1);
  double y1 = homReactionReference.electrolyte->calcDiffusionCoefficient(t);
  double K2 = (1.+epsilon)*K;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K2);
  double y2 = homReactionReference.electrolyte->calcDiffusionCoefficient(t);

  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K);
  return (y2-y1)/(2.*epsilon*K);
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void HomReactionFit::errorHomReactionNotFound(const std::string label)
{
  std::cout << "ERROR IN HomReactionFit.cpp.\nTHE REACTION " << label
    << " WAS NOWHERE FOUND IN THE *.homreactions FILES." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------

