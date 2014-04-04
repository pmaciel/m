//---------------------------------------------------------------------------

#include "HomReaction.h"
#include <math.h>
#include <iostream>

//---------------------------------------------------------------------------

#define cReag(j)    electrolyteSolution->getIonConcentration(reagents[j])
#define cProd(j)    electrolyteSolution->getIonConcentration(products[j])
#define fReag(j)    electrolyteModel->calcActivityCoefficient_MM(reagents[j])
#define fProd(j)    electrolyteModel->calcActivityCoefficient_MM(products[j])

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReaction::HomReaction(ElectrolyteSolution* electrolyteSolution_, ElectrolyteModel* electrolyteModel_, unsigned nReagents_, unsigned nProducts_)
: nReagents(nReagents_), nProducts(nProducts_), electrolyteSolution(electrolyteSolution_), electrolyteModel(electrolyteModel_)
{
  reagents = new unsigned[nReagents];
  products = new unsigned[nProducts];
  stoichReag = new int[nReagents];
  stoichProd = new int[nProducts];
  cReagentsSave = new double[nReagents];
  cProductsSave = new double[nProducts];
}
//---------------------------------------------------------------------------
HomReaction::~HomReaction ()
{
  delete[] reagents;
  delete[] products;
  delete[] stoichReag;
  delete[] stoichProd;
  delete[] cReagentsSave;
  delete[] cProductsSave;
}
//---------------------------------------------------------------------------


//--- METHODS ---------------------------------------------------------------
double HomReaction::calcForwardRateConstant() const
{
  double fReag = 1.;
  for (unsigned j=0; j<nReagents; j++)
  {
    fReag *= pow(fReag(j),-stoichReag[j]);
  }
  return kf*fReag;
}
//---------------------------------------------------------------------------
double HomReaction::calcBackwardRateConstant() const
{
  double fProd = 1.;
  for (unsigned j=0; j<nProducts; j++)
  {
    fProd *= pow(fProd(j),stoichProd[j]);
  }
  return kb*fProd;
}
//---------------------------------------------------------------------------
double HomReaction::calcReactionRate() const
{
  double cReag = 1.;
  double cProd = 1.;
  for (unsigned j=0; j<nReagents; j++) {
    cReag *= pow(cReag(j),-stoichReag[j]);
  }
  for (unsigned j=0; j<nProducts; j++)
  {
    cProd *= pow(cProd(j),stoichProd[j]);
  }
  return calcForwardRateConstant()*cReag - calcBackwardRateConstant()*cProd;
}
//---------------------------------------------------------------------------
double HomReaction::calcMaximumBackwardProgress() const
{
  double progressMax = -cProductsSave[0]/stoichProd[0];
  for (unsigned i=1; i<nProducts; i++)
  {
    double progressMaxTemp = -cProductsSave[i]/stoichProd[i];
    if (progressMaxTemp > progressMax) progressMax = progressMaxTemp;
  }
  return progressMax; //this is a negative number
}
//---------------------------------------------------------------------------
double HomReaction::calcMaximumForwardProgress() const
{
  double progressMax = -cReagentsSave[0]/stoichReag[0];
  for (unsigned i=1; i<nReagents; i++)
  {
    double progressMaxTemp = -cReagentsSave[i]/stoichReag[i];
    if (progressMaxTemp < progressMax) progressMax = progressMaxTemp;
  }
  return progressMax; //this is a positive number
}
//---------------------------------------------------------------------------
double HomReaction::progressToEquilibrium() const
{
  for (unsigned i=0; i<nReagents; i++)
  {
    cReagentsSave[i] = cReag(i);
  }
  for (unsigned i=0; i<nProducts; i++)
  {
    cProductsSave[i] = cProd(i);
  }
  double deviationFromEquilibrium = calcDeviationFromEquilibrium();
  double x = 0.;
  if (deviationFromEquilibrium < 0.)
  {
    double xMax = calcMaximumForwardProgress();
    x = calcEquilibriumProgress(0.,xMax);
    progress(x);
  }
  else if (deviationFromEquilibrium > 0.)
  {
    double xMin = calcMaximumBackwardProgress();
    x = calcEquilibriumProgress(xMin,0.);
    progress(x);
  }
  return x;
}
//---------------------------------------------------------------------------
double HomReaction::calcEquilibriumProgress(double xMin, double xMax) const
{
  double x = 0.5*(xMin+xMax);
  double xPrevious;
  double dx,deviationFromEquilibrium;
  unsigned iter = 0;
  const unsigned iterMax = 20;
  do
  {
    xPrevious = x;
    progress(xPrevious);
    deviationFromEquilibrium = calcDeviationFromEquilibrium();
    dx = 0.;
    for (unsigned i=0; i<nProducts; i++)
    {
      dx += stoichProd[i]*stoichProd[i]/cProd(i);
    }
    for (unsigned i=0; i<nReagents; i++)
    {
      dx += stoichReag[i]*stoichReag[i]/cReag(i);
    }
    dx = -deviationFromEquilibrium/((deviationFromEquilibrium+K)*dx);
    x = xPrevious + dx;
    if (x >= xMax) x = 0.5*(xPrevious + xMax);
    else if (x <= xMin) x = 0.5*(xPrevious + xMin);
    iter++;
  } while ((fabs(dx) > 0.01*fabs(xPrevious)) && (iter < iterMax));
  return x;
}
//---------------------------------------------------------------------------
void HomReaction::progress(double x) const
{
  for (unsigned i=0; i<nReagents; i++)
  {
    electrolyteSolution->setIonConcentration(reagents[i],cReagentsSave[i]+stoichReag[i]*x);
  }
  for (unsigned i=0; i<nProducts; i++)
  {
    electrolyteSolution->setIonConcentration(products[i],cProductsSave[i]+stoichProd[i]*x);
  }
}
//---------------------------------------------------------------------------
double HomReaction::calcActivityProduct() const
{
  double activityProduct = 1.;
  for (unsigned i=0; i<nProducts; i++)
  {
    double activity = fProd(i)*cProd(i);
    activityProduct *= pow(activity,stoichProd[i]);
  }
  for (unsigned i=0; i<nReagents; i++)
  {
    double activity = fReag(i)*cReag(i);
    activityProduct *= pow(activity,stoichReag[i]);
  }
  return activityProduct;
}
//---------------------------------------------------------------------------
double HomReaction::calcDeviationFromEquilibrium() const
{
  return calcActivityProduct() - K;
}
//---------------------------------------------------------------------------
double HomReaction::calcRelativeDeviationFromEquilibrium() const
{
  return log(calcActivityProduct()/K);
}
//---------------------------------------------------------------------------
double HomReaction::calcDeviationFromEquilibriumDerivative(int* stoichList) const
{
  // This function calculates the derivative of this reaction's deviation from
  // equilibrium with respect to the progress of another reaction (df[r]/dx[s]).
  // For this you need both the stoichiometric coefficients of the ions in this reaction
  // (stoichReag and stoichProd) and in the other reaction (signedStoichList)!
  // It is not an exact derivative because we do not take into account the effect that
  // the progresses have on the activity coefficients. But this is not a big problem,
  // because this derivative serves just for the Newton method.

  double derivative = 0.;
  for (unsigned i=0; i<nProducts; i++)
  {
    derivative += stoichProd[i]*stoichList[products[i]]/cProd(i);
  }
  for (unsigned i=0; i<nReagents; i++)
  {
    derivative += stoichReag[i]*stoichList[reagents[i]]/cReag(i);
  }
  derivative *= calcActivityProduct();
  return derivative;
}
//---------------------------------------------------------------------------
int HomReaction::getStoichOf(unsigned i) const
{
  for (unsigned j=0; j<nReagents; j++)
  {
    if (i == reagents[j])
    {
      return stoichReag[j];
    }
  }
  for (unsigned j=0; j<nProducts; j++)
  {
    if (i == products[j])
    {
      return stoichProd[j];
    }
  }
  return 0;
}
//---------------------------------------------------------------------------

