//---------------------------------------------------------------------------

#ifndef BoundaryElementH
#define BoundaryElementH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "ElecReactionTerm.h"
#include "GasReactionTerm.h"

//---------------------------------------------------------------------------

class BoundaryElementContribution
{
public :
  BoundaryElementContribution(
    unsigned nBoundaryElementNodes_,
    unsigned nIons_,
    unsigned nElecReactionsMax_,
    unsigned nGasReactionsMax_,
    ElecReactionTerm* elecReactionTerm_,
    GasReactionTerm* gasReactionTerm_);
  ~BoundaryElementContribution();

  EmptyDoubleVector  calcVec(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
  EmptyDoubleMatrix  calcJac(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);

  DoubleListList  calcElecReactionCurrentDensities(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
  DoubleListList  calcGasReactionRates(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions);

  double  calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
  double  calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions);

protected:
  ElecReactionTerm*  elecReactionTerm;
  GasReactionTerm*  gasReactionTerm;

  unsigned      nBoundaryElementNodes,nIons,nElecReactionsMax,nGasReactionsMax,size;
  EmptyDoubleVector  boundaryElementVec;
  EmptyDoubleMatrix  boundaryElementJac;
  EmptyEmptyDoubleListList  boundaryElementElecReactionCurrentDensity;
  EmptyEmptyDoubleListList  boundaryElementGasReactionRate;
};

//---------------------------------------------------------------------------

inline double BoundaryElementContribution::calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential)
{
  return elecReactionTerm->calcCurrent(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential);
}

//---------------------------------------------------------------------------

inline double BoundaryElementContribution::calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions)
{
  return gasReactionTerm->calcGasGeneration(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nGasReactions, gasReactions);
}

//---------------------------------------------------------------------------

#endif

