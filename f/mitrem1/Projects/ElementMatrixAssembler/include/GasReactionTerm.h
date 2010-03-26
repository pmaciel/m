//---------------------------------------------------------------------------

#ifndef GasReactionTermH
#define GasReactionTermH

//---------------------------------------------------------------------------

#include "BoundaryElementTerm.h"

//---------------------------------------------------------------------------

class GasReactionTerm : public BoundaryElementTerm
{
public :
  GasReactionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_);
  virtual ~GasReactionTerm();

  virtual void calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions) {};
  virtual void calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions) {};

  virtual void calcElecReactionCurrentDensities(EmptyEmptyDoubleListList boundaryElementElecReactionCurrentDensity, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential) {};
  virtual void calcGasReactionRates(EmptyEmptyDoubleListList boundaryElementGasReactionRate, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions);

  virtual double calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential) {return 0.;};
  virtual double calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions);

protected:
  virtual double calcBubbleReactionRateCorrection(DoubleList surfaceGasFractions);
  double** v;
  double** DvDCDissGas;
  bool isBubble;
};

//---------------------------------------------------------------------------

#endif

