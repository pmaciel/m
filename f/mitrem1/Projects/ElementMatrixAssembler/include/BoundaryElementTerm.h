//---------------------------------------------------------------------------

#ifndef BoundaryElementTermH
#define BoundaryElementTermH

//---------------------------------------------------------------------------

#include "Term.h"
#include "BoundaryElementProps.h"

//---------------------------------------------------------------------------

class BoundaryElementTerm : public Term
{
public :
  BoundaryElementTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_);
  virtual ~BoundaryElementTerm();

  virtual void  calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions) = 0;
  virtual void  calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions) = 0;

  virtual void  calcElecReactionCurrentDensities(EmptyEmptyDoubleListList boundaryElementElecReactionCurrentDensity, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential) = 0;
  virtual void  calcGasReactionRates(EmptyEmptyDoubleListList boundaryElementGasReactionRate, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions) = 0;

  virtual double  calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential) = 0;
  virtual double  calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions) = 0;

protected :
  BoundaryElementProps*  boundaryElementProps;
};

//---------------------------------------------------------------------------

#endif

