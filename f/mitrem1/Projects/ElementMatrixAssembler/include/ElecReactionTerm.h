//---------------------------------------------------------------------------

#ifndef ElecReactionTermH
#define ElecReactionTermH

//---------------------------------------------------------------------------

#include "BoundaryElementTerm.h"

//---------------------------------------------------------------------------

class ElecReactionTerm : public BoundaryElementTerm
{
public :
	ElecReactionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_);
	virtual ~ElecReactionTerm();

	virtual void	calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions) {};
	virtual void	calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions) {};

	virtual void	calcElecReactionCurrentDensities(EmptyEmptyDoubleListList boundaryElementElecReactionCurrentDensity, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
	virtual void	calcGasReactionRates(EmptyEmptyDoubleListList boundaryElementGasReactionRate, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions) {};

	virtual double	calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
	virtual double	calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions) {return 0.;};

protected:
	double**	v;
	double**	DvDU;
	double***	DvDCRed;
	double***	DvDCOxi;
};

//---------------------------------------------------------------------------

#endif

