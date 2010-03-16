//---------------------------------------------------------------------------

#ifndef GasReactionTerm_AX_GalerkinH
#define GasReactionTerm_AX_GalerkinH

//---------------------------------------------------------------------------

#include "GasReactionTerm.h"

//---------------------------------------------------------------------------

class GasReactionTerm_AX_Galerkin : public GasReactionTerm
{
public :
	GasReactionTerm_AX_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_);
	virtual ~GasReactionTerm_AX_Galerkin();

	virtual void calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
	virtual void calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
};

//---------------------------------------------------------------------------

#endif

