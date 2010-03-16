//---------------------------------------------------------------------------

#include "GasReactionTerm_1D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
GasReactionTerm_1D_Galerkin::GasReactionTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_)
	: GasReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_, isBubble_)
{
}
//---------------------------------------------------------------------------
GasReactionTerm_1D_Galerkin::~GasReactionTerm_1D_Galerkin()
{
}
//---------------------------------------------------------------------------
void GasReactionTerm_1D_Galerkin::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	// Calculate coefficients
	mitrem->init(concentrations[0],potentials[0],temperatures[0],densities[0]);
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];
		v[0][s] = mitrem->calcGasReactionRate(s);
	}

	// Add to boundary element vector
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];

		// v(r) gives a contribution its dissolved gas
		unsigned i = mitrem->getGasReactionDissolvedGas(s);
		boundaryElementVec[eq(0,i)] += v[0][s];
	}
}
//---------------------------------------------------------------------------
void GasReactionTerm_1D_Galerkin::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	// Calculate coefficients
	mitrem->init(concentrations[0],potentials[0],temperatures[0],densities[0]);
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];
		DvDCDissGas[0][r] = mitrem->calcGasReactionRateDerivativeCDissGas(s);
	}

	// Add to boundary element jacobian
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];

		// v(s) gave a contribution to its dissolved gas
		unsigned i = mitrem->getGasReactionDissolvedGas(s);
		unsigned eq0i = eq(0,i);

		boundaryElementJac[eq0i][var(0,i)] -= DvDCDissGas[0][r];
	}
}
//---------------------------------------------------------------------------

