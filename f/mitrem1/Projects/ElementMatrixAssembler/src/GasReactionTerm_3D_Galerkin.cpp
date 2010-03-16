//---------------------------------------------------------------------------

#include "GasReactionTerm_3D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
GasReactionTerm_3D_Galerkin::GasReactionTerm_3D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_)
	: GasReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_, isBubble_)
{
}
//---------------------------------------------------------------------------
GasReactionTerm_3D_Galerkin::~GasReactionTerm_3D_Galerkin()
{
}
//---------------------------------------------------------------------------
void GasReactionTerm_3D_Galerkin::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++)
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nGasReactions; r++)
		{
			unsigned s = gasReactions[r];
			v[m][s] = mitrem->calcGasReactionRate(s) * calcBubbleReactionRateCorrection(surfaceGasFractions);
		}
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize12 = elementSize/12.;

	// Add to boundary element vector
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];
		for (unsigned m=0; m<nNodes; m++)
		{
			unsigned p = (m+1)%3;
			unsigned q = (m+2)%3;
			double integral = (2.*v[m][s] + v[p][s] + v[q][s])*elementSize12;

			// v(r) gives a contribution its dissolved gas
			unsigned i = mitrem->getGasReactionDissolvedGas(s);
			boundaryElementVec[eq(m,i)] += integral;
		}
	}
}
//---------------------------------------------------------------------------
void GasReactionTerm_3D_Galerkin::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++)
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nGasReactions; r++)
		{
			unsigned s = gasReactions[r];
			DvDCDissGas[m][r] = mitrem->calcGasReactionRateDerivativeCDissGas(s) * calcBubbleReactionRateCorrection(surfaceGasFractions);
		}
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize12 = elementSize/12.;

	// Add to boundary element jacobian
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];

		for (unsigned m=0; m<nNodes; m++)
		{
			unsigned p = (m+1)%3;
			unsigned q = (m+2)%3;

			// v(s) gave a contribution to its dissolved gas
			unsigned i = mitrem->getGasReactionDissolvedGas(s);
			unsigned eqmi = eq(m,i);

			boundaryElementJac[eqmi][var(m,i)] -= 2.*DvDCDissGas[m][r]*elementSize12;
			boundaryElementJac[eqmi][var(p,i)] -=    DvDCDissGas[p][r]*elementSize12;
			boundaryElementJac[eqmi][var(q,i)] -=    DvDCDissGas[q][r]*elementSize12;
		}
	}
}
//---------------------------------------------------------------------------
double GasReactionTerm_3D_Galerkin::calcBubbleReactionRateCorrection(DoubleList surfaceGasFractions)
{
	if (!isBubble)
		return GasReactionTerm::calcBubbleReactionRateCorrection(surfaceGasFractions);

	double surfaceGasFraction = 0.;
	for (unsigned m=0; m<nNodes; ++m)
		surfaceGasFraction += surfaceGasFractions[m];
	surfaceGasFraction /= (double) nNodes;

	return 4. * surfaceGasFraction;
}
//---------------------------------------------------------------------------

