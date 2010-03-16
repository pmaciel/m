//---------------------------------------------------------------------------

#include "GasReactionTerm_AX_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
GasReactionTerm_AX_Galerkin::GasReactionTerm_AX_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_)
	: GasReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_, isBubble_)
{
}
//---------------------------------------------------------------------------
GasReactionTerm_AX_Galerkin::~GasReactionTerm_AX_Galerkin()
{
}
//---------------------------------------------------------------------------
void GasReactionTerm_AX_Galerkin::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	// Calculate coefficients
	double R[2][2];
	for (unsigned m=0; m<nNodes; m++)
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nGasReactions; r++)
		{
			unsigned s = gasReactions[r];
			v[m][s] = mitrem->calcGasReactionRate(s);
		}
		unsigned p = (m+1)%2;
		R[m][m] = 3.*coordinates[m][0] + coordinates[p][0];
		R[m][p] =    coordinates[m][0] + coordinates[p][0];
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize12 = elementSize/12.;

	// Add to boundary element vector
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];
		for (unsigned m=0; m<nNodes; m++)
		{
			unsigned p = (m+1)%2;
			double integral = (R[m][m]*v[m][s] + R[m][p]*v[p][s])*elementSize12;

			// v(r) gives a contribution its dissolved gas
			unsigned i = mitrem->getGasReactionDissolvedGas(s);
			boundaryElementVec[eq(m,i)] += integral;
		}
	}
}
//---------------------------------------------------------------------------
void GasReactionTerm_AX_Galerkin::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	// Calculate coefficients
	double R[2][2];
	for (unsigned m=0; m<nNodes; m++)
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nGasReactions; r++)
		{
			unsigned s = gasReactions[r];
			DvDCDissGas[m][r] = mitrem->calcGasReactionRateDerivativeCDissGas(s);
		}
		unsigned p = (m+1)%2;
		R[m][m] = 3.*coordinates[m][0] + coordinates[p][0];
		R[m][p] =    coordinates[m][0] + coordinates[p][0];
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize12 = elementSize/12.;

	// Add to boundary element jacobian
	for (unsigned r=0; r<nGasReactions; r++)
	{
		unsigned s = gasReactions[r];

		for (unsigned m=0; m<nNodes; m++)
		{
			unsigned p = (m+1)%2;

			// v(s) gave a contribution to its dissolved gas
			unsigned i = mitrem->getGasReactionDissolvedGas(s);
			unsigned eqmi = eq(m,i);

			boundaryElementJac[eqmi][var(m,i)] -= R[m][m]*DvDCDissGas[m][r]*elementSize12;
			boundaryElementJac[eqmi][var(p,i)] -= R[m][p]*DvDCDissGas[p][r]*elementSize12;
		}
	}
}
//---------------------------------------------------------------------------

