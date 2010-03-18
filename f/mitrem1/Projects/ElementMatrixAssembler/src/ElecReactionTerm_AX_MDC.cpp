//---------------------------------------------------------------------------

#include "ElecReactionTerm_AX_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReactionTerm_AX_MDC::ElecReactionTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_) 
	: ElecReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_)
{
}
//---------------------------------------------------------------------------
ElecReactionTerm_AX_MDC::~ElecReactionTerm_AX_MDC()
{
}
//---------------------------------------------------------------------------
void ElecReactionTerm_AX_MDC::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	double surfaceGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		surfaceGasFraction += surfaceGasFractions[m];
	}
	surfaceGasFraction *= 0.5;
	
	// Calculate coefficients
	double R[2][2];
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nElecReactions; r++) 
		{
			unsigned s = elecReactions[r];
			v[m][s] = mitrem->calcElecReactionRate(s,electrodePotential)*(1.-surfaceGasFraction);
		}
		unsigned p = (m+1)%2;
		R[m][m] = 7.*coordinates[m][0] + 2.*coordinates[p][0];
		R[m][p] = 2.*coordinates[m][0] +    coordinates[p][0];
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize24 = elementSize/24.;
	
	// Add to boundary element vector
	for (unsigned r=0; r<nElecReactions; r++) 
	{
		unsigned s = elecReactions[r];
		unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
		unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned p = (m+1)%2;
			double integral = (R[m][m]*v[m][s] + R[m][p]*v[p][s])*elementSize24;

			// v(s) gives a contribution to all its RedAgents
			for (unsigned j=0; j<nAgentsRed; j++) 
			{
				unsigned i = mitrem->getElecReactionAgentsRed(s,j);
				int stoichRed = mitrem->getElecReactionStoichRed(s,j);

				boundaryElementVec[eq(m,i)] -= integral*stoichRed;
			}

			// v(s) gives a contribution to all its OxiAgents
			for (unsigned j=0; j<nAgentsOxi; j++) 
			{
				unsigned i = mitrem->getElecReactionAgentsOxi(s,j);
				int stoichOxi = mitrem->getElecReactionStoichOxi(s,j);

				boundaryElementVec[eq(m,i)] -= integral*stoichOxi;
			}

		}
	}
}
//---------------------------------------------------------------------------
void ElecReactionTerm_AX_MDC::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	double surfaceGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		surfaceGasFraction += surfaceGasFractions[m];
	}
	surfaceGasFraction *= 0.5;
	
	// Calculate coefficients
	double R[2][2];
	for (unsigned m=0; m<nNodes; m++)
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nElecReactions; r++) 
		{
			unsigned s = elecReactions[r];
			unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
			unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);
			
			DvDU[m][r] = mitrem->calcElecReactionRateDerivativeU(s,electrodePotential)*(1.-surfaceGasFraction);
			for (unsigned j=0; j<nAgentsRed; j++) 
			{
				DvDCRed[m][r][j] = mitrem->calcElecReactionRateDerivativeCRed(s,electrodePotential,j)*(1.-surfaceGasFraction);
			}
			for (unsigned j=0; j<nAgentsOxi; j++) 
			{
				DvDCOxi[m][r][j] = mitrem->calcElecReactionRateDerivativeCOxi(s,electrodePotential,j)*(1.-surfaceGasFraction);
			}
		}
		unsigned p = (m+1)%2;
		R[m][m] = 7.*coordinates[m][0] + 2.*coordinates[p][0];
		R[m][p] = 2.*coordinates[m][0] +    coordinates[p][0];
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize24 = elementSize/24.;

	// Add to boundary element jacobian
	for (unsigned r=0; r<nElecReactions; r++) 
	{
		unsigned s = elecReactions[r];
		unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
		unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);

		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned p = (m+1)%2;

			// v(s) gave a contribution to all its agentsRed
			for (unsigned k=0; k<nAgentsRed; k++) 
			{
				unsigned i = mitrem->getElecReactionAgentsRed(s,k);
				unsigned eqmi = eq(m,i);
				int stoichRed = mitrem->getElecReactionStoichRed(s,k);

				// derivatives of v(s) by its agentsRed
				for (unsigned l=0; l<nAgentsRed; l++) 
				{
					unsigned j = mitrem->getElecReactionAgentsRed(s,l);
					boundaryElementJac[eqmi][var(m,j)] += R[m][m]*DvDCRed[m][r][l]*elementSize24*stoichRed;
					boundaryElementJac[eqmi][var(p,j)] += R[m][p]*DvDCRed[p][r][l]*elementSize24*stoichRed;
				}

				// derivatives of v(s) by its agentsOxi
				for (unsigned l=0; l<nAgentsOxi; l++) 
				{
					unsigned j = mitrem->getElecReactionAgentsOxi(s,l);
					boundaryElementJac[eqmi][var(m,j)] += R[m][m]*DvDCOxi[m][r][l]*elementSize24*stoichRed;
					boundaryElementJac[eqmi][var(p,j)] += R[m][p]*DvDCOxi[p][r][l]*elementSize24*stoichRed;
				}

				// derivative of v(s) by U
				boundaryElementJac[eqmi][var(m,nIons)] += R[m][m]*DvDU[m][r]*elementSize24*stoichRed;
				boundaryElementJac[eqmi][var(p,nIons)] += R[m][p]*DvDU[p][r]*elementSize24*stoichRed;
			}

			// v(s) gave a contribution to all its agentsOxi
			for (unsigned k=0; k<nAgentsOxi; k++) 
			{
				unsigned i = mitrem->getElecReactionAgentsOxi(s,k);
				unsigned eqmi = eq(m,i);
				int stoichOxi = mitrem->getElecReactionStoichOxi(s,k);

				// derivatives of v(s) by its agentsRed
				for (unsigned l=0; l<nAgentsRed; l++) 
				{
					unsigned j = mitrem->getElecReactionAgentsRed(s,l);
					boundaryElementJac[eqmi][var(m,j)] += R[m][m]*DvDCRed[m][r][l]*elementSize24*stoichOxi;
					boundaryElementJac[eqmi][var(p,j)] += R[m][p]*DvDCRed[p][r][l]*elementSize24*stoichOxi;
				}

				// derivatives of v(s) by its agentsOxi
				for (unsigned l=0; l<nAgentsOxi; l++) 
				{
					unsigned j = mitrem->getElecReactionAgentsOxi(s,l);
					boundaryElementJac[eqmi][var(m,j)] += R[m][m]*DvDCOxi[m][r][l]*elementSize24*stoichOxi;
					boundaryElementJac[eqmi][var(p,j)] += R[m][p]*DvDCOxi[p][r][l]*elementSize24*stoichOxi;
				}

				// derivative of v(s) by U
				boundaryElementJac[eqmi][var(m,nIons)] += R[m][m]*DvDU[m][r]*elementSize24*stoichOxi;
				boundaryElementJac[eqmi][var(p,nIons)] += R[m][p]*DvDU[p][r]*elementSize24*stoichOxi;
			}
		}
	}
}
//---------------------------------------------------------------------------

