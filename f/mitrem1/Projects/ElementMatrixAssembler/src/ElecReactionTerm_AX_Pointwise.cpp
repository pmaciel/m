//---------------------------------------------------------------------------

#include "ElecReactionTerm_AX_Pointwise.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReactionTerm_AX_Pointwise::ElecReactionTerm_AX_Pointwise(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_) 
	: ElecReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_)
{
}
//---------------------------------------------------------------------------
ElecReactionTerm_AX_Pointwise::~ElecReactionTerm_AX_Pointwise()
{
}
//---------------------------------------------------------------------------
void ElecReactionTerm_AX_Pointwise::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions, double factor)
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
		//unsigned p = (m+1)%2;
		R[m][m] = coordinates[m][0];
/*		R[m][p] =    coordinates[m][0] + coordinates[p][0];*/
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize2 = elementSize/2.;
	
	// Add to boundary element vector
	for (unsigned r=0; r<nElecReactions; r++) 
	{
		unsigned s = elecReactions[r];
		unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
		unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);
		for (unsigned m=0; m<nNodes; m++) 
		{
//			unsigned p = (m+1)%2;
			double integral = (R[m][m]*v[m][s])*elementSize2;

			// v(s) gives a contribution to all its RedAgents
			for (unsigned j=0; j<nAgentsRed; j++) 
			{
				unsigned i = mitrem->getElecReactionAgentsRed(s,j);
				int stoichRed = mitrem->getElecReactionStoichRed(s,j);

				boundaryElementVec[eq(m,i)] -= factor*integral*stoichRed;
			}

			// v(s) gives a contribution to all its OxiAgents
			for (unsigned j=0; j<nAgentsOxi; j++) 
			{
				unsigned i = mitrem->getElecReactionAgentsOxi(s,j);
				int stoichOxi = mitrem->getElecReactionStoichOxi(s,j);

				boundaryElementVec[eq(m,i)] -= factor*integral*stoichOxi;
			}

		}
	}
}
//---------------------------------------------------------------------------
void ElecReactionTerm_AX_Pointwise::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions, double factor)
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
		//unsigned p = (m+1)%2;
		R[m][m] = coordinates[m][0];
//		R[m][p] =    coordinates[m][0] + coordinates[p][0];
	}
	elementSize = boundaryElementProps->calcSize(coordinates);
	double elementSize2 = elementSize/2.;

	// Add to boundary element jacobian
	for (unsigned r=0; r<nElecReactions; r++) 
	{
		unsigned s = elecReactions[r];
		unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
		unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);

		for (unsigned m=0; m<nNodes; m++) 
		{
			//unsigned p = (m+1)%2;

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
					boundaryElementJac[eqmi][var(m,j)] += factor*R[m][m]*DvDCRed[m][r][l]*elementSize2*stoichRed;
//					boundaryElementJac[eqmi][var(p,j)] += factor*R[m][p]*DvDCRed[p][r][l]*elementSize12*stoichRed;
				}

				// derivatives of v(s) by its agentsOxi
				for (unsigned l=0; l<nAgentsOxi; l++) 
				{
					unsigned j = mitrem->getElecReactionAgentsOxi(s,l);
					boundaryElementJac[eqmi][var(m,j)] += factor*R[m][m]*DvDCOxi[m][r][l]*elementSize2*stoichRed;
//					boundaryElementJac[eqmi][var(p,j)] += factor*R[m][p]*DvDCOxi[p][r][l]*elementSize12*stoichRed;
				}

				// derivative of v(s) by U
				boundaryElementJac[eqmi][var(m,nIons)] += factor*R[m][m]*DvDU[m][r]*elementSize2*stoichRed;
//				boundaryElementJac[eqmi][var(p,nIons)] += factor*R[m][p]*DvDU[p][r]*elementSize12*stoichRed;
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
					boundaryElementJac[eqmi][var(m,j)] += factor*R[m][m]*DvDCRed[m][r][l]*elementSize2*stoichOxi;
//					boundaryElementJac[eqmi][var(p,j)] += factor*R[m][p]*DvDCRed[p][r][l]*elementSize12*stoichOxi;
				}

				// derivatives of v(s) by its agentsOxi
				for (unsigned l=0; l<nAgentsOxi; l++) 
				{
					unsigned j = mitrem->getElecReactionAgentsOxi(s,l);
					boundaryElementJac[eqmi][var(m,j)] += factor*R[m][m]*DvDCOxi[m][r][l]*elementSize2*stoichOxi;
//					boundaryElementJac[eqmi][var(p,j)] += factor*R[m][p]*DvDCOxi[p][r][l]*elementSize12*stoichOxi;
				}

				// derivative of v(s) by U
				boundaryElementJac[eqmi][var(m,nIons)] += factor*R[m][m]*DvDU[m][r]*elementSize2*stoichOxi;
//				boundaryElementJac[eqmi][var(p,nIons)] += factor*R[m][p]*DvDU[p][r]*elementSize12*stoichOxi;
			}
		}
	}
}
//---------------------------------------------------------------------------
double ElecReactionTerm_AX_Pointwise::calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential)
{
	double surfaceGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		surfaceGasFraction += surfaceGasFractions[m];
	}
	surfaceGasFraction /= nNodes;
	
	double current = 0.;
	elementSize = boundaryElementProps->calcSize(coordinates);
	
	for (unsigned r=0; r<nElecReactions; r++)
	{
		unsigned s = elecReactions[r];
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned p = (m+1)%2;
			
			// Calculate coefficients
			mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
			v[m][s] = mitrem->calcElecReactionCurrentDensity(s,electrodePotential)*(1.-surfaceGasFraction);

			// Add to current
			current += v[m][s]*(coordinates[m][0]/3. + coordinates[p][0]/6.);
		}
	}
	
	// Total current
	current *= 2.*M_PI*elementSize;

	return current;
}
//---------------------------------------------------------------------------
