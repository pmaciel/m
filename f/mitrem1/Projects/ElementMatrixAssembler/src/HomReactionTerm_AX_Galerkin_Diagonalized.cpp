//---------------------------------------------------------------------------

#include "HomReactionTerm_AX_Galerkin_Diagonalized.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm_AX_Galerkin_Diagonalized::HomReactionTerm_AX_Galerkin_Diagonalized(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: HomReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
HomReactionTerm_AX_Galerkin_Diagonalized::~HomReactionTerm_AX_Galerkin_Diagonalized()
{
}
//---------------------------------------------------------------------------
void HomReactionTerm_AX_Galerkin_Diagonalized::calcMat(EmptyDoubleMatrix elementMat, EmptyDoubleVector elementVec, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors, double factor)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	
	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nHomReactions; r++) 
		{
			kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
			kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize12 = elementSize/12.;
	
	// Add to element matrix
	for (unsigned r=0; r<nHomReactions; r++) 
	{
		unsigned nReagents = mitrem->getHomReactionNReagents(r);
		unsigned nProducts = mitrem->getHomReactionNProducts(r);
		unsigned* reagents = new unsigned[nReagents];
		unsigned* products = new unsigned[nProducts];
		for (unsigned j=0; j<nReagents; j++)
		{
			reagents[j] = mitrem->getHomReactionReagents(r,j);
		}
		for (unsigned j=0; j<nProducts; j++)
		{
			products[j] = mitrem->getHomReactionProducts(r,j);
		}

		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned n = (m+1)%3;
			unsigned p = (m+2)%3;

			// Forward reaction
			if (nReagents == 1) 
			{
				Hfmj[m] = (2.*kf[m][r] +    kf[n][r] +    kf[p][r])*elementSize12*coordinates[m][0];

				elementMat[eq(m,reagents[0])][var(m,reagents[0])] -= factor*Hfmj[m];
				for (unsigned j=0; j<nProducts; j++)
				{
					elementMat[eq(m,products[j])][var(m,reagents[0])] += factor*Hfmj[m];
				}
			}
			else if (nReagents == 2) 
			{
				Hfmjk[m][m] = ( 2.*kf[m][r] +    kf[n][r] +    kf[p][r])*elementSize12*coordinates[m][0];

				for (unsigned j=0; j<nReagents; j++)
				{
					unsigned k = (j+1)%2;
					Hfmj[m] = 0.5*Hfmjk[m][m]*concentrations[m][reagents[k]];

					elementMat[eq(m,reagents[0])][var(m,reagents[j])] -= factor*Hfmj[m];
					elementMat[eq(m,reagents[1])][var(m,reagents[j])] -= factor*Hfmj[m];
					elementMat[eq(m,products[0])][var(m,reagents[j])] += factor*Hfmj[m];
				}				
			}

			// Backward reaction
			if (nProducts == 1) 
			{
				Hbmj[m] = (2.*kb[m][r] +    kb[n][r] +    kb[p][r])*elementSize12*coordinates[m][0];

				elementMat[eq(m,products[0])][var(m,products[0])] -= factor*Hbmj[m];
				for (unsigned j=0; j<nReagents; j++)
				{
					elementMat[eq(m,reagents[j])][var(m,products[0])] += factor*Hbmj[m];
				}
			}
			else if (nProducts == 2) 
			{
				Hbmjk[m][m] = ( 2.*kb[m][r] +    kb[n][r] +    kb[p][r])*elementSize12*coordinates[m][0];

				for (unsigned j=0; j<nProducts; j++)
				{
					unsigned k = (j+1)%2;
					Hbmj[m] = 0.5*Hbmjk[m][m]*concentrations[m][products[k]];

					elementMat[eq(m,products[0])][var(m,products[j])] -= factor*Hbmj[m];
					elementMat[eq(m,products[1])][var(m,products[j])] -= factor*Hbmj[m];
					elementMat[eq(m,reagents[0])][var(m,products[j])] += factor*Hbmj[m];
				}				
			}
		}
		delete[] products;
		delete[] reagents;
	}
}
//---------------------------------------------------------------------------
void HomReactionTerm_AX_Galerkin_Diagonalized::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors, double factor)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	
	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nHomReactions; r++) 
		{
			kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
			kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize12 = elementSize/12.;
	
	// Add to element matrix
	for (unsigned r=0; r<nHomReactions; r++) 
	{
		unsigned nReagents = mitrem->getHomReactionNReagents(r);
		unsigned nProducts = mitrem->getHomReactionNProducts(r);
		unsigned* reagents = new unsigned[nReagents];
		unsigned* products = new unsigned[nProducts];
		for (unsigned j=0; j<nReagents; j++)
		{
			reagents[j] = mitrem->getHomReactionReagents(r,j);
		}
		for (unsigned j=0; j<nProducts; j++)
		{
			products[j] = mitrem->getHomReactionProducts(r,j);
		}

		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned n = (m+1)%3;
			unsigned p = (m+2)%3;

			// Forward reaction
			if (nReagents == 2) 
			{
				Hfmjk[m][m] = ( 2.*kf[m][r] +    kf[n][r] +    kf[p][r])*elementSize12*coordinates[m][0];

				for (unsigned j=0; j<nReagents; j++)
				{
					unsigned k = (j+1)%2;
					Hfmj[m] = 0.5*Hfmjk[m][m]*concentrations[m][reagents[k]];

					elementJac[eq(m,reagents[0])][var(m,reagents[j])] -= factor*Hfmj[m];
					elementJac[eq(m,reagents[1])][var(m,reagents[j])] -= factor*Hfmj[m];
					elementJac[eq(m,products[0])][var(m,reagents[j])] += factor*Hfmj[m];
				}				
			}

			// Backward reaction
			if (nProducts == 2) 
			{
				Hbmjk[m][m] = ( 2.*kb[m][r] +    kb[n][r] +    kb[p][r])*elementSize12*coordinates[m][0];

				for (unsigned j=0; j<nProducts; j++)
				{
					unsigned k = (j+1)%2;
					Hbmj[m] = 0.5*Hbmjk[m][m]*concentrations[m][products[k]];

					elementJac[eq(m,products[0])][var(m,products[j])] -= factor*Hbmj[m];
					elementJac[eq(m,products[1])][var(m,products[j])] -= factor*Hbmj[m];
					elementJac[eq(m,reagents[0])][var(m,products[j])] += factor*Hbmj[m];
				}				
			}
		}
		delete[] products;
		delete[] reagents;
	}
}
//---------------------------------------------------------------------------

