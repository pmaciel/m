//---------------------------------------------------------------------------

#include "HomReactionTerm_1D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm_1D_Galerkin::HomReactionTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: HomReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
HomReactionTerm_1D_Galerkin::~HomReactionTerm_1D_Galerkin()
{
}
//---------------------------------------------------------------------------
void HomReactionTerm_1D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction *= 0.5;
	
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
	double elementSize60 = elementSize/60.;
	
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
			unsigned n = (m+1)%2;

			// Forward reaction
			if (nReagents == 1) 
			{
				Hfmj[m] = (3.*kf[m][r] + kf[n][r])*elementSize12;
				Hfmj[n] = (   kf[m][r] + kf[n][r])*elementSize12;

				elementMat[eq(m,reagents[0])][var(m,reagents[0])] -= Hfmj[m];
				elementMat[eq(m,reagents[0])][var(n,reagents[0])] -= Hfmj[n];
				for (unsigned j=0; j<nProducts; j++)
				{
					elementMat[eq(m,products[j])][var(m,reagents[0])] += Hfmj[m];
					elementMat[eq(m,products[j])][var(n,reagents[0])] += Hfmj[n];
				}
			}
			else if (nReagents == 2) 
			{
				Hfmjk[m][m] = (12.*kf[m][r] + 3.*kf[n][r])*elementSize60;
				Hfmjk[n][n] = ( 2.*kf[m][r] + 3.*kf[n][r])*elementSize60;
				Hfmjk[m][n] = ( 3.*kf[m][r] + 2.*kf[n][r])*elementSize60;
				Hfmjk[n][m] = Hfmjk[m][n];

				for (unsigned j=0; j<nReagents; j++)
				{
					unsigned k = (j+1)%2;
					Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]]);
					Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]]);

					elementMat[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
					elementMat[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
					elementMat[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
					elementMat[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
					for (unsigned h = 0; h < nProducts; h++)
					{
						elementMat[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
						elementMat[eq(m,products[h])][var(n,reagents[j])] += Hfmj[n];
					}
				}				
			}

			// Backward reaction
			if (nProducts == 1) 
			{
				Hbmj[m] = (3.*kb[m][r] + kb[n][r])*elementSize12;
				Hbmj[n] = (   kb[m][r] + kb[n][r])*elementSize12;

				elementMat[eq(m,products[0])][var(m,products[0])] -= Hbmj[m];
				elementMat[eq(m,products[0])][var(n,products[0])] -= Hbmj[n];
				for (unsigned j=0; j<nReagents; j++)
				{
					elementMat[eq(m,reagents[j])][var(m,products[0])] += Hbmj[m];
					elementMat[eq(m,reagents[j])][var(n,products[0])] += Hbmj[n];
				}
			}
			else if (nProducts == 2) 
			{
				Hbmjk[m][m] = (12.*kb[m][r] + 3.*kb[n][r])*elementSize60;
				Hbmjk[n][n] = ( 2.*kb[m][r] + 3.*kb[n][r])*elementSize60;
				Hbmjk[m][n] = ( 3.*kb[m][r] + 2.*kb[n][r])*elementSize60;
				Hbmjk[n][m] = Hbmjk[m][n];

				for (unsigned j=0; j<nProducts; j++)
				{
					unsigned k = (j+1)%2;
					Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]]);
					Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]]);

					elementMat[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
					elementMat[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
					elementMat[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
					elementMat[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
					for (unsigned h = 0; h < nReagents; h++)
					{
						elementMat[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
						elementMat[eq(m,reagents[h])][var(n,products[j])] += Hbmj[n];
					}
				}				
			}
		}
		delete[] products;
		delete[] reagents;
	}
}
//---------------------------------------------------------------------------
void HomReactionTerm_1D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction *= 0.5;
	
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
	double elementSize60 = elementSize/60.;
	
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
			unsigned n = (m+1)%2;

			// Forward reaction
			if (nReagents == 2) 
			{
				Hfmjk[m][m] = (12.*kf[m][r] + 3.*kf[n][r])*elementSize60;
				Hfmjk[n][n] = ( 2.*kf[m][r] + 3.*kf[n][r])*elementSize60;
				Hfmjk[m][n] = ( 3.*kf[m][r] + 2.*kf[n][r])*elementSize60;
				Hfmjk[n][m] = Hfmjk[m][n];

				for (unsigned j=0; j<nReagents; j++)
				{
					unsigned k = (j+1)%2;
					Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]]);
					Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]]);

					elementJac[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
					elementJac[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
					elementJac[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
					elementJac[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
					for (unsigned h = 0; h < nProducts; h++)
					{
						elementJac[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
						elementJac[eq(m,products[h])][var(n,reagents[j])] += Hfmj[n];
					}
				}				
			}

			// Backward reaction
			if (nProducts == 2) 
			{
				Hbmjk[m][m] = (12.*kb[m][r] + 3.*kb[n][r])*elementSize60;
				Hbmjk[n][n] = ( 2.*kb[m][r] + 3.*kb[n][r])*elementSize60;
				Hbmjk[m][n] = ( 3.*kb[m][r] + 2.*kb[n][r])*elementSize60;
				Hbmjk[n][m] = Hbmjk[m][n];

				for (unsigned j=0; j<nProducts; j++)
				{
					unsigned k = (j+1)%2;
					Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]]);
					Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]]);

					elementJac[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
					elementJac[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
					elementJac[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
					elementJac[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
					for (unsigned h = 0; h < nReagents; h++)
					{
						elementJac[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
						elementJac[eq(m,reagents[h])][var(n,products[j])] += Hbmj[n];
					}
				}				
			}
		}
		delete[] products;
		delete[] reagents;
	}
}
//---------------------------------------------------------------------------

