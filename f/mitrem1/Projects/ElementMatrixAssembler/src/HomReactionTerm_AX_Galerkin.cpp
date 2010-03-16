//---------------------------------------------------------------------------

#include "HomReactionTerm_AX_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm_AX_Galerkin::HomReactionTerm_AX_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: HomReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
HomReactionTerm_AX_Galerkin::~HomReactionTerm_AX_Galerkin()
{
}
//---------------------------------------------------------------------------
void HomReactionTerm_AX_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	
	// Calculate coefficients
	double R1[3][3][3];
	double R2[3][3][3][3];
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nHomReactions; r++) 
		{
			kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
			kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
		}
		unsigned n = (m+1)%3;
		unsigned p = (m+2)%3;
		R1[m][m][m] = 12.*coordinates[m][0] + 3.*coordinates[n][0] + 3.*coordinates[p][0];
		R1[m][n][n] =  2.*coordinates[m][0] + 3.*coordinates[n][0] +    coordinates[p][0];
		R1[m][p][p] =  2.*coordinates[m][0] +    coordinates[n][0] + 3.*coordinates[p][0];
		R1[m][m][n] =  3.*coordinates[m][0] + 2.*coordinates[n][0] +    coordinates[p][0];
		R1[m][n][m] = R1[m][m][n];
		R1[m][m][p] =  3.*coordinates[m][0] +    coordinates[n][0] + 2.*coordinates[p][0];
		R1[m][p][m] = R1[m][m][p];
		R1[m][n][p] =     coordinates[m][0] +    coordinates[n][0] +    coordinates[p][0];
		R1[m][p][n] = R1[m][n][p];
		R2[m][m][m][m] = 60.*coordinates[m][0] + 12.*coordinates[n][0] + 12.*coordinates[p][0];
		R2[m][n][n][n] =  6.*coordinates[m][0] + 12.*coordinates[n][0] +  3.*coordinates[p][0];
		R2[m][p][p][p] =  6.*coordinates[m][0] +  3.*coordinates[n][0] + 12.*coordinates[p][0];
		R2[m][m][m][n] = 12.*coordinates[m][0] +  6.*coordinates[n][0] +  3.*coordinates[p][0];
		R2[m][m][n][m] = R2[m][m][m][n];
		R2[m][n][m][m] = R2[m][m][m][n];
		R2[m][m][m][p] = 12.*coordinates[m][0] +  3.*coordinates[n][0] +  6.*coordinates[p][0];
		R2[m][m][p][m] = R2[m][m][m][p];
		R2[m][p][m][m] = R2[m][m][m][p];
		R2[m][m][n][n] =  6.*coordinates[m][0] +  6.*coordinates[n][0] +  2.*coordinates[p][0];
		R2[m][n][m][n] = R2[m][m][n][n];
		R2[m][n][n][m] = R2[m][m][n][n];
		R2[m][m][p][p] =  6.*coordinates[m][0] +  2.*coordinates[n][0] +  6.*coordinates[p][0];
		R2[m][p][m][p] = R2[m][m][p][p];
		R2[m][p][p][m] = R2[m][m][p][p];
		R2[m][m][n][p] =  3.*coordinates[m][0] +  2.*coordinates[n][0] +  2.*coordinates[p][0];
		R2[m][m][p][n] = R2[m][m][n][p];
		R2[m][n][m][p] = R2[m][m][n][p];
		R2[m][p][m][n] = R2[m][m][n][p];
		R2[m][n][p][m] = R2[m][m][n][p];
		R2[m][p][n][m] = R2[m][m][n][p];
		R2[m][n][n][p] =  2.*coordinates[m][0] +  3.*coordinates[n][0] +  2.*coordinates[p][0];
		R2[m][n][p][n] = R2[m][n][n][p];
		R2[m][p][n][n] = R2[m][n][n][p];
		R2[m][p][p][n] =  2.*coordinates[m][0] +  2.*coordinates[n][0] +  3.*coordinates[p][0];
		R2[m][p][n][p] = R2[m][p][p][n];
		R2[m][n][p][p] = R2[m][p][p][n];
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize180 = elementSize/180.;
	//double elementSize2520 = elementSize/2520.; BUG?
	double elementSize1260 = elementSize/1260;
	
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
				Hfmj[m] = (R1[m][m][m]*kf[m][r] + R1[m][m][n]*kf[n][r] + R1[m][m][p]*kf[p][r])*elementSize180;
				Hfmj[n] = (R1[m][n][m]*kf[m][r] + R1[m][n][n]*kf[n][r] + R1[m][n][p]*kf[p][r])*elementSize180;
				Hfmj[p] = (R1[m][p][m]*kf[m][r] + R1[m][p][n]*kf[n][r] + R1[m][p][p]*kf[p][r])*elementSize180;

				elementMat[eq(m,reagents[0])][var(m,reagents[0])] -= Hfmj[m];
				elementMat[eq(m,reagents[0])][var(n,reagents[0])] -= Hfmj[n];
				elementMat[eq(m,reagents[0])][var(p,reagents[0])] -= Hfmj[p];
				for (unsigned j=0; j<nProducts; j++)
				{
					elementMat[eq(m,products[j])][var(m,reagents[0])] += Hfmj[m];
					elementMat[eq(m,products[j])][var(n,reagents[0])] += Hfmj[n];
					elementMat[eq(m,products[j])][var(p,reagents[0])] += Hfmj[p];
				}
			}
			else if (nReagents == 2) 
			{
				Hfmjk[m][m] = (R2[m][m][m][m]*kf[m][r] + R2[m][m][m][n]*kf[n][r] + R2[m][m][m][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[n][n] = (R2[m][n][n][m]*kf[m][r] + R2[m][n][n][n]*kf[n][r] + R2[m][n][n][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[p][p] = (R2[m][p][p][m]*kf[m][r] + R2[m][p][p][n]*kf[n][r] + R2[m][p][p][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[m][n] = (R2[m][m][n][m]*kf[m][r] + R2[m][m][n][n]*kf[n][r] + R2[m][m][n][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[n][m] = Hfmjk[m][n];
				Hfmjk[m][p] = (R2[m][m][p][m]*kf[m][r] + R2[m][m][p][n]*kf[n][r] + R2[m][m][p][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[p][m] = Hfmjk[m][p];
				Hfmjk[n][p] = (R2[m][n][p][m]*kf[m][r] + R2[m][n][p][n]*kf[n][r] + R2[m][n][p][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[p][n] = Hfmjk[n][p];

				for (unsigned j=0; j<nReagents; j++)
				{
					unsigned k = (j+1)%2;
					Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]] + Hfmjk[m][p]*concentrations[p][reagents[k]]);
					Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]] + Hfmjk[n][p]*concentrations[p][reagents[k]]);
					Hfmj[p] = 0.5*(Hfmjk[p][m]*concentrations[m][reagents[k]] + Hfmjk[p][n]*concentrations[n][reagents[k]] + Hfmjk[p][p]*concentrations[p][reagents[k]]);

					elementMat[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
					elementMat[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
					elementMat[eq(m,reagents[0])][var(p,reagents[j])] -= Hfmj[p];
					elementMat[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
					elementMat[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
					elementMat[eq(m,reagents[1])][var(p,reagents[j])] -= Hfmj[p];
					elementMat[eq(m,products[0])][var(m,reagents[j])] += Hfmj[m];
					elementMat[eq(m,products[0])][var(n,reagents[j])] += Hfmj[n];
					elementMat[eq(m,products[0])][var(p,reagents[j])] += Hfmj[p];
				}				
			}

			// Backward reaction
			if (nProducts == 1) 
			{
				Hbmj[m] = (R1[m][m][m]*kb[m][r] + R1[m][m][n]*kb[n][r] + R1[m][m][p]*kb[p][r])*elementSize180;
				Hbmj[n] = (R1[m][n][m]*kb[m][r] + R1[m][n][n]*kb[n][r] + R1[m][n][p]*kb[p][r])*elementSize180;
				Hbmj[p] = (R1[m][p][m]*kb[m][r] + R1[m][p][n]*kb[n][r] + R1[m][p][p]*kb[p][r])*elementSize180;

				elementMat[eq(m,products[0])][var(m,products[0])] -= Hbmj[m];
				elementMat[eq(m,products[0])][var(n,products[0])] -= Hbmj[n];
				elementMat[eq(m,products[0])][var(p,products[0])] -= Hbmj[p];
				for (unsigned j=0; j<nReagents; j++)
				{
					elementMat[eq(m,reagents[j])][var(m,products[0])] += Hbmj[m];
					elementMat[eq(m,reagents[j])][var(n,products[0])] += Hbmj[n];
					elementMat[eq(m,reagents[j])][var(p,products[0])] += Hbmj[p];
				}
			}
			else if (nProducts == 2) 
			{
				Hbmjk[m][m] = (R2[m][m][m][m]*kb[m][r] + R2[m][m][m][n]*kb[n][r] + R2[m][m][m][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[n][n] = (R2[m][n][n][m]*kb[m][r] + R2[m][n][n][n]*kb[n][r] + R2[m][n][n][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[p][p] = (R2[m][p][p][m]*kb[m][r] + R2[m][p][p][n]*kb[n][r] + R2[m][p][p][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[m][n] = (R2[m][m][n][m]*kb[m][r] + R2[m][m][n][n]*kb[n][r] + R2[m][m][n][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[n][m] = Hbmjk[m][n];
				Hbmjk[m][p] = (R2[m][m][p][m]*kb[m][r] + R2[m][m][p][n]*kb[n][r] + R2[m][m][p][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[p][m] = Hbmjk[m][p];
				Hbmjk[n][p] = (R2[m][n][p][m]*kb[m][r] + R2[m][n][p][n]*kb[n][r] + R2[m][n][p][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[p][n] = Hbmjk[n][p];

				for (unsigned j=0; j<nProducts; j++)
				{
					unsigned k = (j+1)%2;
					Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]] + Hbmjk[m][p]*concentrations[p][products[k]]);
					Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]] + Hbmjk[n][p]*concentrations[p][products[k]]);
					Hbmj[p] = 0.5*(Hbmjk[p][m]*concentrations[m][products[k]] + Hbmjk[p][n]*concentrations[n][products[k]] + Hbmjk[p][p]*concentrations[p][products[k]]);

					elementMat[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
					elementMat[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
					elementMat[eq(m,products[0])][var(p,products[j])] -= Hbmj[p];
					elementMat[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
					elementMat[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
					elementMat[eq(m,products[1])][var(p,products[j])] -= Hbmj[p];
					elementMat[eq(m,reagents[0])][var(m,products[j])] += Hbmj[m];
					elementMat[eq(m,reagents[0])][var(n,products[j])] += Hbmj[n];
					elementMat[eq(m,reagents[0])][var(p,products[j])] += Hbmj[p];
				}				
			}
		}
		delete[] products;
		delete[] reagents;
	}
}
//---------------------------------------------------------------------------
void HomReactionTerm_AX_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	
	// Calculate coefficients
	double R2[3][3][3][3];
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned r=0; r<nHomReactions; r++) 
		{
			kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
			kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
		}
		unsigned n = (m+1)%3;
		unsigned p = (m+2)%3;
		R2[m][m][m][m] = 60.*coordinates[m][0] + 12.*coordinates[n][0] + 12.*coordinates[p][0];
		R2[m][n][n][n] =  6.*coordinates[m][0] + 12.*coordinates[n][0] +  3.*coordinates[p][0];
		R2[m][p][p][p] =  6.*coordinates[m][0] +  3.*coordinates[n][0] + 12.*coordinates[p][0];
		R2[m][m][m][n] = 12.*coordinates[m][0] +  6.*coordinates[n][0] +  3.*coordinates[p][0];
		R2[m][m][n][m] = R2[m][m][m][n];
		R2[m][n][m][m] = R2[m][m][m][n];
		R2[m][m][m][p] = 12.*coordinates[m][0] +  3.*coordinates[n][0] +  6.*coordinates[p][0];
		R2[m][m][p][m] = R2[m][m][m][p];
		R2[m][p][m][m] = R2[m][m][m][p];
		R2[m][m][n][n] =  6.*coordinates[m][0] +  6.*coordinates[n][0] +  2.*coordinates[p][0];
		R2[m][n][m][n] = R2[m][m][n][n];
		R2[m][n][n][m] = R2[m][m][n][n];
		R2[m][m][p][p] =  6.*coordinates[m][0] +  2.*coordinates[n][0] +  6.*coordinates[p][0];
		R2[m][p][m][p] = R2[m][m][p][p];
		R2[m][p][p][m] = R2[m][m][p][p];
		R2[m][m][n][p] =  3.*coordinates[m][0] +  2.*coordinates[n][0] +  2.*coordinates[p][0];
		R2[m][m][p][n] = R2[m][m][n][p];
		R2[m][n][m][p] = R2[m][m][n][p];
		R2[m][p][m][n] = R2[m][m][n][p];
		R2[m][n][p][m] = R2[m][m][n][p];
		R2[m][p][n][m] = R2[m][m][n][p];
		R2[m][n][n][p] =  2.*coordinates[m][0] +  3.*coordinates[n][0] +  2.*coordinates[p][0];
		R2[m][n][p][n] = R2[m][n][n][p];
		R2[m][p][n][n] = R2[m][n][n][p];
		R2[m][p][p][n] =  2.*coordinates[m][0] +  2.*coordinates[n][0] +  3.*coordinates[p][0];
		R2[m][p][n][p] = R2[m][p][p][n];
		R2[m][n][p][p] = R2[m][p][p][n];
	}
	elementSize = elementProps->calcSize(coordinates);
	//double elementSize2520 = elementSize/2520.; BUG?
	double elementSize1260 = elementSize/1260;
	
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
				Hfmjk[m][m] = (R2[m][m][m][m]*kf[m][r] + R2[m][m][m][n]*kf[n][r] + R2[m][m][m][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[n][n] = (R2[m][n][n][m]*kf[m][r] + R2[m][n][n][n]*kf[n][r] + R2[m][n][n][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[p][p] = (R2[m][p][p][m]*kf[m][r] + R2[m][p][p][n]*kf[n][r] + R2[m][p][p][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[m][n] = (R2[m][m][n][m]*kf[m][r] + R2[m][m][n][n]*kf[n][r] + R2[m][m][n][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[n][m] = Hfmjk[m][n];
				Hfmjk[m][p] = (R2[m][m][p][m]*kf[m][r] + R2[m][m][p][n]*kf[n][r] + R2[m][m][p][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[p][m] = Hfmjk[m][p];
				Hfmjk[n][p] = (R2[m][n][p][m]*kf[m][r] + R2[m][n][p][n]*kf[n][r] + R2[m][n][p][p]*kf[p][r])*elementSize1260;//elementSize2520; BUG?
				Hfmjk[p][n] = Hfmjk[n][p];

				for (unsigned j=0; j<nReagents; j++)
				{
					unsigned k = (j+1)%2;
					Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]] + Hfmjk[m][p]*concentrations[p][reagents[k]]);
					Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]] + Hfmjk[n][p]*concentrations[p][reagents[k]]);
					Hfmj[p] = 0.5*(Hfmjk[p][m]*concentrations[m][reagents[k]] + Hfmjk[p][n]*concentrations[n][reagents[k]] + Hfmjk[p][p]*concentrations[p][reagents[k]]);

					elementJac[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
					elementJac[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
					elementJac[eq(m,reagents[0])][var(p,reagents[j])] -= Hfmj[p];
					elementJac[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
					elementJac[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
					elementJac[eq(m,reagents[1])][var(p,reagents[j])] -= Hfmj[p];
					elementJac[eq(m,products[0])][var(m,reagents[j])] += Hfmj[m];
					elementJac[eq(m,products[0])][var(n,reagents[j])] += Hfmj[n];
					elementJac[eq(m,products[0])][var(p,reagents[j])] += Hfmj[p];
				}				
			}

			// Backward reaction
			if (nProducts == 2) 
			{
				Hbmjk[m][m] = (R2[m][m][m][m]*kb[m][r] + R2[m][m][m][n]*kb[n][r] + R2[m][m][m][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[n][n] = (R2[m][n][n][m]*kb[m][r] + R2[m][n][n][n]*kb[n][r] + R2[m][n][n][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[p][p] = (R2[m][p][p][m]*kb[m][r] + R2[m][p][p][n]*kb[n][r] + R2[m][p][p][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[m][n] = (R2[m][m][n][m]*kb[m][r] + R2[m][m][n][n]*kb[n][r] + R2[m][m][n][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[n][m] = Hbmjk[m][n];
				Hbmjk[m][p] = (R2[m][m][p][m]*kb[m][r] + R2[m][m][p][n]*kb[n][r] + R2[m][m][p][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[p][m] = Hbmjk[m][p];
				Hbmjk[n][p] = (R2[m][n][p][m]*kb[m][r] + R2[m][n][p][n]*kb[n][r] + R2[m][n][p][p]*kb[p][r])*elementSize1260;//elementSize2520; BUG?
				Hbmjk[p][n] = Hbmjk[n][p];

				for (unsigned j=0; j<nProducts; j++)
				{
					unsigned k = (j+1)%2;
					Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]] + Hbmjk[m][p]*concentrations[p][products[k]]);
					Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]] + Hbmjk[n][p]*concentrations[p][products[k]]);
					Hbmj[p] = 0.5*(Hbmjk[p][m]*concentrations[m][products[k]] + Hbmjk[p][n]*concentrations[n][products[k]] + Hbmjk[p][p]*concentrations[p][products[k]]);

					elementJac[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
					elementJac[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
					elementJac[eq(m,products[0])][var(p,products[j])] -= Hbmj[p];
					elementJac[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
					elementJac[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
					elementJac[eq(m,products[1])][var(p,products[j])] -= Hbmj[p];
					elementJac[eq(m,reagents[0])][var(m,products[j])] += Hbmj[m];
					elementJac[eq(m,reagents[0])][var(n,products[j])] += Hbmj[n];
					elementJac[eq(m,reagents[0])][var(p,products[j])] += Hbmj[p];
				}				
			}
		}
		delete[] products;
		delete[] reagents;
	}
}
//---------------------------------------------------------------------------

