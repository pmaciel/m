//---------------------------------------------------------------------------

#include "ConvectionTerm_1D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ConvectionTerm_1D_Galerkin::ConvectionTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ConvectionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	Nu = new double*[nNodes];
	for (unsigned m=0; m<nNodes; m++)
	{
		Nu[m] = new double[nDimensions];
	}
}
//---------------------------------------------------------------------------
ConvectionTerm_1D_Galerkin::~ConvectionTerm_1D_Galerkin()
{
	for (unsigned m=0; m<nNodes; m++)
	{
		delete[] Nu[m];
	}
	delete[] Nu;
}
//---------------------------------------------------------------------------
void ConvectionTerm_1D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions)
{
	double voidFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		voidFraction += voidFractions[m];
	}
	voidFraction *= 0.5;
	
	// Calculate coefficients
	unsigned nTargetNodes = 0;
	unsigned targetNode = 0;
	//unsigned nOtherNodes = 0;
	//unsigned otherNodes[2];
	/*for (unsigned d=0; d<nDimensions; d++)
	{
		averageVelocity[d] = 0;
		for (unsigned m=0; m<nNodes; m++)
		{
			averageVelocity[d] += velocities[m][d];
		}
		averageVelocity[d] *= 0.5;
	}*/

	for (unsigned m=0; m<nNodes; m++)
	{
		for (unsigned d=0; d<nDimensions; d++)
		{
			Nu[m][d] = 2*velocities[m][d] + velocities[(m+1)%nNodes][d];
		}
	}

	// Add to element matrix
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			for (unsigned j=0; j<nIons; j++) 
			{
				//Dmij[0] = (D[0][i][j] + D[1][i][j])*normals[m][0]/elementSize2;
				for (unsigned n=0; n<nNodes; n++) 
				{
					elementMat[eqmi][var(n,j)] += 1.0/6*Nu[j][0]*normals[m][0];
				}
			}
		}
	}

	
	for (unsigned m=0; m<nNodes; m++) 
	{
		normal = elementProps->calcNormal(m,coordinates);
		k[m] = normal[0]*averageVelocity[0];
		if (k[m] > 0.) 
		{
			targetNode = m;
			nTargetNodes++;
		}
		//else 
		//{
		//	otherNodes[nOtherNodes] = m;
		//	nOtherNodes++;
		//}
	}

	// Add to element matrix
	if (nTargetNodes == 1) 
	{
		for (unsigned i=0; i<nIons; i++) 
		{
			for (unsigned m=0; m<nNodes; m++) 
			{
				elementMat[eq(targetNode,i)][var(m,i)] -= k[m]*(1.-voidFraction);
			}
		}
	}
}
//---------------------------------------------------------------------------
void ConvectionTerm_1D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions)
{
}
//---------------------------------------------------------------------------

