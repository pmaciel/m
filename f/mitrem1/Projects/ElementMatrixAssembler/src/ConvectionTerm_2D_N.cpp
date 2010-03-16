//---------------------------------------------------------------------------

#include "ConvectionTerm_2D_N.h"

#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ConvectionTerm_2D_N::ConvectionTerm_2D_N(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ConvectionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ConvectionTerm_2D_N::~ConvectionTerm_2D_N()
{
}
//---------------------------------------------------------------------------
void ConvectionTerm_2D_N::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	
	// Calculate coefficients
	unsigned nTargetNodes = 0;
	unsigned targetNodes[2];
	unsigned nOtherNodes = 0;
	unsigned otherNodes[3];
	for (unsigned d=0; d<nDimensions; d++)
	{
		averageVelocity[d] = 0;
		for (unsigned m=0; m<nNodes; m++)
		{
			averageVelocity[d] += velocities[m][d];
		}
		averageVelocity[d] /= 3.;
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		normal = elementProps->calcNormal(m,coordinates);
		k[m] = 0.5*(normal[0]*averageVelocity[0] + normal[1]*averageVelocity[1]);
		if (k[m] > 0.) 
		{
			targetNodes[nTargetNodes] = m;
			nTargetNodes++;
		}
		else 
		{
			otherNodes[nOtherNodes] = m;
			nOtherNodes++;
		}
	}

	// Add to element matrix
	if (nTargetNodes == 1) 
	{
		for (unsigned i=0; i<nIons; i++) 
		{
			for (unsigned m=0; m<nNodes; m++) 
			{
				elementMat[eq(targetNodes[0],i)][var(m,i)] -= k[m]*(1.-volumeGasFraction);
			}
		}
	}
	else if (nTargetNodes == 2) 
	{
		for (unsigned i=0; i<nIons; i++) 
		{
			elementMat[eq(targetNodes[0],i)][var(targetNodes[0],i)] -= k[targetNodes[0]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[0],i)][var(otherNodes[0],i)]  += k[targetNodes[0]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(targetNodes[1],i)] -= k[targetNodes[1]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(otherNodes[0],i)]  += k[targetNodes[1]]*(1.-volumeGasFraction);
		}
	}
}
//---------------------------------------------------------------------------
void ConvectionTerm_2D_N::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

