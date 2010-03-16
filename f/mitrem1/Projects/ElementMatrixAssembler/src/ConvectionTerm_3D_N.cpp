//---------------------------------------------------------------------------

#include "ConvectionTerm_3D_N.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ConvectionTerm_3D_N::ConvectionTerm_3D_N(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ConvectionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ConvectionTerm_3D_N::~ConvectionTerm_3D_N()
{
}
//---------------------------------------------------------------------------
void ConvectionTerm_3D_N::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction *= 0.25;
	
	// Calculate coefficients
	unsigned nTargetNodes = 0;
	unsigned targetNodes[3];
	unsigned nOtherNodes = 0;
	unsigned otherNodes[4];
	for (unsigned d=0; d<nDimensions; d++)
	{
		averageVelocity[d] = 0;
		for (unsigned m=0; m<nNodes; m++)
		{
			averageVelocity[d] += velocities[m][d];
		}
		averageVelocity[d] *= 0.25;
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		normal = elementProps->calcNormal(m,coordinates);
		k[m] = (normal[0]*averageVelocity[0] + normal[1]*averageVelocity[1] + normal[2]*averageVelocity[2])/3.;
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
		double factorTarget0 = k[targetNodes[0]]/(k[targetNodes[0]] + k[targetNodes[1]]);
		double factorTarget1 = k[targetNodes[1]]/(k[targetNodes[0]] + k[targetNodes[1]]);
		for (unsigned i=0; i<nIons; i++) 
		{
			elementMat[eq(targetNodes[0],i)][var(targetNodes[0],i)] -= k[targetNodes[0]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[0],i)][var(otherNodes[0],i)]  -= k[otherNodes[0]]*factorTarget0*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[0],i)][var(otherNodes[1],i)]  -= k[otherNodes[1]]*factorTarget0*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(targetNodes[1],i)] -= k[targetNodes[1]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(otherNodes[0],i)]  -= k[otherNodes[0]]*factorTarget1*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(otherNodes[1],i)]  -= k[otherNodes[1]]*factorTarget1*(1.-volumeGasFraction);
		}
	}
	else if (nTargetNodes == 3) 
	{
		for (unsigned i=0; i<nIons; i++) 
		{
			elementMat[eq(targetNodes[0],i)][var(targetNodes[0],i)] -= k[targetNodes[0]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[0],i)][var(otherNodes[0],i)]  += k[targetNodes[0]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(targetNodes[1],i)] -= k[targetNodes[1]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[1],i)][var(otherNodes[0],i)]  += k[targetNodes[1]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[2],i)][var(targetNodes[2],i)] -= k[targetNodes[2]]*(1.-volumeGasFraction);
			elementMat[eq(targetNodes[2],i)][var(otherNodes[0],i)]  += k[targetNodes[2]]*(1.-volumeGasFraction);
		}
	}
}
//---------------------------------------------------------------------------
void ConvectionTerm_3D_N::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

