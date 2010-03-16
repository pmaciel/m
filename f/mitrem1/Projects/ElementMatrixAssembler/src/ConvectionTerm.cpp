//---------------------------------------------------------------------------

#include "ConvectionTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ConvectionTerm::ConvectionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ElementTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	Vm = new double[nDimensions];
	k = new double[nNodes];
	averageVelocity = new double[nDimensions];
}
//---------------------------------------------------------------------------
ConvectionTerm::~ConvectionTerm()
{
	delete[] Vm;
	delete[] k;
	delete[] averageVelocity;
}
//---------------------------------------------------------------------------
void ConvectionTerm::calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= nNodes;
	
	// Calculate coefficients
	for (unsigned d=0; d<nDimensions; d++)
	{
		averageVelocity[d] = 0;
		for (unsigned m=0; m<nNodes; m++)
		{
			averageVelocity[d] += velocities[m][d];
		}
		averageVelocity[d] /= nNodes;
	}

	// Add to element current density
	for (unsigned i=0; i<nIons; i++) 
	{
		double ziF = F_CONST*mitrem->getIonChargeNumber(i);
		double ci = 0.;
		for (unsigned m=0; m<nNodes; m++) 
		{
			ci += concentrations[m][i];
		}
		ci /= nNodes;
		for (unsigned c=0; c<nDimensions; c++) 
		{			
			elementCurr[i][c] += ziF*averageVelocity[c]*ci*(1.-volumeGasFraction);
		}
	}
}
//---------------------------------------------------------------------------

