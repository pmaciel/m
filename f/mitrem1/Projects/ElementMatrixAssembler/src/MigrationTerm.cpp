//---------------------------------------------------------------------------

#include "MigrationTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm::MigrationTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ElementTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	W = new double*[nNodes];
	Wmi = new EmptyDoubleVector[nNodes];
	for (unsigned n=0; n<nNodes; n++)
	{
		W[n] = new double[nIons];
		Wmi[n] = new double[nDimensions];
	}
	Mmi = new double[nDimensions];
	gradU = new double[nDimensions];
}
//---------------------------------------------------------------------------
MigrationTerm::~MigrationTerm()
{
	delete[] gradU;
	delete[] Mmi;
	for (unsigned n=0; n<nNodes; n++)
	{
		delete[] Wmi[n];
		delete[] W[n];
	}
	delete[] Wmi;
	delete[] W;
}
//---------------------------------------------------------------------------
void MigrationTerm::calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= nNodes;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);
	
	// Calculate coefficients
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] = 0;
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		normal = elementProps->calcNormal(m,coordinates);
		for (unsigned c=0; c<nDimensions; c++)
		{
			gradU[c] += normal[c]*potentials[m];
		}
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned i=0; i<nIons; i++) 
		{
			W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] /= nDimensions*elementSize;
	}
	
	// Add to element current density
	for (unsigned i=0; i<nIons; i++) 
	{
		double ziF = F_CONST*mitrem->getIonChargeNumber(i);
		double Wi = 0.;
		double ci = 0.;
		for (unsigned m=0; m<nNodes; m++) 
		{
			Wi += W[m][i];
			ci += concentrations[m][i];
		}
		Wi /= nNodes;
		ci /= nNodes;
		for (unsigned c=0; c<nDimensions; c++) 
		{			
			elementCurr[i][c] -= ziF*Wi*ci*gradU[c];
		}
	}	
}
//---------------------------------------------------------------------------

