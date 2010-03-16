//---------------------------------------------------------------------------

#include "MagneticTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MagneticTerm::MagneticTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ElementTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	W = new double*[nNodes];
	cCoefficients = new double*[nNodes];
	for (unsigned n=0; n<nNodes; n++)
	{
		W[n] = new double[nIons];
		cCoefficients[n] = new double[nDimensions];
	}
	vCrossB = new double[nDimensions];
}
//---------------------------------------------------------------------------
MagneticTerm::~MagneticTerm()
{
	delete[] vCrossB;
	for (unsigned n=0; n<nNodes; n++)
	{
		delete[] W[n];
		delete[] cCoefficients[n];
	}
	delete[] W;
	delete[] cCoefficients;
}
//---------------------------------------------------------------------------
void MagneticTerm::calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= nNodes;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);

	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned i=0; i<nIons; i++) 
		{
			W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		}
	}

	// PRODUCT OF AVERAGES or AVERAGE OF PRODUCTS? We chose product of averages...
	if (nDimensions == 3)
	{
		// Calculate average velocity and average magnetic field vector
		double* averageVelocity = new double[3];
		double* averageMagneticFieldVector = new double[3];
		for (unsigned c=0; c<3; c++)
		{
			averageVelocity[c] = 0;
			averageMagneticFieldVector[c] = 0;
		}
		for (unsigned c=0; c<3; c++)
		{
			for (unsigned m=0; m<nNodes; m++) 
			{
				averageVelocity[c] += velocities[m][c];
				averageMagneticFieldVector[c] += magneticFieldVectors[m][c];
			}
			averageVelocity[c] /= nNodes;
			averageMagneticFieldVector[c] /= nNodes;
		}
		// Calculate cross product of v and B
		for (unsigned c=0; c<3; c++)
		{
			unsigned d = (c+1)%3;
			unsigned e = (c+2)%3;
			vCrossB[c] = averageVelocity[d]*averageMagneticFieldVector[e] - averageVelocity[e]*averageMagneticFieldVector[d];
		}
		delete[] averageVelocity;
		delete[] averageMagneticFieldVector;
	}

	else if (nDimensions == 2)
	{
		// Calculate average velocity and average magnetic field vector
		double* averageVelocity = new double[2];
		double averageMagneticFieldVector = 0;
		averageVelocity[0] = 0;
		averageVelocity[1] = 0;
		for (unsigned c=0; c<2; c++)
		{
			for (unsigned m=0; m<nNodes; m++) 
			{
				averageVelocity[c] += velocities[m][c];
			}
			averageVelocity[c] /= nNodes;
		}
		for (unsigned m=0; m<nNodes; m++) 
		{
			averageMagneticFieldVector += magneticFieldVectors[m][2];
		}
		averageMagneticFieldVector /= nNodes;
		// Calculate cross product of v and B
		vCrossB[0] = averageVelocity[1]*averageMagneticFieldVector;
		vCrossB[1] = - averageVelocity[0]*averageMagneticFieldVector;
		delete[] averageVelocity;
	}

	else if (nDimensions == 1)
	{
		vCrossB[0] = 0;
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
			elementCurr[i][c] += ziF*Wi*ci*vCrossB[c];
		}
	}	
}
//---------------------------------------------------------------------------

