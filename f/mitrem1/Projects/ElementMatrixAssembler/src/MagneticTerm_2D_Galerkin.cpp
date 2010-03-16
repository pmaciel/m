//---------------------------------------------------------------------------

#include "MagneticTerm_2D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MagneticTerm_2D_Galerkin::MagneticTerm_2D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MagneticTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MagneticTerm_2D_Galerkin::~MagneticTerm_2D_Galerkin()
{
}
//---------------------------------------------------------------------------
void MagneticTerm_2D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);

	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		normal = elementProps->calcNormal(m,coordinates);
		for (unsigned c=0; c<nDimensions; c++)
		{
			normals[m][c] = normal[c];
		}		
		for (unsigned i=0; i<nIons; i++) 
		{
			W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		}
	}

	// Add to element matrix
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			for (unsigned n=0; n<nNodes; n++) 
			{
				unsigned p = (n+1)%3;
				unsigned q = (n+2)%3;
				for (unsigned c=0; c<nDimensions; c++)
				{
					unsigned d = (c+1)%2;
					double sign = 1.-2.*(c%2); //if(c==0) sign=+1; else if(c==1) sign=-1;
					cCoefficients[n][c] = sign*(magneticFieldVectors[n][2]*(
																	3.*velocities[n][d]*(4.*W[n][i] + W[p][i] + W[q][i])
																	+ velocities[p][d]*(3.*W[n][i] + 2.*W[p][i] + W[q][i])
																	+ velocities[q][d]*(3.*W[n][i] + W[p][i] + 2.*W[q][i]))
																+ magneticFieldVectors[p][2]*(
																	velocities[n][d]*(3.*W[n][i] + 2.*W[p][i] + W[q][i])
																	+ velocities[p][d]*(2.*W[n][i] + 3.*W[p][i] + W[q][i])
																	+ velocities[q][d]*(W[n][i] + W[p][i] + W[q][i]))
																+ magneticFieldVectors[q][2]*(
																	velocities[n][d]*(3.*W[n][i] + W[p][i] + 2.*W[q][i])
																	+ velocities[p][d]*(W[n][i] + W[p][i] + W[q][i])
																	+ velocities[q][d]*(2.*W[n][i] + W[p][i] + 3.*W[q][i]))) / 360.;
				}
				elementMat[eqmi][var(n,i)] += cCoefficients[n][0]*normals[m][0] + cCoefficients[n][1]*normals[m][1];
			}
		}
	}
}
//---------------------------------------------------------------------------
void MagneticTerm_2D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

