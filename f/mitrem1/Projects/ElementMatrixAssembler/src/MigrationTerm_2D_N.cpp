//---------------------------------------------------------------------------

#include "MigrationTerm_2D_N.h"
#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_2D_N::MigrationTerm_2D_N(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	gradc = new double*[nIons];
	for (unsigned i=0; i<nIons; i++)
	{
		gradc[i] = new double[nDimensions];
	}
}
//---------------------------------------------------------------------------
MigrationTerm_2D_N::~MigrationTerm_2D_N()
{
	for (unsigned i=0; i<nIons; i++)
	{
		delete[] gradc[i];
	}
	delete[] gradc;
}
//---------------------------------------------------------------------------
void MigrationTerm_2D_N::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);
	
	// Calculate coefficients
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned c=0; c<nDimensions; c++)
		{
			gradc[i][c] = 0.;
		}
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		normal = elementProps->calcNormal(m,coordinates);
		for (unsigned i=0; i<nIons; i++) 
		{
			for (unsigned c=0; c<nDimensions; c++)
			{
				gradc[i][c] += normal[c]*concentrations[m][i];
			}
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned c=0; c<nDimensions; c++)
		{
			gradc[i][c] /= nDimensions*elementSize;
		}
	}


	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] = 0;
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		normal = elementProps->calcNormal(m,coordinates);
		for (unsigned c=0; c<nDimensions; c++)
		{
			normals[m][c] = normal[c];
			gradU[c] += normal[c]*potentials[m];
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize2 = 2.*elementSize;
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] /= elementSize2;
	}


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
	elementSize = elementProps->calcSize(coordinates);
	double elementSize48 = 48.*elementSize;

//	std::cout << "Matrix" << std::endl;
//	std::cout << "extra" << std::endl;

	// Add to element matrix
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			unsigned p = (m+1)%3;
			unsigned q = (m+2)%3;	
			for (unsigned c=0; c<nDimensions; c++)
			{
				Wmi[m][c] = (2.*W[m][i] + W[p][i] + W[q][i]);//*normals[m][c]/24.;

				Mmi[c] = ((2.*W[m][i] + W[p][i] + W[q][i])*concentrations[m][i]
						 + (W[m][i] + 2.*W[p][i] + W[q][i])*concentrations[p][i]
						 + (W[m][i] + W[p][i] + 2.*W[q][i])*concentrations[q][i])
						 *normals[m][c]/elementSize48;

			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] -= Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1];
//				std::cout << "Mat[" << eqmi << "][" << var(n,nIons) << "]: " << -(Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1]) << std::endl;
				elementMat[eqmi][var(n,nIons)] -= gradc[i][0]*Wmi[m][0]*normals[n][0]/24. + gradc[i][1]*Wmi[m][1]*normals[n][1]/24.;
//				std::cout << "Mat[" << eqmi << "][" << var(n,nIons) << "]: " << -(gradc[i][0]*Wmi[m][0]*normals[n][0]/24. + gradc[i][1]*Wmi[m][1]*normals[n][1]/24.) << std::endl;
			}
		}
	}
}
//---------------------------------------------------------------------------
void MigrationTerm_2D_N::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
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
			normals[m][c] = normal[c];
			gradU[c] += normal[c]*potentials[m];
		}
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned i=0; i<nIons; i++) 
		{
			W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize2 = 2.*elementSize;
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] /= elementSize2;
	}

//	std::cout << "Jacobian" << std::endl;
//	std::cout << "klassiek" << std::endl;
	// Add to element jacobian
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			unsigned p = (m+1)%3;
			unsigned q = (m+2)%3;
			for (unsigned c=0; c<nDimensions; c++)
			{
				Wmi[m][c] = (2.*W[m][i] + W[p][i] + W[q][i]);//*normals[m][c]/24.;
				Wmi[p][c] = (W[m][i] + 2.*W[p][i] + W[q][i]);//*normals[m][c]/24.;
				Wmi[q][c] = (W[m][i] + W[p][i] + 2.*W[q][i]);//*normals[m][c]/24.;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0]*normals[m][0]/24. + gradU[1]*Wmi[n][1]*normals[m][1]/24.;
//				std::cout << "Jac[" << eqmi << "][" << var(n,i) << "]: " << -(gradU[0]*Wmi[n][0]*normals[m][0]/24. + gradU[1]*Wmi[n][1]*normals[m][1]/24.) << std::endl;
				elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[m][0]*normals[n][0]/24. + gradU[1]*Wmi[m][1]*normals[n][1]/24.;
//				std::cout << "Jac[" << eqmi << "][" << var(n,i) << "]: " << -(gradU[0]*Wmi[m][0]*normals[n][0]/24. + gradU[1]*Wmi[m][1]*normals[n][1]/24.) << std::endl;
			}
		}
	}
}
//---------------------------------------------------------------------------

