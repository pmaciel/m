//---------------------------------------------------------------------------

#include "MigrationTerm_AX_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_AX_Galerkin::MigrationTerm_AX_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_AX_Galerkin::~MigrationTerm_AX_Galerkin()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_AX_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);
	
	// Calculate coefficients
	double R[3][3];
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
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		R[m][m] = 6.*coordinates[m][0] + 2.*coordinates[p][0] + 2.*coordinates[q][0];
		R[m][p] = 2.*coordinates[m][0] + 2.*coordinates[p][0] +    coordinates[q][0];
		R[m][q] = 2.*coordinates[m][0] +    coordinates[p][0] + 2.*coordinates[q][0];
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize240 = 240.*elementSize;

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
				Mmi[c] = ((R[m][m]*W[m][i] + R[m][p]*W[p][i] + R[m][q]*W[q][i])*concentrations[m][i]
						 + (R[p][m]*W[m][i] + R[p][p]*W[p][i] + R[p][q]*W[q][i])*concentrations[p][i]
						 + (R[q][m]*W[m][i] + R[q][p]*W[p][i] + R[q][q]*W[q][i])*concentrations[q][i])
						 *normals[m][c]/elementSize240;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] -= Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1];
			}
		}
	}
}
//---------------------------------------------------------------------------
void MigrationTerm_AX_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);

	// Calculate coefficients
	double R[3][3];
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
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		R[m][m] = 6.*coordinates[m][0] + 2.*coordinates[p][0] + 2.*coordinates[q][0];
		R[m][p] = 2.*coordinates[m][0] + 2.*coordinates[p][0] +    coordinates[q][0];
		R[m][q] = 2.*coordinates[m][0] +    coordinates[p][0] + 2.*coordinates[q][0];
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize2 = 2.*elementSize;
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] /= elementSize2;
	}

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
				Wmi[m][c] = (R[m][m]*W[m][i] + R[m][p]*W[p][i] + R[m][q]*W[q][i])*normals[m][c]/120.;
				Wmi[p][c] = (R[p][m]*W[m][i] + R[p][p]*W[p][i] + R[p][q]*W[q][i])*normals[m][c]/120.;
				Wmi[q][c] = (R[q][m]*W[m][i] + R[q][p]*W[p][i] + R[q][q]*W[q][i])*normals[m][c]/120.;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0] + gradU[1]*Wmi[n][1];
			}
		}
	}
}
//---------------------------------------------------------------------------

