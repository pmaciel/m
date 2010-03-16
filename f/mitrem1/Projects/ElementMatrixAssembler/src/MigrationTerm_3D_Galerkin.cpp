//---------------------------------------------------------------------------

#include "MigrationTerm_3D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_3D_Galerkin::MigrationTerm_3D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_3D_Galerkin::~MigrationTerm_3D_Galerkin()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_3D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction *= 0.25;
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
	elementSize = elementProps->calcSize(coordinates);
	double elementSize180 = 180.*elementSize;

	// Add to element matrix
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			unsigned p = (m+1)%4;
			unsigned q = (m+2)%4;
			unsigned r = (m+3)%4;
			for (unsigned c=0; c<nDimensions; c++)
			{
				Mmi[c] = ((2.*W[m][i] + W[p][i] + W[q][i] + W[r][i])*concentrations[m][i]
						 + (W[m][i] + 2.*W[p][i] + W[q][i] + W[r][i])*concentrations[p][i]
						 + (W[m][i] + W[p][i] + 2.*W[q][i] + W[r][i])*concentrations[q][i]
						 + (W[m][i] + W[p][i] + W[q][i] + 2.*W[r][i])*concentrations[r][i])
						 *normals[m][c]/elementSize180;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] -= Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1] + Mmi[2]*normals[n][2];
			}
		}
	}
}
//---------------------------------------------------------------------------
void MigrationTerm_3D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction *= 0.25;
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
	double elementSize3 = 3.*elementSize;
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] /= elementSize3;
	}

	// Add to element jacobian
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			unsigned p = (m+1)%4;
			unsigned q = (m+2)%4;
			unsigned r = (m+3)%4;
			for (unsigned c=0; c<nDimensions; c++)
			{
				Wmi[m][c] = (2.*W[m][i] + W[p][i] + W[q][i] + W[r][i])*normals[m][c]/60.;
				Wmi[p][c] = (W[m][i] + 2.*W[p][i] + W[q][i] + W[r][i])*normals[m][c]/60.;
				Wmi[q][c] = (W[m][i] + W[p][i] + 2.*W[q][i] + W[r][i])*normals[m][c]/60.;
				Wmi[r][c] = (W[m][i] + W[p][i] + W[q][i] + 2.*W[r][i])*normals[m][c]/60.;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0] + gradU[1]*Wmi[n][1] + gradU[2]*Wmi[n][2];
			}
		}
	}
}
//---------------------------------------------------------------------------

