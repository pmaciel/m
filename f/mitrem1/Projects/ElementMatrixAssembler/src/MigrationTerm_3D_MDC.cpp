//---------------------------------------------------------------------------

#include "MigrationTerm_3D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_3D_MDC::MigrationTerm_3D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_3D_MDC::~MigrationTerm_3D_MDC()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_3D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
	double elementSize15552 = 15552.*elementSize;

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
				Mmi[c] = (230*normals[m][c]*W[m][i]*concentrations[m][i]
						 + (119.*normals[m][c] - 37.*normals[p][c])*W[p][i]*concentrations[p][i]
						 + (119.*normals[m][c] - 37.*normals[q][c])*W[q][i]*concentrations[q][i]
						 + (119.*normals[m][c] - 37.*normals[r][c])*W[r][i]*concentrations[r][i]
						 + (92.*normals[m][c] - 46.*normals[p][c])*(W[m][i]*concentrations[p][i] + W[p][i]*concentrations[m][i])
						 + (92.*normals[m][c] - 46.*normals[q][c])*(W[m][i]*concentrations[q][i] + W[q][i]*concentrations[m][i])
						 + (92.*normals[m][c] - 46.*normals[r][c])*(W[m][i]*concentrations[r][i] + W[r][i]*concentrations[m][i])
						 + (69.*normals[m][c] + 13.*normals[r][c])*(W[p][i]*concentrations[q][i] + W[q][i]*concentrations[p][i])
						 + (69.*normals[m][c] + 13.*normals[q][c])*(W[r][i]*concentrations[p][i] + W[p][i]*concentrations[r][i])
						 + (69.*normals[m][c] + 13.*normals[p][c])*(W[q][i]*concentrations[r][i] + W[r][i]*concentrations[q][i]))/elementSize15552;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] -= (Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1] + Mmi[2]*normals[n][2]);
			}
		}
	}
}
//---------------------------------------------------------------------------
void MigrationTerm_3D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
				Wmi[m][c] = (230.*normals[m][c]*W[m][i]
							+ (119.*normals[m][c] - 37.*normals[p][c])*W[p][i]
							+ (119.*normals[m][c] - 37.*normals[q][c])*W[q][i]
							+ (119.*normals[m][c] - 37.*normals[r][c])*W[r][i])/5184.;		
				Wmi[p][c] = ((119.*normals[m][c] - 37.*normals[p][c])*W[m][i]
							+ (92.*normals[m][c] - 46.*normals[p][c])*W[p][i]
							+ (69.*normals[m][c] + 13.*normals[r][c])*W[q][i]
							+ (69.*normals[m][c] + 13.*normals[q][c])*W[r][i])/5184.;
				Wmi[q][c] = ((119.*normals[m][c] - 37.*normals[q][c])*W[m][i]
							+ (69.*normals[m][c] + 13.*normals[r][c])*W[p][i]
							+ (92.*normals[m][c] - 46.*normals[q][c])*W[q][i]
							+ (69.*normals[m][c] + 13.*normals[q][c])*W[r][i])/5184.;	
				Wmi[r][c] = ((119.*normals[m][c] - 37.*normals[r][c])*W[m][i]					  
							+ (69.*normals[m][c] + 13.*normals[q][c])*W[p][i]
							+ (69.*normals[m][c] + 13.*normals[p][c])*W[q][i]
							+ (92.*normals[m][c] - 46.*normals[r][c])*W[r][i])/5184.;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0] + gradU[1]*Wmi[n][1] + gradU[2]*Wmi[n][2];
			}
		}
	}
}
//---------------------------------------------------------------------------

