//---------------------------------------------------------------------------

#include "MigrationTerm_AX_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_AX_MDC::MigrationTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_AX_MDC::~MigrationTerm_AX_MDC()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_AX_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);
	
	// Calculate coefficients
	double R[3][3][3][2];
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
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		for (unsigned c=0; c<nDimensions; c++)
		{
			R[m][m][m][c] = 195.*normals[m][c]*coordinates[m][1]
							+ (109.*normals[m][c] - 43.*normals[p][c])*coordinates[p][1]
							+ (109.*normals[m][c] - 43.*normals[q][c])*coordinates[q][1];
			R[m][p][p][c] = (89.*normals[m][c] - 53.*normals[p][c])*coordinates[m][1]
							+ (81.*normals[m][c] - 57.*normals[p][c])*coordinates[p][1]
							+ (46.*normals[m][c] - 10.*normals[p][c])*coordinates[q][1];
			R[m][q][q][c] = (89.*normals[m][c] - 53.*normals[q][c])*coordinates[m][1]
							+ (46.*normals[m][c] - 10.*normals[q][c])*coordinates[p][1]
							+ (81.*normals[m][c] - 57.*normals[q][c])*coordinates[q][1];
			R[m][m][p][c] = (109.*normals[m][c] - 43.*normals[p][c])*coordinates[m][1]
							+ (89.*normals[m][c] - 53.*normals[p][c])*coordinates[p][1]
							+ 66.*normals[m][c]*coordinates[q][1];
			R[m][p][m][c] = R[m][m][p][c];
			R[m][m][q][c] = (109.*normals[m][c] - 43.*normals[q][c])*coordinates[m][1]
							+ 66.*normals[m][c]*coordinates[p][1]
							+ (89.*normals[m][c] - 53.*normals[q][c])*coordinates[q][1];
			R[m][q][m][c] = R[m][m][p][c];
			R[m][p][q][c] = 66.*normals[m][c]*coordinates[m][1]
							+ (46.*normals[m][c] - 10.*normals[p][c])*coordinates[p][1]
							+ (46.*normals[m][c] - 10.*normals[q][c])*coordinates[q][1];
			R[m][q][p][c] = R[m][p][q][c];
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize10368 = 10368.*elementSize;

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
				Mmi[c] = (R[m][m][m][c]*W[m][i]*concentrations[m][i]
						 + R[m][p][p][c]*W[p][i]*concentrations[p][i]
						 + R[m][q][q][c]*W[q][i]*concentrations[q][i]
						 + R[m][m][p][c]*(W[m][i]*concentrations[p][i] + W[p][i]*concentrations[m][i])
						 + R[m][m][q][c]*(W[m][i]*concentrations[q][i] + W[q][i]*concentrations[m][i])
						 + R[m][p][q][c]*(W[p][i]*concentrations[q][i] + W[q][i]*concentrations[p][i]))/elementSize10368;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] -= (Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1]);
			}
		}
	}
}
//---------------------------------------------------------------------------
void MigrationTerm_AX_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);
	
	// Calculate coefficients
	double R[3][3][3][2];
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
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		for (unsigned c=0; c<nDimensions; c++)
		{
			R[m][m][m][c] = 195.*normals[m][c]*coordinates[m][1]
							+ (109.*normals[m][c] - 43.*normals[p][c])*coordinates[p][1]
							+ (109.*normals[m][c] - 43.*normals[q][c])*coordinates[q][1];
			R[m][p][p][c] = (89.*normals[m][c] - 53.*normals[p][c])*coordinates[m][1]
							+ (81.*normals[m][c] - 57.*normals[p][c])*coordinates[p][1]
							+ (46.*normals[m][c] - 10.*normals[p][c])*coordinates[q][1];
			R[m][q][q][c] = (89.*normals[m][c] - 53.*normals[q][c])*coordinates[m][1]
							+ (46.*normals[m][c] - 10.*normals[q][c])*coordinates[p][1]
							+ (81.*normals[m][c] - 57.*normals[q][c])*coordinates[q][1];
			R[m][m][p][c] = (109.*normals[m][c] - 43.*normals[p][c])*coordinates[m][1]
							+ (89.*normals[m][c] - 53.*normals[p][c])*coordinates[p][1]
							+ 66.*normals[m][c]*coordinates[q][1];
			R[m][p][m][c] = R[m][m][p][c];
			R[m][m][q][c] = (109.*normals[m][c] - 43.*normals[q][c])*coordinates[m][1]
							+ 66.*normals[m][c]*coordinates[p][1]
							+ (89.*normals[m][c] - 53.*normals[q][c])*coordinates[q][1];
			R[m][q][m][c] = R[m][m][p][c];
			R[m][p][q][c] = 66.*normals[m][c]*coordinates[m][1]
							+ (46.*normals[m][c] - 10.*normals[p][c])*coordinates[p][1]
							+ (46.*normals[m][c] - 10.*normals[q][c])*coordinates[q][1];
			R[m][q][p][c] = R[m][p][q][c];
		}
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
				Wmi[m][c] = (R[m][m][m][c]*W[m][i] + R[m][m][p][c]*W[p][i] + R[m][m][q][c]*W[q][i])/5184.;
				Wmi[p][c] = (R[m][p][m][c]*W[m][i] + R[m][p][p][c]*W[p][i] + R[m][p][q][c]*W[q][i])/5184.;
				Wmi[q][c] = (R[m][q][m][c]*W[m][i] + R[m][q][p][c]*W[p][i] + R[m][q][q][c]*W[q][i])/5184.;
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0] + gradU[1]*Wmi[n][1];
			}
		}
	}
}
//---------------------------------------------------------------------------

