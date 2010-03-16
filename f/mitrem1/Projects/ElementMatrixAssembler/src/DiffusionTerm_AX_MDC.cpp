//---------------------------------------------------------------------------

#include "DiffusionTerm_AX_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
DiffusionTerm_AX_MDC::DiffusionTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: DiffusionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
DiffusionTerm_AX_MDC::~DiffusionTerm_AX_MDC()
{
}
//---------------------------------------------------------------------------
void DiffusionTerm_AX_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	// CHECK THIS, because it is probably incorrect for axisymmetric cases
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	double Bruggeman = pow(1.-volumeGasFraction,1.5);
	
	// Calculate coefficients
	double R[3][3][2];
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
			for (unsigned j=0; j<nIons; j++) 
			{
				D[m][i][j] = mitrem->calcTransportDiffusionFactor(i,j)*Bruggeman;
			}
		}
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		for (unsigned c=0; c<nDimensions; c++)
		{
			R[m][m][c] = 19.*normals[m][c]*coordinates[m][1]
						 + (11.*normals[m][c] - 4.*normals[p][c])*coordinates[p][1]
						 + (11.*normals[m][c] - 4.*normals[q][c])*coordinates[q][1];
			R[m][p][c] = (11.*normals[m][c] - 4.*normals[p][c])*coordinates[m][1]
						 + (9.*normals[m][c] - 5.*normals[p][c])*coordinates[p][1]
						 + 7.*normals[m][c]*coordinates[q][1];
			R[m][q][c] = (11.*normals[m][c] - 4.*normals[q][c])*coordinates[m][1]
						 + 7.*normals[m][c]*coordinates[p][1]
						 + (9.*normals[m][c] - 5.*normals[q][c])*coordinates[q][1];
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize432 = 432.*elementSize;
	
	// Add to element matrix
	for (unsigned i=0; i<nIons; i++) 
	{
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned p = (m+1)%3;
			unsigned q = (m+2)%3;
			unsigned eqmi = eq(m,i);
			for (unsigned j=0; j<nIons; j++) 
			{
				for (unsigned c=0; c<nDimensions; c++)
				{		
					Dmij[c] = (R[m][m][c]*D[m][i][j] + R[m][p][c]*D[p][i][j] + R[m][q][c]*D[q][i][j])/elementSize432;
				}
				for (unsigned n=0; n<nNodes; n++) 
				{
					elementMat[eqmi][var(n,j)] -= (Dmij[0]*normals[n][0] + Dmij[1]*normals[n][1]);
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
void DiffusionTerm_AX_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

