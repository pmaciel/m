//---------------------------------------------------------------------------

#include "DiffusionTerm_2D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
DiffusionTerm_2D_MDC::DiffusionTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: DiffusionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
DiffusionTerm_2D_MDC::~DiffusionTerm_2D_MDC()
{
}
//---------------------------------------------------------------------------
void DiffusionTerm_2D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
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
			for (unsigned j=0; j<nIons; j++) 
			{
				D[m][i][j] = mitrem->calcTransportDiffusionFactor(i,j)*Bruggeman;
			}
		}
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize48 = 48.*elementSize;
	
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
					Dmij[c] = (5.*normals[m][c]*D[m][i][j]
							  + (3.*normals[m][c] - normals[p][c])*D[p][i][j]
							  + (3.*normals[m][c] - normals[q][c])*D[q][i][j])/elementSize48;
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
void DiffusionTerm_2D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

