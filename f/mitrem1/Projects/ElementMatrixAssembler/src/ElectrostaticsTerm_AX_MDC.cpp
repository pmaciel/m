//---------------------------------------------------------------------------

#include "ElectrostaticsTerm_AX_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrostaticsTerm_AX_MDC::ElectrostaticsTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ElectrostaticsTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ElectrostaticsTerm_AX_MDC::~ElectrostaticsTerm_AX_MDC()
{
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_AX_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
	// Calculate coefficients	
	double Ru[3][2];
	double Rz[3][3];
	for (unsigned m=0; m<nNodes; m++)
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		normal = elementProps->calcNormal(m,coordinates);
		for (unsigned c=0; c<nDimensions; c++)
		{
			normals[m][c] = normal[c];
		}
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		Rz[m][m] = 170.*coordinates[m][1] + 47.*coordinates[p][1] + 47.*coordinates[q][1];
		Rz[m][p] =  47.*coordinates[m][1] + 23.*coordinates[p][1] + 14.*coordinates[q][1];
		Rz[m][q] =  47.*coordinates[m][1] + 14.*coordinates[p][1] + 23.*coordinates[q][1];
	}
	for (unsigned m=0; m<nNodes; m++)
	{
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		for (unsigned c=0; c<nDimensions; c++)
		{
			Ru[m][c] = 5.*normals[m][c]*coordinates[m][1]
					   + (3.*normals[m][c] - normals[p][c])*coordinates[p][1] 
					   + (3.*normals[m][c] - normals[q][c])*coordinates[q][1];
		}
	}
	K = mitrem->calcElectrostaticsPotentialFactor();
	for (unsigned i=0; i<nIons; i++) 
	{
		Z[i] = mitrem->calcElectrostaticsConcentrationFactor(i);
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize1296 = elementSize/1296.;
	double elementSize48 = 48.*elementSize;

	// Add to element matrix
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned eqmnIons = eq(m,nIons);
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		for (unsigned j=0; j<nIons; j++) 
		{
			elementMat[eqmnIons][var(m,j)] += Rz[m][m]*Z[j]*elementSize1296;
			elementMat[eqmnIons][var(p,j)] += Rz[m][p]*Z[j]*elementSize1296;
			elementMat[eqmnIons][var(q,j)] += Rz[m][q]*Z[j]*elementSize1296;
		}
		for (unsigned n=0; n<nNodes; n++) 
		{
			elementMat[eqmnIons][var(n,nIons)] -= (Ru[m][0]*normals[n][0] + Ru[m][1]*normals[n][1])*K/elementSize48;
		}
	}
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_AX_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

