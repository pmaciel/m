//---------------------------------------------------------------------------

#include "ElectrostaticsTerm_AX_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrostaticsTerm_AX_Galerkin::ElectrostaticsTerm_AX_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ElectrostaticsTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ElectrostaticsTerm_AX_Galerkin::~ElectrostaticsTerm_AX_Galerkin()
{
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_AX_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{	
	// Calculate coefficients	
	double Ru;
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
		Rz[m][m] = 6.*coordinates[m][0] + 2.*coordinates[p][0] + 2.*coordinates[q][0];
		Rz[m][p] = 2.*coordinates[m][0] + 2.*coordinates[p][0] +    coordinates[q][0];
		Rz[m][q] = 2.*coordinates[m][0] +    coordinates[p][0] + 2.*coordinates[q][0];
	}
	Ru = coordinates[0][0] + coordinates[1][0] + coordinates[2][0];
	K = mitrem->calcElectrostaticsPotentialFactor();
	for (unsigned i=0; i<nIons; i++) 
	{
		Z[i] = mitrem->calcElectrostaticsConcentrationFactor(i);
	}
	elementSize = elementProps->calcSize(coordinates);
	double elementSize60 = elementSize/60.;
	double elementSize12 = 12.*elementSize;

	// Add to element matrix
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned eqmnIons = eq(m,nIons);
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		for (unsigned j=0; j<nIons; j++) 
		{
			elementMat[eqmnIons][var(m,j)] += Rz[m][m]*Z[j]*elementSize60;
			elementMat[eqmnIons][var(p,j)] += Rz[m][p]*Z[j]*elementSize60;
			elementMat[eqmnIons][var(q,j)] += Rz[m][q]*Z[j]*elementSize60;
		}
		for (unsigned n=0; n<nNodes; n++) 
		{
			double normalProductmn = normals[m][0]*normals[n][0] + normals[m][1]*normals[n][1];
			elementMat[eqmnIons][var(n,nIons)] -= normalProductmn*K*Ru/elementSize12;
		}
	}
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_AX_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

