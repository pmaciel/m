//---------------------------------------------------------------------------

#include "TimeTerm_AX_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_AX_Galerkin::TimeTerm_AX_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_AX_Galerkin::~TimeTerm_AX_Galerkin()
{
}
//---------------------------------------------------------------------------
void TimeTerm_AX_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	double volumeGasFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		volumeGasFraction += volumeGasFractions[m];
	}
	volumeGasFraction /= 3.;
	
	// Calculate coefficients
	double R[3][3];
	elementSize = elementProps->calcSize(coordinates);
	double elementSize60 = elementSize/60.;
	
	// Add to element matrix
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		R[m][m] = 6.*coordinates[m][0] + 2.*coordinates[p][0] + 2.*coordinates[q][0];
		R[m][p] = 2.*coordinates[m][0] + 2.*coordinates[p][0] +    coordinates[q][0];
		R[m][q] = 2.*coordinates[m][0] +    coordinates[p][0] + 2.*coordinates[q][0];
		for (unsigned i=0; i<nIons; i++) 
		{
			unsigned eqmi = eq(m,i);
			elementMat[eqmi][var(m,i)] += R[m][m]*elementSize60*(1.-volumeGasFraction);
			elementMat[eqmi][var(p,i)] += R[m][p]*elementSize60*(1.-volumeGasFraction);
			elementMat[eqmi][var(q,i)] += R[m][q]*elementSize60*(1.-volumeGasFraction);
		}
	}
}
//---------------------------------------------------------------------------
void TimeTerm_AX_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

