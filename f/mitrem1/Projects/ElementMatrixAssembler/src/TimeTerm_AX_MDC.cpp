//---------------------------------------------------------------------------

#include "TimeTerm_AX_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_AX_MDC::TimeTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_AX_MDC::~TimeTerm_AX_MDC()
{
}
//---------------------------------------------------------------------------
void TimeTerm_AX_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
	double elementSize1296 = elementSize/1296.;
	
	// Add to element matrix
	for (unsigned m=0; m<nNodes; m++) 
	{
		unsigned p = (m+1)%3;
		unsigned q = (m+2)%3;
		R[m][m] = 170.*coordinates[m][1] + 47.*coordinates[p][1] + 47.*coordinates[q][1];
		R[m][p] =  47.*coordinates[m][1] + 23.*coordinates[p][1] + 14.*coordinates[q][1];
		R[m][q] =  47.*coordinates[m][1] + 14.*coordinates[p][1] + 23.*coordinates[q][1];
		for (unsigned i=0; i<nIons; i++) 
		{
			unsigned eqmi = eq(m,i);
			elementMat[eqmi][var(m,i)] += R[m][m]*elementSize1296*(1.-volumeGasFraction);
			elementMat[eqmi][var(p,i)] += R[m][p]*elementSize1296*(1.-volumeGasFraction);
			elementMat[eqmi][var(q,i)] += R[m][q]*elementSize1296*(1.-volumeGasFraction);
		}
	}
}
//---------------------------------------------------------------------------
void TimeTerm_AX_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

