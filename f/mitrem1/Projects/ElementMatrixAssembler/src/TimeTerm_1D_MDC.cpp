//---------------------------------------------------------------------------

#include "TimeTerm_1D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_1D_MDC::TimeTerm_1D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_1D_MDC::~TimeTerm_1D_MDC()
{
}
//---------------------------------------------------------------------------
void TimeTerm_1D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.5;

  // Calculate coefficients
  elementSize = elementProps->calcSize(coordinates);
  double elementSize8 = elementSize/8.;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned p = (m+1)%2;
    for (unsigned i=0; i<nIons; i++)
    {
      unsigned eqmi = eq(m,i);
      elementMat[eqmi][var(m,i)] += 3.*elementSize8*(1.-volumeGasFraction);
      elementMat[eqmi][var(p,i)] +=    elementSize8*(1.-volumeGasFraction);
    }
  }
}
//---------------------------------------------------------------------------
void TimeTerm_1D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

