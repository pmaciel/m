//---------------------------------------------------------------------------

#include "TimeTerm_2D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_2D_MDC::TimeTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_2D_MDC::~TimeTerm_2D_MDC()
{
}
//---------------------------------------------------------------------------
void TimeTerm_2D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;

  // Calculate coefficients
  elementSize = elementProps->calcSize(coordinates);
  double elementSize108 = elementSize/108.;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned p = (m+1)%3;
    unsigned q = (m+2)%3;
    for (unsigned i=0; i<nIons; i++)
    {
      unsigned eqmi = eq(m,i);
      elementMat[eqmi][var(m,i)] += 22.*elementSize108*(1.-volumeGasFraction);
      elementMat[eqmi][var(p,i)] +=  7.*elementSize108*(1.-volumeGasFraction);
      elementMat[eqmi][var(q,i)] +=  7.*elementSize108*(1.-volumeGasFraction);
    }
  }
}
//---------------------------------------------------------------------------
void TimeTerm_2D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

