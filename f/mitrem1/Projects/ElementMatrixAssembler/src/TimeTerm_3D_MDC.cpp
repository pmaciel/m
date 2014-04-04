//---------------------------------------------------------------------------

#include "TimeTerm_3D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_3D_MDC::TimeTerm_3D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_3D_MDC::~TimeTerm_3D_MDC()
{
}
//---------------------------------------------------------------------------
void TimeTerm_3D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.25;

  // Calculate coefficients
  elementSize = elementProps->calcSize(coordinates);
  double elementSize576 = elementSize/576.;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned p = (m+1)%4;
    unsigned q = (m+2)%4;
    unsigned r = (m+3)%4;
    for (unsigned i=0; i<nIons; i++)
    {
      unsigned eqmi = eq(m,i);
      elementMat[eqmi][var(m,i)] += 75.*elementSize576*(1.-volumeGasFraction);
      elementMat[eqmi][var(p,i)] += 23.*elementSize576*(1.-volumeGasFraction);
      elementMat[eqmi][var(q,i)] += 23.*elementSize576*(1.-volumeGasFraction);
      elementMat[eqmi][var(r,i)] += 23.*elementSize576*(1.-volumeGasFraction);
    }
  }
}
//---------------------------------------------------------------------------
void TimeTerm_3D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

