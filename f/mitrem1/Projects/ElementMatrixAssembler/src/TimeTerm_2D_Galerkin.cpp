//---------------------------------------------------------------------------

#include "TimeTerm_2D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_2D_Galerkin::TimeTerm_2D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_2D_Galerkin::~TimeTerm_2D_Galerkin()
{
}
//---------------------------------------------------------------------------
void TimeTerm_2D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;

  // Calculate coefficients
  elementSize = elementProps->calcSize(coordinates);
  double elementSize12 = elementSize/12.;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned p = (m+1)%3;
    unsigned q = (m+2)%3;
    for (unsigned i=0; i<nIons; i++)
    {
      unsigned eqmi = eq(m,i);
      elementMat[eqmi][var(m,i)] += 2.*elementSize12*(1.-volumeGasFraction);
      elementMat[eqmi][var(p,i)] +=    elementSize12*(1.-volumeGasFraction);
      elementMat[eqmi][var(q,i)] +=    elementSize12*(1.-volumeGasFraction);
    }
  }
}
//---------------------------------------------------------------------------
void TimeTerm_2D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

