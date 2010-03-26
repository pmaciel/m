//---------------------------------------------------------------------------

#include "TimeTerm_3D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_3D_Galerkin::TimeTerm_3D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
TimeTerm_3D_Galerkin::~TimeTerm_3D_Galerkin()
{
}
//---------------------------------------------------------------------------
void TimeTerm_3D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.25;

  // Calculate coefficients
  elementSize = elementProps->calcSize(coordinates);
  double elementSize20 = elementSize/20.;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned p = (m+1)%4;
    unsigned q = (m+2)%4;
    unsigned r = (m+3)%4;
    for (unsigned i=0; i<nIons; i++)
    {
      unsigned eqmi = eq(m,i);
      elementMat[eqmi][var(m,i)] += 2*elementSize20*(1.-volumeGasFraction);
      elementMat[eqmi][var(p,i)] +=   elementSize20*(1.-volumeGasFraction);
      elementMat[eqmi][var(q,i)] +=   elementSize20*(1.-volumeGasFraction);
      elementMat[eqmi][var(r,i)] +=   elementSize20*(1.-volumeGasFraction);
    }
  }
}
//---------------------------------------------------------------------------
void TimeTerm_3D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

