//---------------------------------------------------------------------------

#include "ConvectionTerm_1D_N.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ConvectionTerm_1D_N::ConvectionTerm_1D_N(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ConvectionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ConvectionTerm_1D_N::~ConvectionTerm_1D_N()
{
}
//---------------------------------------------------------------------------
void ConvectionTerm_1D_N::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.5;

  // Calculate coefficients
  unsigned nTargetNodes = 0;
  unsigned targetNode = 0;
  //unsigned nOtherNodes = 0;
  //unsigned otherNodes[2];
  for (unsigned d=0; d<nDimensions; d++)
  {
    averageVelocity[d] = 0;
    for (unsigned m=0; m<nNodes; m++)
    {
      averageVelocity[d] += velocities[m][d];
    }
    averageVelocity[d] *= 0.5;
  }
  for (unsigned m=0; m<nNodes; m++)
  {
    normal = elementProps->calcNormal(m,coordinates);
    k[m] = normal[0]*averageVelocity[0];
    if (k[m] > 0.)
    {
      targetNode = m;
      nTargetNodes++;
    }
    //else
    //{
    //  otherNodes[nOtherNodes] = m;
    //  nOtherNodes++;
    //}
  }

  // Add to element matrix
  if (nTargetNodes == 1)
  {
    for (unsigned i=0; i<nIons; i++)
    {
      for (unsigned m=0; m<nNodes; m++)
      {
        elementMat[eq(targetNode,i)][var(m,i)] -= k[m]*(1.-volumeGasFraction);
      }
    }
  }
}
//---------------------------------------------------------------------------
void ConvectionTerm_1D_N::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

