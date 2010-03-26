//---------------------------------------------------------------------------

#include "DiffusionTerm_1D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
DiffusionTerm_1D_MDC::DiffusionTerm_1D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : DiffusionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
DiffusionTerm_1D_MDC::~DiffusionTerm_1D_MDC()
{
}
//---------------------------------------------------------------------------
void DiffusionTerm_1D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.5;
  double Bruggeman = pow(1.-volumeGasFraction,1.5);

  // Calculate coefficients
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    normal = elementProps->calcNormal(m,coordinates);
    for (unsigned c=0; c<nDimensions; c++)
    {
      normals[m][c] = normal[c];
    }
    for (unsigned i=0; i<nIons; i++)
    {
      for (unsigned j=0; j<nIons; j++)
      {
        D[m][i][j] = mitrem->calcTransportDiffusionFactor(i,j)*Bruggeman;
      }
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize2 = 2.*elementSize;

  // Add to element matrix
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned eqmi = eq(m,i);
      for (unsigned j=0; j<nIons; j++)
      {
        Dmij[0] = (D[0][i][j] + D[1][i][j])*normals[m][0]/elementSize2;
        for (unsigned n=0; n<nNodes; n++)
        {
          elementMat[eqmi][var(n,j)] -= Dmij[0]*normals[n][0];
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void DiffusionTerm_1D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

