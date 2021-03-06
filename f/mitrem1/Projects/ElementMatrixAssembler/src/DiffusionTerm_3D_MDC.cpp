//---------------------------------------------------------------------------

#include "DiffusionTerm_3D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
DiffusionTerm_3D_MDC::DiffusionTerm_3D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : DiffusionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
DiffusionTerm_3D_MDC::~DiffusionTerm_3D_MDC()
{
}
//---------------------------------------------------------------------------
void DiffusionTerm_3D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.25;
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
  double elementSize324 = 324.*elementSize;

  // Add to element matrix
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned p = (m+1)%4;
      unsigned q = (m+2)%4;
      unsigned r = (m+3)%4;
      unsigned eqmi = eq(m,i);
      for (unsigned j=0; j<nIons; j++)
      {
        for (unsigned c=0; c<nDimensions; c++)
        {
          Dmij[c] = (13.*normals[m][c]*D[m][i][j]
                + (7.*normals[m][c] - 2.*normals[p][c])*D[p][i][j]
                + (7.*normals[m][c] - 2.*normals[q][c])*D[q][i][j]
                + (7.*normals[m][c] - 2.*normals[r][c])*D[r][i][j])/elementSize324;
        }
        for (unsigned n=0; n<nNodes; n++)
        {
          elementMat[eqmi][var(n,j)] -= (Dmij[0]*normals[n][0] + Dmij[1]*normals[n][1] + Dmij[2]*normals[n][2]);
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void DiffusionTerm_3D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

