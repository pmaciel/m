//---------------------------------------------------------------------------

#include "MigrationTerm_1D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_1D_MDC::MigrationTerm_1D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_1D_MDC::~MigrationTerm_1D_MDC()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_1D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
      W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize4 = 4.*elementSize;

  // Add to element matrix
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned eqmi = eq(m,i);
      Mmi[0] = (W[0][i] + W[1][i])*(concentrations[0][i] + concentrations[1][i])*normals[m][0]/elementSize4;
      for (unsigned n=0; n<nNodes; n++)
      {
        elementMat[eqmi][var(n,nIons)] -= Mmi[0]*normals[n][0];
      }
    }
  }
}
//---------------------------------------------------------------------------
void MigrationTerm_1D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.5;
  double Bruggeman = pow(1.-volumeGasFraction,1.5);

  // Calculate coefficients
  for (unsigned c=0; c<nDimensions; c++)
  {
    gradU[c] = 0;
  }
  for (unsigned m=0; m<nNodes; m++)
  {
    normal = elementProps->calcNormal(m,coordinates);
    for (unsigned c=0; c<nDimensions; c++)
    {
      normals[m][c] = normal[c];
      gradU[c] += normal[c]*potentials[m];
    }
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned i=0; i<nIons; i++)
    {
      W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  for (unsigned c=0; c<nDimensions; c++)
  {
    gradU[c] /= elementSize;
  }

  // Add to element jacobian
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned eqmi = eq(m,i);
      Wmi[0][0] = 0.25*(W[0][i] + W[1][i])*normals[m][0]; // equals Wmi[1][0]
      for (unsigned n=0; n<nNodes; n++)
      {
        elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[0][0];
      }
    }
  }
}
//---------------------------------------------------------------------------

