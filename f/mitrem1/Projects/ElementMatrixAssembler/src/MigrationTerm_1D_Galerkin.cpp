//---------------------------------------------------------------------------

#include "MigrationTerm_1D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_1D_Galerkin::MigrationTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_1D_Galerkin::~MigrationTerm_1D_Galerkin()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_1D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
  double elementSize6 = 6.*elementSize;

  // Add to element matrix
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned eqmi = eq(m,i);
      Mmi[0] = ((2.*W[0][i] + W[1][i])*concentrations[0][i] + (W[0][i] + 2.*W[1][i])*concentrations[1][i])*normals[m][0]/elementSize6;
      for (unsigned n=0; n<nNodes; n++)
      {
        elementMat[eqmi][var(n,nIons)] -= Mmi[0]*normals[n][0];
      }
    }
  }
}
//---------------------------------------------------------------------------
void MigrationTerm_1D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
      gradU[c] += normal[c]*potentials[m];
    }
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned i=0; i<nIons; i++)
    {
      W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  //double elementSize6 = 6.*elementSize;
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
      unsigned p = (m+1)%2;
      Wmi[m][0] = (2.*W[m][i] + W[p][i])*normals[m][0]/6.;
      Wmi[p][0] = (W[m][i] + 2.*W[p][i])*normals[m][0]/6.;
      for (unsigned n=0; n<nNodes; n++)
      {
        elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0];
      }
    }
  }
}
//---------------------------------------------------------------------------

