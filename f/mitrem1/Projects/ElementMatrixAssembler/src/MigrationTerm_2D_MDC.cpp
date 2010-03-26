//---------------------------------------------------------------------------

#include "MigrationTerm_2D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_2D_MDC::MigrationTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_2D_MDC::~MigrationTerm_2D_MDC()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_2D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;
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
  double elementSize432 = 432.*elementSize;

  // Add to element matrix
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned eqmi = eq(m,i);
      unsigned p = (m+1)%3;
      unsigned q = (m+2)%3;
      for (unsigned c=0; c<nDimensions; c++)
      {
        Mmi[c] = (19.*normals[m][c]*W[m][i]*concentrations[m][i]
             + (9.*normals[m][c] - 5.*normals[p][c])*W[p][i]*concentrations[p][i]
             + (9.*normals[m][c] - 5.*normals[q][c])*W[q][i]*concentrations[q][i]
             + (11.*normals[m][c] - 4.*normals[p][c])*(W[m][i]*concentrations[p][i] + W[p][i]*concentrations[m][i])
             + (11.*normals[m][c] - 4.*normals[q][c])*(W[m][i]*concentrations[q][i] + W[q][i]*concentrations[m][i])
             + 7.*normals[m][c]*(W[p][i]*concentrations[q][i] + W[q][i]*concentrations[p][i]))/elementSize432;
      }
      for (unsigned n=0; n<nNodes; n++)
      {
        elementMat[eqmi][var(n,nIons)] -= (Mmi[0]*normals[n][0] + Mmi[1]*normals[n][1]);
      }
    }
  }
}
//---------------------------------------------------------------------------
void MigrationTerm_2D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;
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
  double elementSize2 = 2.*elementSize;
  for (unsigned c=0; c<nDimensions; c++)
  {
    gradU[c] /= elementSize2;
  }

  // Add to element jacobian
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned eqmi = eq(m,i);
      unsigned p = (m+1)%3;
      unsigned q = (m+2)%3;
      for (unsigned c=0; c<nDimensions; c++)
      {
        Wmi[m][c] = (19.*normals[m][c]*W[m][i]
              + (11.*normals[m][c] - 4.*normals[p][c])*W[p][i]
              + (11.*normals[m][c] - 4.*normals[q][c])*W[q][i])/216.;
        Wmi[p][c] = ((11.*normals[m][c] - 4.*normals[p][c])*W[m][i]
              + (9.*normals[m][c] - 5.*normals[p][c])*W[p][i]
              + 7.*normals[m][c]*W[q][i])/216.;
        Wmi[q][c] = ((11.*normals[m][c] - 4.*normals[q][c])*W[m][i]
              + 7.*normals[m][c]*W[p][i]
              + (9.*normals[m][c] - 5.*normals[q][c])*W[q][i])/216.;
      }
      for (unsigned n=0; n<nNodes; n++)
      {
        elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0] + gradU[1]*Wmi[n][1];
      }
    }
  }
}
//---------------------------------------------------------------------------

