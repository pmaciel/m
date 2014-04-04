//---------------------------------------------------------------------------

#include "DiffusionTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
DiffusionTerm::DiffusionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ElementTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
  D = new double**[nNodes];
  for (unsigned m=0; m<nNodes; m++)
  {
    D[m] = new double*[nIons];
    for (unsigned i=0; i<nIons; i++)
    {
      D[m][i] = new double[nIons];
    }
  }
  Dmij = new double[nDimensions];
  gradc = new double*[nIons];
  for (unsigned i=0; i<nIons; i++)
  {
    gradc[i] = new double[nDimensions];
  }
}
//---------------------------------------------------------------------------
DiffusionTerm::~DiffusionTerm()
{
  for (unsigned m=0; m<nNodes; m++)
  {
    for (unsigned i=0; i<nIons; i++)
    {
      delete[] D[m][i];
    }
    delete[] D[m];
  }
  delete[] D;
  delete[] Dmij;
  for (unsigned i=0; i<nIons; i++)
  {
    delete[] gradc[i];
  }
  delete[] gradc;
}
//---------------------------------------------------------------------------
void DiffusionTerm::calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= nNodes;
  double Bruggeman = pow(1.-volumeGasFraction,1.5);

  // Calculate coefficients
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      gradc[i][c] = 0.;
    }
  }
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    normal = elementProps->calcNormal(m,coordinates);
    for (unsigned i=0; i<nIons; i++)
    {
      for (unsigned c=0; c<nDimensions; c++)
      {
        gradc[i][c] += normal[c]*concentrations[m][i];
      }
      for (unsigned j=0; j<nIons; j++)
      {
        D[m][i][j] = mitrem->calcTransportDiffusionFactor(i,j)*Bruggeman;
      }
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      gradc[i][c] /= nDimensions*elementSize;
    }
  }

  // Add to element current density
  for (unsigned i=0; i<nIons; i++)
  {
    double ziF = F_CONST*mitrem->getIonChargeNumber(i);
    for (unsigned j=0; j<nIons; j++)
    {
      double Dij = 0.;
      for (unsigned m=0; m<nNodes; m++)
      {
        Dij += D[m][i][j];
      }
      Dij /= nNodes;
      for (unsigned c=0; c<nDimensions; c++)
      {
        elementCurr[i][c] -= ziF*Dij*gradc[j][c];
      }
    }
  }
}
//---------------------------------------------------------------------------

