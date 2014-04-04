//---------------------------------------------------------------------------

#include "ConvectionTerm_2D_N_Migration.h"

#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ConvectionTerm_2D_N_Migration::ConvectionTerm_2D_N_Migration(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ConvectionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
  gradc = new double*[nIons];
  for (unsigned i=0; i<nIons; i++)
  {
    gradc[i] = new double[nDimensions];
  }
  gradU = new double[nDimensions];
}
//---------------------------------------------------------------------------
ConvectionTerm_2D_N_Migration::~ConvectionTerm_2D_N_Migration()
{
  for (unsigned i=0; i<nIons; i++)
  {
    delete[] gradc[i];
  }
  delete[] gradc;
  delete[] gradU;
}
//---------------------------------------------------------------------------
void ConvectionTerm_2D_N_Migration::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
      gradU[c] += normal[c]*potentials[m];
    }
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
  }
  elementSize = elementProps->calcSize(coordinates);
  for (unsigned c=0; c<nDimensions; c++)
  {
    gradU[c] /= nDimensions*elementSize;
  }

  unsigned nTargetNodes;
  unsigned targetNodes[2];
  unsigned nOtherNodes;
  unsigned otherNodes[3];

  for (unsigned i=0; i<nIons; i++)
  {
    nTargetNodes = 0;
    nOtherNodes = 0;

    for (unsigned d=0; d<nDimensions; d++)
    {
      averageVelocity[d] = 0;
      for (unsigned m=0; m<nNodes; m++)
      {
        averageVelocity[d] += velocities[m][d]*(1.-volumeGasFraction);
      }
      averageVelocity[d] /= 3.;
      averageVelocity[d] -= mitrem->calcTransportMigrationFactor(i)*Bruggeman*gradU[d];
    }
    for (unsigned m=0; m<nNodes; m++)
    {
      normal = elementProps->calcNormal(m,coordinates);
      k[m] = 0.5*(normal[0]*averageVelocity[0] + normal[1]*averageVelocity[1]);
      if (k[m] > 0.)
      {
        targetNodes[nTargetNodes] = m;
        nTargetNodes++;
      }
      else
      {
        otherNodes[nOtherNodes] = m;
        nOtherNodes++;
      }
    }

  // Add to element matrix
    if (nTargetNodes == 1)
    {
      for (unsigned m=0; m<nNodes; m++)
      {
        elementMat[eq(targetNodes[0],i)][var(m,i)] -= k[m];
      }
    }
    else if (nTargetNodes == 2)
    {
      elementMat[eq(targetNodes[0],i)][var(targetNodes[0],i)] -= k[targetNodes[0]];
      elementMat[eq(targetNodes[0],i)][var(otherNodes[0],i)]  += k[targetNodes[0]];
      elementMat[eq(targetNodes[1],i)][var(targetNodes[1],i)] -= k[targetNodes[1]];
      elementMat[eq(targetNodes[1],i)][var(otherNodes[0],i)]  += k[targetNodes[1]];
    }
  }
}
//---------------------------------------------------------------------------
void ConvectionTerm_2D_N_Migration::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;
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

  unsigned nTargetNodes;
  unsigned targetNodes[2];
  unsigned nOtherNodes;
  unsigned otherNodes[3];
  for (unsigned i=0; i<nIons; i++)
  {
    nTargetNodes = 0;
    nOtherNodes = 0;
    for (unsigned d=0; d<nDimensions; d++)
    {
      averageVelocity[d] = 0;
      averageVelocity[d] -= mitrem->calcTransportMigrationFactor(i)*Bruggeman*gradc[i][d];
    }
    for (unsigned m=0; m<nNodes; m++)
    {
      normal = elementProps->calcNormal(m,coordinates);
      k[m] = 0.5*(normal[0]*averageVelocity[0] + normal[1]*averageVelocity[1]);
      if (k[m] > 0.)
      {
        targetNodes[nTargetNodes] = m;
        nTargetNodes++;
      }
      else
      {
        otherNodes[nOtherNodes] = m;
        nOtherNodes++;
      }
    }

    // Add to element matrix
    if (nTargetNodes == 1)
    {
      for (unsigned m=0; m<nNodes; m++)
      {
        elementJac[eq(targetNodes[0],i)][var(m,nIons)] -= k[m];
      }
    }
    else if (nTargetNodes == 2)
    {
      elementJac[eq(targetNodes[0],i)][var(targetNodes[0],nIons)] -= k[targetNodes[0]];
      elementJac[eq(targetNodes[0],i)][var(otherNodes[0],nIons)]  += k[targetNodes[0]];
      elementJac[eq(targetNodes[1],i)][var(targetNodes[1],nIons)] -= k[targetNodes[1]];
      elementJac[eq(targetNodes[1],i)][var(otherNodes[0],nIons)]  += k[targetNodes[1]];
    }
  }
}
//---------------------------------------------------------------------------

