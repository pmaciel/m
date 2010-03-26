//---------------------------------------------------------------------------

#include "ElectrostaticsTerm_1D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrostaticsTerm_1D_MDC::ElectrostaticsTerm_1D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ElectrostaticsTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ElectrostaticsTerm_1D_MDC::~ElectrostaticsTerm_1D_MDC()
{
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_1D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
  // Calculate coefficients
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    normal = elementProps->calcNormal(m,coordinates);
    for (unsigned c=0; c<nDimensions; c++)
    {
      normals[m][c] = normal[c];
    }
  }
  K = mitrem->calcElectrostaticsPotentialFactor();
  for (unsigned i=0; i<nIons; i++)
  {
    Z[i] = mitrem->calcElectrostaticsConcentrationFactor(i);
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize8 = elementSize/8.;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned eqmnIons = eq(m,nIons);
    unsigned p = (m+1)%2;
    for (unsigned j=0; j<nIons; j++)
    {
      elementMat[eqmnIons][var(m,j)] += 3.*Z[j]*elementSize8;
      elementMat[eqmnIons][var(p,j)] +=    Z[j]*elementSize8;
    }
    for (unsigned n=0; n<nNodes; n++)
    {
      double normalProductmn = normals[m][0]*normals[n][0];
      elementMat[eqmnIons][var(n,nIons)] -= normalProductmn*K/elementSize;
    }
  }
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_1D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

