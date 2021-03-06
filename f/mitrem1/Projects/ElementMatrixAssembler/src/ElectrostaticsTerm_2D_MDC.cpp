//---------------------------------------------------------------------------

#include "ElectrostaticsTerm_2D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrostaticsTerm_2D_MDC::ElectrostaticsTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ElectrostaticsTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ElectrostaticsTerm_2D_MDC::~ElectrostaticsTerm_2D_MDC()
{
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_2D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
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
  double elementSize108 = elementSize/108.;
  double elementSize4 = 4.*elementSize;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned eqmnIons = eq(m,nIons);
    unsigned p = (m+1)%3;
    unsigned q = (m+2)%3;
    for (unsigned j=0; j<nIons; j++)
    {
      elementMat[eqmnIons][var(m,j)] += 22.*Z[j]*elementSize108;
      elementMat[eqmnIons][var(p,j)] +=  7.*Z[j]*elementSize108;
      elementMat[eqmnIons][var(q,j)] +=  7.*Z[j]*elementSize108;
    }
    for (unsigned n=0; n<nNodes; n++)
    {
      double normalProductmn = normals[m][0]*normals[n][0] + normals[m][1]*normals[n][1];
      elementMat[eqmnIons][var(n,nIons)] -= normalProductmn*K/elementSize4;
    }
  }
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_2D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

