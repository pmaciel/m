//---------------------------------------------------------------------------

#include "ElectrostaticsTerm_3D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrostaticsTerm_3D_MDC::ElectrostaticsTerm_3D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ElectrostaticsTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
ElectrostaticsTerm_3D_MDC::~ElectrostaticsTerm_3D_MDC()
{
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_3D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
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
  double elementSize576 = elementSize/576.;
  double elementSize12 = 12.*elementSize;

  // Add to element matrix
  for (unsigned m=0; m<nNodes; m++)
  {
    unsigned eqmnIons = eq(m,nIons);
    unsigned p = (m+1)%4;
    unsigned q = (m+2)%4;
    unsigned r = (m+3)%4;
    for (unsigned j=0; j<nIons; j++)
    {
      elementMat[eqmnIons][var(m,j)] += 75.*Z[j]*elementSize576;
      elementMat[eqmnIons][var(p,j)] += 23.*Z[j]*elementSize576;
      elementMat[eqmnIons][var(q,j)] += 23.*Z[j]*elementSize576;
      elementMat[eqmnIons][var(r,j)] += 23.*Z[j]*elementSize576;
    }
    for (unsigned n=0; n<nNodes; n++)
    {
      double normalProductmn = normals[m][0]*normals[n][0] + normals[m][1]*normals[n][1] + normals[m][2]*normals[n][2];
      elementMat[eqmnIons][var(n,nIons)] -= normalProductmn*K/elementSize12;
    }
  }
}
//---------------------------------------------------------------------------
void ElectrostaticsTerm_3D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

