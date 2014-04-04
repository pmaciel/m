//---------------------------------------------------------------------------

#include "ElecReactionTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReactionTerm::ElecReactionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_)
  : BoundaryElementTerm(nDimensions_, nNodes_, nVariables_, mitrem_, boundaryElementProps_)
{
  unsigned nElecReactions = mitrem->getNElecReactions();
  v = new double*[nNodes];
  DvDU = new double*[nNodes];
  DvDCRed = new double**[nNodes];
  DvDCOxi = new double**[nNodes];
  for (unsigned m=0; m<nNodes; m++)
  {
    v[m] = new double[nElecReactions];
    DvDU[m] = new double[nElecReactions];
    DvDCRed[m] = new double*[nElecReactions];
    DvDCOxi[m] = new double*[nElecReactions];
    for (unsigned r=0; r<nElecReactions; r++)
    {
      unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(r);
      unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(r);
      DvDCRed[m][r] = new double[nAgentsRed];
      DvDCOxi[m][r] = new double[nAgentsOxi];
    }
  }
}
//---------------------------------------------------------------------------
ElecReactionTerm::~ElecReactionTerm()
{
  unsigned nElecReactions = mitrem->getNElecReactions();
  for (unsigned m=0; m<nNodes; m++)
  {
    for (unsigned r=0; r<nElecReactions; r++)
    {
      delete[] DvDCRed[m][r];
      delete[] DvDCOxi[m][r];
    }
    delete[] v[m];
    delete[] DvDU[m];
    delete[] DvDCRed[m];
    delete[] DvDCOxi[m];
  }
  delete[] v;
  delete[] DvDU;
  delete[] DvDCRed;
  delete[] DvDCOxi;
}
//---------------------------------------------------------------------------
void ElecReactionTerm::calcElecReactionCurrentDensities(EmptyEmptyDoubleListList boundaryElementElecReactionCurrentDensity, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential)
{
  double surfaceGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    surfaceGasFraction += surfaceGasFractions[m];
  }
  surfaceGasFraction /= nNodes;

  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    for (unsigned m=0; m<nNodes; m++)
    {
      // Calculate coefficients
      mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
      v[m][s] = mitrem->calcElecReactionCurrentDensity(s,electrodePotential)*(1.-surfaceGasFraction);

      boundaryElementElecReactionCurrentDensity[m][s] = v[m][s];
    }
  }
}
//---------------------------------------------------------------------------
double ElecReactionTerm::calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential)
{
  double surfaceGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    surfaceGasFraction += surfaceGasFractions[m];
  }
  surfaceGasFraction /= nNodes;

  double current = 0.;
  elementSize = boundaryElementProps->calcSize(coordinates);

  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    for (unsigned m=0; m<nNodes; m++)
    {
      // Calculate coefficients
      mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
      v[m][s] = mitrem->calcElecReactionCurrentDensity(s,electrodePotential)*(1.-surfaceGasFraction);

      // Add to current
      current += v[m][s];
    }
  }

  // Average current density
  current /= nNodes;
  // Total current
  current *= elementSize;

  return current;
}
//---------------------------------------------------------------------------

