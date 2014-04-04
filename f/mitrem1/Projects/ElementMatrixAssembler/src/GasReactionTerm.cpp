//---------------------------------------------------------------------------

#include "GasReactionTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
GasReactionTerm::GasReactionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_)
  : BoundaryElementTerm(nDimensions_, nNodes_, nVariables_, mitrem_, boundaryElementProps_),
    isBubble(isBubble_)
{
  unsigned nGasReactions = mitrem->getNGasReactions();
  v = new double*[nNodes];
  DvDCDissGas = new double*[nNodes];
  for (unsigned m=0; m<nNodes; m++)
  {
    v[m] = new double[nGasReactions];
    DvDCDissGas[m] = new double[nGasReactions];
  }
}
//---------------------------------------------------------------------------
GasReactionTerm::~GasReactionTerm()
{
  //unsigned nGasReactions = mitrem->getNGasReactions();
  for (unsigned m=0; m<nNodes; m++)
  {
    delete[] v[m];
    delete[] DvDCDissGas[m];
  }
  delete[] v;
  delete[] DvDCDissGas;
}
//---------------------------------------------------------------------------
void GasReactionTerm::calcGasReactionRates(EmptyEmptyDoubleListList boundaryElementGasReactionRate, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions)
{
  for (unsigned r=0; r<nGasReactions; r++)
  {
    unsigned s = gasReactions[r];
    for (unsigned m=0; m<nNodes; m++)
    {
      // Calculate coefficients
      mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
      v[m][s] = mitrem->calcGasReactionRate(s) * calcBubbleReactionRateCorrection(surfaceGasFractions);

      boundaryElementGasReactionRate[m][s] = v[m][s];
    }
  }
}
//---------------------------------------------------------------------------
double GasReactionTerm::calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions)
{
  double gasGeneration = 0.;
  elementSize = boundaryElementProps->calcSize(coordinates);

  for (unsigned r=0; r<nGasReactions; r++)
  {
    unsigned s = gasReactions[r];
    for (unsigned m=0; m<nNodes; m++)
    {
      // Calculate coefficients
      mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
      v[m][s] = mitrem->calcGasReactionRate(s) * calcBubbleReactionRateCorrection(surfaceGasFractions);

      // Add to gas generation
      gasGeneration += v[m][s];
    }
  }

  // Average gas generation density
  gasGeneration /= nNodes;
  // Total gas generation
  gasGeneration *= elementSize*R_CONST*mitrem->getSolutionTemperature()/101300.;

  return gasGeneration;
}
//---------------------------------------------------------------------------
double GasReactionTerm::calcBubbleReactionRateCorrection(DoubleList surfaceGasFractions)
{
  // Heidi and Steven think it is not necessary to correct the gas reaction rates for the surface coverage
  /*
    double surfaceGasFraction = 0.;
    for (unsigned m=0; m<nNodes; m++)
      surfaceGasFraction += surfaceGasFractions[m];
    surfaceGasFraction /= nNodes;
  */
  return 1.;
}
//---------------------------------------------------------------------------
