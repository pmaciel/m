//---------------------------------------------------------------------------

#include "ElecReactionTerm_1D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReactionTerm_1D_Galerkin::ElecReactionTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_)
  : ElecReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_)
{
}
//---------------------------------------------------------------------------
ElecReactionTerm_1D_Galerkin::~ElecReactionTerm_1D_Galerkin()
{
}
//---------------------------------------------------------------------------
void ElecReactionTerm_1D_Galerkin::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
  // Calculate coefficients
  mitrem->init(concentrations[0],potentials[0],temperatures[0],densities[0]);
  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    v[0][s] = mitrem->calcElecReactionRate(s,electrodePotential)*(1.-surfaceGasFractions[0]);
  }

  // Add to boundary element vector
  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
    unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);

    // v(s) gives a contribution to all its agentsRed
    for (unsigned j=0; j<nAgentsRed; j++)
    {
      unsigned i = mitrem->getElecReactionAgentsRed(s,j);
      int stoichRed = mitrem->getElecReactionStoichRed(s,j);

      boundaryElementVec[eq(0,i)] -= v[0][s]*stoichRed;
    }

    // v(s) gives a contribution to all its agentsOxi
    for (unsigned j=0; j<nAgentsOxi; j++)
    {
      unsigned i = mitrem->getElecReactionAgentsOxi(s,j);
      int stoichOxi = mitrem->getElecReactionStoichOxi(s,j);

      boundaryElementVec[eq(0,i)] -= v[0][s]*stoichOxi;
    }

  }
}
//---------------------------------------------------------------------------
void ElecReactionTerm_1D_Galerkin::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
  // Calculate coefficients
  mitrem->init(concentrations[0],potentials[0],temperatures[0],densities[0]);
  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
    unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);
    DvDU[0][r] = mitrem->calcElecReactionRateDerivativeU(s,electrodePotential)*(1.-surfaceGasFractions[0]);
    for (unsigned j=0; j<nAgentsRed; j++)
    {
      DvDCRed[0][r][j] = mitrem->calcElecReactionRateDerivativeCRed(s,electrodePotential,j)*(1.-surfaceGasFractions[0]);
    }
    for (unsigned j=0; j<nAgentsOxi; j++)
    {
      DvDCOxi[0][r][j] = mitrem->calcElecReactionRateDerivativeCOxi(s,electrodePotential,j)*(1.-surfaceGasFractions[0]);
    }
  }

  // Add to boundary element jacobian
  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
    unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);

    // v(s) gave a contribution to all its agentsRed
    for (unsigned k=0; k<nAgentsRed; k++)
    {
      unsigned i = mitrem->getElecReactionAgentsRed(s,k);
      unsigned eq0i = eq(0,i);
      int stoichRed = mitrem->getElecReactionStoichRed(s,k);

      // derivatives of v(s) by its agentsRed
      for (unsigned l=0; l<nAgentsRed; l++)
      {
        unsigned j = mitrem->getElecReactionAgentsRed(s,l);

        boundaryElementJac[eq0i][var(0,j)] += DvDCRed[0][r][l]*stoichRed;
      }

      // derivatives of v(s) by its agentsOxi
      for (unsigned l=0; l<nAgentsOxi; l++)
      {
        unsigned j = mitrem->getElecReactionAgentsOxi(s,l);

        boundaryElementJac[eq0i][var(0,j)] += DvDCOxi[0][r][l]*stoichRed;
      }

      // derivative of v(s) by U
      boundaryElementJac[eq0i][var(0,nIons)] += DvDU[0][r]*stoichRed;
    }

    // v(s) gave a contribution to all its agentsOxi
    for (unsigned k=0; k<nAgentsOxi; k++)
    {
      unsigned i = mitrem->getElecReactionAgentsOxi(s,k);
      unsigned eq0i = eq(0,i);
      int stoichOxi = mitrem->getElecReactionStoichOxi(s,k);

      // derivatives of v(s) by its agentsRed
      for (unsigned l=0; l<nAgentsRed; l++)
      {
        unsigned j = mitrem->getElecReactionAgentsRed(s,l);

        boundaryElementJac[eq0i][var(0,j)] += DvDCRed[0][r][l]*stoichOxi;
      }

      // derivatives of v(s) by its agentsOxi
      for (unsigned l=0; l<nAgentsOxi; l++)
      {
        unsigned j = mitrem->getElecReactionAgentsOxi(s,l);

        boundaryElementJac[eq0i][var(0,j)] += DvDCOxi[0][r][l]*stoichOxi;
      }

      // derivative of v(s) by U
      boundaryElementJac[eq0i][var(0,nIons)] += DvDU[0][r]*stoichOxi;
    }
  }
}
//---------------------------------------------------------------------------

