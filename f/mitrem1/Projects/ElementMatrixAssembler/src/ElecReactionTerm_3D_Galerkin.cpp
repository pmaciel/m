//---------------------------------------------------------------------------

#include "ElecReactionTerm_3D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReactionTerm_3D_Galerkin::ElecReactionTerm_3D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_)
  : ElecReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, boundaryElementProps_)
{
}
//---------------------------------------------------------------------------
ElecReactionTerm_3D_Galerkin::~ElecReactionTerm_3D_Galerkin()
{
}
//---------------------------------------------------------------------------
void ElecReactionTerm_3D_Galerkin::calcVec(double* boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
  double surfaceGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    surfaceGasFraction += surfaceGasFractions[m];
  }
  surfaceGasFraction /= 3.;

  // Calculate coefficients
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned r=0; r<nElecReactions; r++)
    {
      unsigned s = elecReactions[r];
      v[m][s] = mitrem->calcElecReactionRate(s,electrodePotential)*(1.-surfaceGasFraction);
    }
  }
  elementSize = boundaryElementProps->calcSize(coordinates);
  double elementSize12 = elementSize/12.;

  // Add to boundary element vector
  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
    unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);
    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned p = (m+1)%3;
      unsigned q = (m+2)%3;
      double integral = (2.*v[m][s] + v[p][s] + v[q][s])*elementSize12;

      // v(r) gives a contribution to all its RedAgents
      for (unsigned j=0; j<nAgentsRed; j++)
      {
        unsigned i = mitrem->getElecReactionAgentsRed(s,j);
        int stoichRed = mitrem->getElecReactionStoichRed(s,j);

        boundaryElementVec[eq(m,i)] -= integral*stoichRed;
      }

      // v(r) gives a contribution to all its OxiAgents
      for (unsigned j=0; j<nAgentsOxi; j++)
      {
        unsigned i = mitrem->getElecReactionAgentsOxi(s,j);
        int stoichOxi = mitrem->getElecReactionStoichOxi(s,j);

        boundaryElementVec[eq(m,i)] -= integral*stoichOxi;
      }

    }
  }
}
//---------------------------------------------------------------------------
void ElecReactionTerm_3D_Galerkin::calcJac(double** boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
  double surfaceGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    surfaceGasFraction += surfaceGasFractions[m];
  }
  surfaceGasFraction /= 3.;

  // Calculate coefficients
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned r=0; r<nElecReactions; r++)
    {
      unsigned s = elecReactions[r];
      unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
      unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);

      DvDU[m][r] = mitrem->calcElecReactionRateDerivativeU(s,electrodePotential)*(1.-surfaceGasFraction);
      for (unsigned j=0; j<nAgentsRed; j++)
      {
        DvDCRed[m][r][j] = mitrem->calcElecReactionRateDerivativeCRed(s,electrodePotential,j)*(1.-surfaceGasFraction);
      }
      for (unsigned j=0; j<nAgentsOxi; j++)
      {
        DvDCOxi[m][r][j] = mitrem->calcElecReactionRateDerivativeCOxi(s,electrodePotential,j)*(1.-surfaceGasFraction);
      }
    }
  }
  elementSize = boundaryElementProps->calcSize(coordinates);
  double elementSize12 = elementSize/12.;

  // Add to boundary element jacobian
  for (unsigned r=0; r<nElecReactions; r++)
  {
    unsigned s = elecReactions[r];
    unsigned nAgentsRed = mitrem->getElecReactionNAgentsRed(s);
    unsigned nAgentsOxi = mitrem->getElecReactionNAgentsOxi(s);

    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned p = (m+1)%3;
      unsigned q = (m+2)%3;

      // v(s) gave a contribution to all its agentsRed
      for (unsigned k=0; k<nAgentsRed; k++)
      {
        unsigned i = mitrem->getElecReactionAgentsRed(s,k);
        unsigned eqmi = eq(m,i);
        int stoichRed = mitrem->getElecReactionStoichRed(s,k);

        // derivatives of v(s) by its agentsRed
        for (unsigned l=0; l<nAgentsRed; l++)
        {
          unsigned j = mitrem->getElecReactionAgentsRed(s,l);
          boundaryElementJac[eqmi][var(m,j)] += 2.*DvDCRed[m][r][l]*elementSize12*stoichRed;
          boundaryElementJac[eqmi][var(p,j)] +=    DvDCRed[p][r][l]*elementSize12*stoichRed;
          boundaryElementJac[eqmi][var(q,j)] +=    DvDCRed[q][r][l]*elementSize12*stoichRed;
        }

        // derivatives of v(s) by its agentsOxi
        for (unsigned l=0; l<nAgentsOxi; l++)
        {
          unsigned j = mitrem->getElecReactionAgentsOxi(s,l);
          boundaryElementJac[eqmi][var(m,j)] += 2.*DvDCOxi[m][r][l]*elementSize12*stoichRed;
          boundaryElementJac[eqmi][var(p,j)] +=    DvDCOxi[p][r][l]*elementSize12*stoichRed;
          boundaryElementJac[eqmi][var(q,j)] +=    DvDCOxi[q][r][l]*elementSize12*stoichRed;
        }

        // derivative of v(s) by U
        boundaryElementJac[eqmi][var(m,nIons)] += 2.*DvDU[m][r]*elementSize12*stoichRed;
        boundaryElementJac[eqmi][var(p,nIons)] +=    DvDU[p][r]*elementSize12*stoichRed;
        boundaryElementJac[eqmi][var(q,nIons)] +=    DvDU[q][r]*elementSize12*stoichRed;
      }

      // v(s) gave a contribution to all its agentsOxi
      for (unsigned k=0; k<nAgentsOxi; k++)
      {
        unsigned i = mitrem->getElecReactionAgentsOxi(s,k);
        unsigned eqmi = eq(m,i);
        int stoichOxi = mitrem->getElecReactionStoichOxi(s,k);

        // derivatives of v(s) by its agentsRed
        for (unsigned l=0; l<nAgentsRed; l++)
        {
          unsigned j = mitrem->getElecReactionAgentsRed(s,l);
          boundaryElementJac[eqmi][var(m,j)] += 2.*DvDCRed[m][r][l]*elementSize12*stoichOxi;
          boundaryElementJac[eqmi][var(p,j)] +=    DvDCRed[p][r][l]*elementSize12*stoichOxi;
          boundaryElementJac[eqmi][var(q,j)] +=    DvDCRed[q][r][l]*elementSize12*stoichOxi;
        }

        // derivatives of v(s) by its agentsOxi
        for (unsigned l=0; l<nAgentsOxi; l++)
        {
          unsigned j = mitrem->getElecReactionAgentsOxi(s,l);
          boundaryElementJac[eqmi][var(m,j)] += 2.*DvDCOxi[m][r][l]*elementSize12*stoichOxi;
          boundaryElementJac[eqmi][var(p,j)] +=    DvDCOxi[p][r][l]*elementSize12*stoichOxi;
          boundaryElementJac[eqmi][var(q,j)] +=    DvDCOxi[q][r][l]*elementSize12*stoichOxi;
        }

        // derivative of v(s) by U
        boundaryElementJac[eqmi][var(m,nIons)] += 2.*DvDU[m][r]*elementSize12*stoichOxi;
        boundaryElementJac[eqmi][var(p,nIons)] +=    DvDU[p][r]*elementSize12*stoichOxi;
        boundaryElementJac[eqmi][var(q,nIons)] +=    DvDU[q][r]*elementSize12*stoichOxi;
      }
    }
  }
}
//---------------------------------------------------------------------------

