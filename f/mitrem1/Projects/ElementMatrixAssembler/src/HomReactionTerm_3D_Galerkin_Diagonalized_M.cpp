//---------------------------------------------------------------------------

#include "HomReactionTerm_3D_Galerkin_Diagonalized.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm_3D_Galerkin_Diagonalized::HomReactionTerm_3D_Galerkin_Diagonalized(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : HomReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
HomReactionTerm_3D_Galerkin_Diagonalized::~HomReactionTerm_3D_Galerkin_Diagonalized()
{
}
//---------------------------------------------------------------------------
void HomReactionTerm_3D_Galerkin_Diagonalized::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.25;

  // Calculate coefficients
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned r=0; r<nHomReactions; r++)
    {
      kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
      kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize20 = elementSize/20.;

  // Add to element matrix
  for (unsigned r=0; r<nHomReactions; r++)
  {
    unsigned nReagents = mitrem->getHomReactionNReagents(r);
    unsigned nProducts = mitrem->getHomReactionNProducts(r);
    unsigned* reagents = new unsigned[nReagents];
    unsigned* products = new unsigned[nProducts];
    for (unsigned j=0; j<nReagents; j++)
    {
      reagents[j] = mitrem->getHomReactionReagents(r,j);
    }
    for (unsigned j=0; j<nProducts; j++)
    {
      products[j] = mitrem->getHomReactionProducts(r,j);
    }

    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned n = (m+1)%4;
      unsigned p = (m+2)%4;
      unsigned q = (m+3)%4;

      // Forward reaction
      if (nReagents == 1)
      {
        Hfmj[m] = (2.*kf[m][r] +    kf[n][r] +    kf[p][r] +    kf[q][r])*elementSize20;

        elementMat[eq(m,reagents[0])][var(m,reagents[0])] -= Hfmj[m];
        for (unsigned j=0; j<nProducts; j++)
        {
          elementMat[eq(m,products[j])][var(m,reagents[0])] += Hfmj[m];
        }
      }
      else if (nReagents == 2)
      {
        Hfmjk[m][m] = ( 2.*kf[m][r] +    kf[n][r] +    kf[p][r] +    kf[q][r])*elementSize20;

        for (unsigned j=0; j<nReagents; j++)
        {
          unsigned k = (j+1)%2;
          Hfmj[m] = 0.5*Hfmjk[m][m]*concentrations[m][reagents[k]];

          elementMat[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
          elementMat[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
          for (unsigned h = 0; h < nProducts; h++)
          {
            elementMat[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
          }

        }
      }

      // Backward reaction
      if (nProducts == 1)
      {
        Hbmj[m] = (2.*kb[m][r] +    kb[n][r] +    kb[p][r] +    kb[q][r])*elementSize20;

        elementMat[eq(m,products[0])][var(m,products[0])] -= Hbmj[m];
        for (unsigned j=0; j<nReagents; j++)
        {
          elementMat[eq(m,reagents[j])][var(m,products[0])] += Hbmj[m];
        }
      }
      else if (nProducts == 2)
      {
        Hbmjk[m][m] = ( 2.*kb[m][r] +    kb[n][r] +    kb[p][r] +    kb[q][r])*elementSize20;

        for (unsigned j=0; j<nProducts; j++)
        {
          unsigned k = (j+1)%2;
          Hbmj[m] = 0.5*Hbmjk[m][m]*concentrations[m][products[k]];

          elementMat[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
          elementMat[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
          for (unsigned h = 0; h < nReagents; h++)
          {
            elementMat[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
          }

        }
      }
    }
    delete[] products;
    delete[] reagents;
  }
}
//---------------------------------------------------------------------------
void HomReactionTerm_3D_Galerkin_Diagonalized::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction *= 0.25;

  // Calculate coefficients
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned r=0; r<nHomReactions; r++)
    {
      kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
      kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
    }
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize20 = elementSize/20.;

  // Add to element matrix
  for (unsigned r=0; r<nHomReactions; r++)
  {
    unsigned nReagents = mitrem->getHomReactionNReagents(r);
    unsigned nProducts = mitrem->getHomReactionNProducts(r);
    unsigned* reagents = new unsigned[nReagents];
    unsigned* products = new unsigned[nProducts];
    for (unsigned j=0; j<nReagents; j++)
    {
      reagents[j] = mitrem->getHomReactionReagents(r,j);
    }
    for (unsigned j=0; j<nProducts; j++)
    {
      products[j] = mitrem->getHomReactionProducts(r,j);
    }

    for (unsigned m=0; m<nNodes; m++)
    {
      unsigned n = (m+1)%4;
      unsigned p = (m+2)%4;
      unsigned q = (m+3)%4;

      // Forward reaction
      if (nReagents == 2)
      {
        Hfmjk[m][m] = ( 2.*kf[m][r] +    kf[n][r] +    kf[p][r] +    kf[q][r])*elementSize20;

        for (unsigned j=0; j<nReagents; j++)
        {
          unsigned k = (j+1)%2;
          Hfmj[m] = 0.5*Hfmjk[m][m]*concentrations[m][reagents[k]];

          elementJac[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
          elementJac[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
          for (unsigned h = 0; h < nProducts; h++)
          {
            elementJac[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
          }

        }
      }

      // Backward reaction
      if (nProducts == 2)
      {
        Hbmjk[m][m] = ( 2.*kb[m][r] +    kb[n][r] +    kb[p][r] +    kb[q][r])*elementSize20;

        for (unsigned j=0; j<nProducts; j++)
        {
          unsigned k = (j+1)%2;
          Hbmj[m] = 0.5*Hbmjk[m][m]*concentrations[m][products[k]];

          elementJac[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
          elementJac[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
          for (unsigned h = 0; h < nReagents; h++)
          {
            elementJac[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
          }
        }
      }
    }

    delete[] products;
    delete[] reagents;
  }
}
//---------------------------------------------------------------------------

