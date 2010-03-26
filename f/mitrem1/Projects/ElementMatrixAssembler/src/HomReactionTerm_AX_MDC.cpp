//---------------------------------------------------------------------------

#include "HomReactionTerm_AX_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm_AX_MDC::HomReactionTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : HomReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
HomReactionTerm_AX_MDC::~HomReactionTerm_AX_MDC()
{
}
//---------------------------------------------------------------------------
void HomReactionTerm_AX_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;

  // Calculate coefficients
  double R1[3][3][3];
  double R2[3][3][3][3];
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned r=0; r<nHomReactions; r++)
    {
      kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
      kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
    }
    unsigned n = (m+1)%3;
    unsigned p = (m+2)%3;
    R1[m][m][m] = 1150.*coordinates[m][0] + 275.*coordinates[n][0] + 275.*coordinates[p][0];
    R1[m][n][n] =  123.*coordinates[m][0] +  73.*coordinates[n][0] +  34.*coordinates[p][0];
    R1[m][p][p] =  123.*coordinates[m][0] +  34.*coordinates[n][0] +  73.*coordinates[p][0];
    R1[m][m][n] =  275.*coordinates[m][0] + 123.*coordinates[n][0] +  72.*coordinates[p][0];
    R1[m][n][m] = R1[m][m][n];
    R1[m][m][p] =  275.*coordinates[m][0] +  72.*coordinates[n][0] + 123.*coordinates[p][0];
    R1[m][p][m] = R1[m][m][p];
    R1[m][n][p] =   72.*coordinates[m][0] +  34.*coordinates[n][0] +  34.*coordinates[p][0];
    R1[m][p][n] = R1[m][n][p];
    R2[m][m][m][m] = 16660.*coordinates[m][0] + 3748.*coordinates[n][0] + 3748.*coordinates[p][0];
    R2[m][n][n][n] =   757.*coordinates[m][0] +  526.*coordinates[n][0] +  211.*coordinates[p][0];
    R2[m][p][p][p] =   757.*coordinates[m][0] +  211.*coordinates[n][0] +  526.*coordinates[p][0];
    R2[m][m][m][n] =  3748.*coordinates[m][0] + 1516.*coordinates[n][0] +  862.*coordinates[p][0];
    R2[m][m][n][m] = R2[m][m][m][n];
    R2[m][n][m][m] = R2[m][m][m][n];
    R2[m][m][m][p] =  3748.*coordinates[m][0] +  862.*coordinates[n][0] + 1516.*coordinates[p][0];
    R2[m][m][p][m] = R2[m][m][m][p];
    R2[m][p][m][m] = R2[m][m][m][p];
    R2[m][m][n][n] =  1516.*coordinates[m][0] +  757.*coordinates[n][0] +  349.*coordinates[p][0];
    R2[m][n][m][n] = R2[m][m][n][n];
    R2[m][n][n][m] = R2[m][m][n][n];
    R2[m][m][p][p] =  1516.*coordinates[m][0] +  349.*coordinates[n][0] +  757.*coordinates[p][0];
    R2[m][p][m][p] = R2[m][m][p][p];
    R2[m][p][p][m] = R2[m][m][p][p];
    R2[m][m][n][p] =   862.*coordinates[m][0] +  349.*coordinates[n][0] +  349.*coordinates[p][0];
    R2[m][m][p][n] = R2[m][m][n][p];
    R2[m][n][m][p] = R2[m][m][n][p];
    R2[m][p][m][n] = R2[m][m][n][p];
    R2[m][n][p][m] = R2[m][m][n][p];
    R2[m][p][n][m] = R2[m][m][n][p];
    R2[m][n][n][p] =   349.*coordinates[m][0] +  211.*coordinates[n][0] +  160.*coordinates[p][0];
    R2[m][n][p][n] = R2[m][n][n][p];
    R2[m][p][n][n] = R2[m][n][n][p];
    R2[m][p][p][n] =   349.*coordinates[m][0] +  160.*coordinates[n][0] +  211.*coordinates[p][0];
    R2[m][p][n][p] = R2[m][p][p][n];
    R2[m][n][p][p] = R2[m][p][p][n];
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize12960 = elementSize/12960.;
  double elementSize77760 = elementSize/77760.;

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
      unsigned n = (m+1)%3;
      unsigned p = (m+2)%3;

      // Forward reaction
      if (nReagents == 1)
      {
        Hfmj[m] = (R1[m][m][m]*kf[m][r] + R1[m][m][n]*kf[n][r] + R1[m][m][p]*kf[p][r])*elementSize12960;
        Hfmj[n] = (R1[m][n][m]*kf[m][r] + R1[m][n][n]*kf[n][r] + R1[m][n][p]*kf[p][r])*elementSize12960;
        Hfmj[p] = (R1[m][p][m]*kf[m][r] + R1[m][p][n]*kf[n][r] + R1[m][p][p]*kf[p][r])*elementSize12960;

        elementMat[eq(m,reagents[0])][var(m,reagents[0])] -= Hfmj[m];
        elementMat[eq(m,reagents[0])][var(n,reagents[0])] -= Hfmj[n];
        elementMat[eq(m,reagents[0])][var(p,reagents[0])] -= Hfmj[p];
        for (unsigned j=0; j<nProducts; j++)
        {
          elementMat[eq(m,products[j])][var(m,reagents[0])] += Hfmj[m];
          elementMat[eq(m,products[j])][var(n,reagents[0])] += Hfmj[n];
          elementMat[eq(m,products[j])][var(p,reagents[0])] += Hfmj[p];
        }
      }
      else if (nReagents == 2)
      {
        Hfmjk[m][m] = (R2[m][m][m][m]*kf[m][r] + R2[m][m][m][n]*kf[n][r] + R2[m][m][m][p]*kf[p][r])*elementSize77760;
        Hfmjk[n][n] = (R2[m][n][n][m]*kf[m][r] + R2[m][n][n][n]*kf[n][r] + R2[m][n][n][p]*kf[p][r])*elementSize77760;
        Hfmjk[p][p] = (R2[m][p][p][m]*kf[m][r] + R2[m][p][p][n]*kf[n][r] + R2[m][p][p][p]*kf[p][r])*elementSize77760;
        Hfmjk[m][n] = (R2[m][m][n][m]*kf[m][r] + R2[m][m][n][n]*kf[n][r] + R2[m][m][n][p]*kf[p][r])*elementSize77760;
        Hfmjk[n][m] = Hfmjk[m][n];
        Hfmjk[m][p] = (R2[m][m][p][m]*kf[m][r] + R2[m][m][p][n]*kf[n][r] + R2[m][m][p][p]*kf[p][r])*elementSize77760;
        Hfmjk[p][m] = Hfmjk[m][p];
        Hfmjk[n][p] = (R2[m][n][p][m]*kf[m][r] + R2[m][n][p][n]*kf[n][r] + R2[m][n][p][p]*kf[p][r])*elementSize77760;
        Hfmjk[p][n] = Hfmjk[n][p];

        for (unsigned j=0; j<nReagents; j++)
        {
          unsigned k = (j+1)%2;
          Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]] + Hfmjk[m][p]*concentrations[p][reagents[k]]);
          Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]] + Hfmjk[n][p]*concentrations[p][reagents[k]]);
          Hfmj[p] = 0.5*(Hfmjk[p][m]*concentrations[m][reagents[k]] + Hfmjk[p][n]*concentrations[n][reagents[k]] + Hfmjk[p][p]*concentrations[p][reagents[k]]);

          elementMat[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
          elementMat[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
          elementMat[eq(m,reagents[0])][var(p,reagents[j])] -= Hfmj[p];
          elementMat[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
          elementMat[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
          elementMat[eq(m,reagents[1])][var(p,reagents[j])] -= Hfmj[p];
          for (unsigned h = 0; h < nProducts; h++)
          {
            elementMat[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
            elementMat[eq(m,products[h])][var(n,reagents[j])] += Hfmj[n];
            elementMat[eq(m,products[h])][var(p,reagents[j])] += Hfmj[p];
          }
        }
      }

      // Backward reaction
      if (nProducts == 1)
      {
        Hbmj[m] = (R1[m][m][m]*kb[m][r] + R1[m][m][n]*kb[n][r] + R1[m][m][p]*kb[p][r])*elementSize12960;
        Hbmj[n] = (R1[m][n][m]*kb[m][r] + R1[m][n][n]*kb[n][r] + R1[m][n][p]*kb[p][r])*elementSize12960;
        Hbmj[p] = (R1[m][p][m]*kb[m][r] + R1[m][p][n]*kb[n][r] + R1[m][p][p]*kb[p][r])*elementSize12960;

        elementMat[eq(m,products[0])][var(m,products[0])] -= Hbmj[m];
        elementMat[eq(m,products[0])][var(n,products[0])] -= Hbmj[n];
        elementMat[eq(m,products[0])][var(p,products[0])] -= Hbmj[p];
        for (unsigned j=0; j<nReagents; j++)
        {
          elementMat[eq(m,reagents[j])][var(m,products[0])] += Hbmj[m];
          elementMat[eq(m,reagents[j])][var(n,products[0])] += Hbmj[n];
          elementMat[eq(m,reagents[j])][var(p,products[0])] += Hbmj[p];
        }
      }
      else if (nProducts == 2)
      {
        Hbmjk[m][m] = (R2[m][m][m][m]*kb[m][r] + R2[m][m][m][n]*kb[n][r] + R2[m][m][m][p]*kb[p][r])*elementSize77760;
        Hbmjk[n][n] = (R2[m][n][n][m]*kb[m][r] + R2[m][n][n][n]*kb[n][r] + R2[m][n][n][p]*kb[p][r])*elementSize77760;
        Hbmjk[p][p] = (R2[m][p][p][m]*kb[m][r] + R2[m][p][p][n]*kb[n][r] + R2[m][p][p][p]*kb[p][r])*elementSize77760;
        Hbmjk[m][n] = (R2[m][m][n][m]*kb[m][r] + R2[m][m][n][n]*kb[n][r] + R2[m][m][n][p]*kb[p][r])*elementSize77760;
        Hbmjk[n][m] = Hbmjk[m][n];
        Hbmjk[m][p] = (R2[m][m][p][m]*kb[m][r] + R2[m][m][p][n]*kb[n][r] + R2[m][m][p][p]*kb[p][r])*elementSize77760;
        Hbmjk[p][m] = Hbmjk[m][p];
        Hbmjk[n][p] = (R2[m][n][p][m]*kb[m][r] + R2[m][n][p][n]*kb[n][r] + R2[m][n][p][p]*kb[p][r])*elementSize77760;
        Hbmjk[p][n] = Hbmjk[n][p];

        for (unsigned j=0; j<nProducts; j++)
        {
          unsigned k = (j+1)%2;
          Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]] + Hbmjk[m][p]*concentrations[p][products[k]]);
          Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]] + Hbmjk[n][p]*concentrations[p][products[k]]);
          Hbmj[p] = 0.5*(Hbmjk[p][m]*concentrations[m][products[k]] + Hbmjk[p][n]*concentrations[n][products[k]] + Hbmjk[p][p]*concentrations[p][products[k]]);

          elementMat[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
          elementMat[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
          elementMat[eq(m,products[0])][var(p,products[j])] -= Hbmj[p];
          elementMat[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
          elementMat[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
          elementMat[eq(m,products[1])][var(p,products[j])] -= Hbmj[p];
          for (unsigned h = 0; h < nReagents; h++)
          {
            elementMat[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
            elementMat[eq(m,reagents[h])][var(n,products[j])] += Hbmj[n];
            elementMat[eq(m,reagents[h])][var(p,products[j])] += Hbmj[p];
          }
        }
      }
    }
    delete[] products;
    delete[] reagents;
  }
}
//---------------------------------------------------------------------------
void HomReactionTerm_AX_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  double volumeGasFraction = 0.;
  for (unsigned m=0; m<nNodes; m++)
  {
    volumeGasFraction += volumeGasFractions[m];
  }
  volumeGasFraction /= 3.;

  // Calculate coefficients
  double R2[3][3][3][3];
  for (unsigned m=0; m<nNodes; m++)
  {
    mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
    for (unsigned r=0; r<nHomReactions; r++)
    {
      kf[m][r] = mitrem->calcHomReactionForwardRateConstant(r)*(1.-volumeGasFraction);
      kb[m][r] = mitrem->calcHomReactionBackwardRateConstant(r)*(1.-volumeGasFraction);
    }
    unsigned n = (m+1)%3;
    unsigned p = (m+2)%3;
    R2[m][m][m][m] = 16660.*coordinates[m][0] + 3748.*coordinates[n][0] + 3748.*coordinates[p][0];
    R2[m][n][n][n] =   757.*coordinates[m][0] +  526.*coordinates[n][0] +  211.*coordinates[p][0];
    R2[m][p][p][p] =   757.*coordinates[m][0] +  211.*coordinates[n][0] +  526.*coordinates[p][0];
    R2[m][m][m][n] =  3748.*coordinates[m][0] + 1516.*coordinates[n][0] +  862.*coordinates[p][0];
    R2[m][m][n][m] = R2[m][m][m][n];
    R2[m][n][m][m] = R2[m][m][m][n];
    R2[m][m][m][p] =  3748.*coordinates[m][0] +  862.*coordinates[n][0] + 1516.*coordinates[p][0];
    R2[m][m][p][m] = R2[m][m][m][p];
    R2[m][p][m][m] = R2[m][m][m][p];
    R2[m][m][n][n] =  1516.*coordinates[m][0] +  757.*coordinates[n][0] +  349.*coordinates[p][0];
    R2[m][n][m][n] = R2[m][m][n][n];
    R2[m][n][n][m] = R2[m][m][n][n];
    R2[m][m][p][p] =  1516.*coordinates[m][0] +  349.*coordinates[n][0] +  757.*coordinates[p][0];
    R2[m][p][m][p] = R2[m][m][p][p];
    R2[m][p][p][m] = R2[m][m][p][p];
    R2[m][m][n][p] =   862.*coordinates[m][0] +  349.*coordinates[n][0] +  349.*coordinates[p][0];
    R2[m][m][p][n] = R2[m][m][n][p];
    R2[m][n][m][p] = R2[m][m][n][p];
    R2[m][p][m][n] = R2[m][m][n][p];
    R2[m][n][p][m] = R2[m][m][n][p];
    R2[m][p][n][m] = R2[m][m][n][p];
    R2[m][n][n][p] =   349.*coordinates[m][0] +  211.*coordinates[n][0] +  160.*coordinates[p][0];
    R2[m][n][p][n] = R2[m][n][n][p];
    R2[m][p][n][n] = R2[m][n][n][p];
    R2[m][p][p][n] =   349.*coordinates[m][0] +  160.*coordinates[n][0] +  211.*coordinates[p][0];
    R2[m][p][n][p] = R2[m][p][p][n];
    R2[m][n][p][p] = R2[m][p][p][n];
  }
  elementSize = elementProps->calcSize(coordinates);
  double elementSize77760 = elementSize/77760.;

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
      unsigned n = (m+1)%3;
      unsigned p = (m+2)%3;

      // Forward reaction
      if (nReagents == 2)
      {
        Hfmjk[m][m] = (R2[m][m][m][m]*kf[m][r] + R2[m][m][m][n]*kf[n][r] + R2[m][m][m][p]*kf[p][r])*elementSize77760;
        Hfmjk[n][n] = (R2[m][n][n][m]*kf[m][r] + R2[m][n][n][n]*kf[n][r] + R2[m][n][n][p]*kf[p][r])*elementSize77760;
        Hfmjk[p][p] = (R2[m][p][p][m]*kf[m][r] + R2[m][p][p][n]*kf[n][r] + R2[m][p][p][p]*kf[p][r])*elementSize77760;
        Hfmjk[m][n] = (R2[m][m][n][m]*kf[m][r] + R2[m][m][n][n]*kf[n][r] + R2[m][m][n][p]*kf[p][r])*elementSize77760;
        Hfmjk[n][m] = Hfmjk[m][n];
        Hfmjk[m][p] = (R2[m][m][p][m]*kf[m][r] + R2[m][m][p][n]*kf[n][r] + R2[m][m][p][p]*kf[p][r])*elementSize77760;
        Hfmjk[p][m] = Hfmjk[m][p];
        Hfmjk[n][p] = (R2[m][n][p][m]*kf[m][r] + R2[m][n][p][n]*kf[n][r] + R2[m][n][p][p]*kf[p][r])*elementSize77760;
        Hfmjk[p][n] = Hfmjk[n][p];

        for (unsigned j=0; j<nReagents; j++)
        {
          unsigned k = (j+1)%2;
          Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]] + Hfmjk[m][p]*concentrations[p][reagents[k]]);
          Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]] + Hfmjk[n][p]*concentrations[p][reagents[k]]);
          Hfmj[p] = 0.5*(Hfmjk[p][m]*concentrations[m][reagents[k]] + Hfmjk[p][n]*concentrations[n][reagents[k]] + Hfmjk[p][p]*concentrations[p][reagents[k]]);

          elementJac[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
          elementJac[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
          elementJac[eq(m,reagents[0])][var(p,reagents[j])] -= Hfmj[p];
          elementJac[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
          elementJac[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
          elementJac[eq(m,reagents[1])][var(p,reagents[j])] -= Hfmj[p];
          for (unsigned h = 0; h < nProducts; h++)
          {
            elementJac[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
            elementJac[eq(m,products[h])][var(n,reagents[j])] += Hfmj[n];
            elementJac[eq(m,products[h])][var(p,reagents[j])] += Hfmj[p];
          }
        }
      }

      // Backward reaction
      if (nProducts == 2)
      {
        Hbmjk[m][m] = (R2[m][m][m][m]*kb[m][r] + R2[m][m][m][n]*kb[n][r] + R2[m][m][m][p]*kb[p][r])*elementSize77760;
        Hbmjk[n][n] = (R2[m][n][n][m]*kb[m][r] + R2[m][n][n][n]*kb[n][r] + R2[m][n][n][p]*kb[p][r])*elementSize77760;
        Hbmjk[p][p] = (R2[m][p][p][m]*kb[m][r] + R2[m][p][p][n]*kb[n][r] + R2[m][p][p][p]*kb[p][r])*elementSize77760;
        Hbmjk[m][n] = (R2[m][m][n][m]*kb[m][r] + R2[m][m][n][n]*kb[n][r] + R2[m][m][n][p]*kb[p][r])*elementSize77760;
        Hbmjk[n][m] = Hbmjk[m][n];
        Hbmjk[m][p] = (R2[m][m][p][m]*kb[m][r] + R2[m][m][p][n]*kb[n][r] + R2[m][m][p][p]*kb[p][r])*elementSize77760;
        Hbmjk[p][m] = Hbmjk[m][p];
        Hbmjk[n][p] = (R2[m][n][p][m]*kb[m][r] + R2[m][n][p][n]*kb[n][r] + R2[m][n][p][p]*kb[p][r])*elementSize77760;
        Hbmjk[p][n] = Hbmjk[n][p];

        for (unsigned j=0; j<nProducts; j++)
        {
          unsigned k = (j+1)%2;
          Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]] + Hbmjk[m][p]*concentrations[p][products[k]]);
          Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]] + Hbmjk[n][p]*concentrations[p][products[k]]);
          Hbmj[p] = 0.5*(Hbmjk[p][m]*concentrations[m][products[k]] + Hbmjk[p][n]*concentrations[n][products[k]] + Hbmjk[p][p]*concentrations[p][products[k]]);

          elementJac[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
          elementJac[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
          elementJac[eq(m,products[0])][var(p,products[j])] -= Hbmj[p];
          elementJac[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
          elementJac[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
          elementJac[eq(m,products[1])][var(p,products[j])] -= Hbmj[p];
          for (unsigned h = 0; h < nReagents; h++)
          {
            elementJac[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
            elementJac[eq(m,reagents[h])][var(n,products[j])] += Hbmj[n];
            elementJac[eq(m,reagents[h])][var(p,products[j])] += Hbmj[p];
          }
        }
      }
    }
    delete[] products;
    delete[] reagents;
  }
}
//---------------------------------------------------------------------------

