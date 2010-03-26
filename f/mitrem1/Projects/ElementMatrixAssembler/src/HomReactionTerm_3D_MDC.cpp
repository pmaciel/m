//---------------------------------------------------------------------------

#include "HomReactionTerm_3D_MDC.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm_3D_MDC::HomReactionTerm_3D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : HomReactionTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
HomReactionTerm_3D_MDC::~HomReactionTerm_3D_MDC()
{
}
//---------------------------------------------------------------------------
void HomReactionTerm_3D_MDC::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
  double elementSize17280 = elementSize/17280.;
  double elementSize414720 = elementSize/414720.;

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
        Hfmj[m] = (1245.*kf[m][r] + 335.*kf[n][r] + 335.*kf[p][r] + 335.*kf[q][r])*elementSize17280;
        Hfmj[n] = ( 335.*kf[m][r] + 161.*kf[n][r] +  97.*kf[p][r] +  97.*kf[q][r])*elementSize17280;
        Hfmj[p] = ( 335.*kf[m][r] +  97.*kf[n][r] + 161.*kf[p][r] +  97.*kf[q][r])*elementSize17280;
        Hfmj[q] = ( 335.*kf[m][r] +  97.*kf[n][r] +  97.*kf[p][r] + 161.*kf[q][r])*elementSize17280;

        elementMat[eq(m,reagents[0])][var(m,reagents[0])] -= Hfmj[m];
        elementMat[eq(m,reagents[0])][var(n,reagents[0])] -= Hfmj[n];
        elementMat[eq(m,reagents[0])][var(p,reagents[0])] -= Hfmj[p];
        elementMat[eq(m,reagents[0])][var(q,reagents[0])] -= Hfmj[q];
        for (unsigned j=0; j<nProducts; j++)
        {
          elementMat[eq(m,products[j])][var(m,reagents[0])] += Hfmj[m];
          elementMat[eq(m,products[j])][var(n,reagents[0])] += Hfmj[n];
          elementMat[eq(m,products[j])][var(p,reagents[0])] += Hfmj[p];
          elementMat[eq(m,products[j])][var(q,reagents[0])] += Hfmj[q];
        }
      }
      else if (nReagents == 2)
      {
        Hfmjk[m][m] = (17535.*kf[m][r] + 4115.*kf[n][r] + 4115.*kf[p][r] + 4115.*kf[q][r])*elementSize414720;
        Hfmjk[n][n] = ( 1815.*kf[m][r] + 1067.*kf[n][r] +  491.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[p][p] = ( 1815.*kf[m][r] +  491.*kf[n][r] + 1067.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[q][q] = ( 1815.*kf[m][r] +  491.*kf[n][r] +  491.*kf[p][r] + 1067.*kf[q][r])*elementSize414720;
        Hfmjk[m][n] = ( 4115.*kf[m][r] + 1815.*kf[n][r] + 1055.*kf[p][r] + 1055.*kf[q][r])*elementSize414720;
        Hfmjk[n][m] = Hfmjk[m][n];
        Hfmjk[m][p] = ( 4115.*kf[m][r] + 1055.*kf[n][r] + 1815.*kf[p][r] + 1055.*kf[q][r])*elementSize414720;
        Hfmjk[p][m] = Hfmjk[m][p];
        Hfmjk[m][q] = ( 4115.*kf[m][r] + 1055.*kf[n][r] + 1055.*kf[p][r] + 1815.*kf[q][r])*elementSize414720;
        Hfmjk[q][m] = Hfmjk[m][p];
        Hfmjk[n][p] = ( 1055.*kf[m][r] +  491.*kf[n][r] +  491.*kf[p][r] +  291.*kf[q][r])*elementSize414720;
        Hfmjk[p][n] = Hfmjk[n][p];
        Hfmjk[n][q] = ( 1055.*kf[m][r] +  491.*kf[n][r] +  291.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[q][n] = Hfmjk[n][q];
        Hfmjk[p][q] = ( 1055.*kf[m][r] +  291.*kf[n][r] +  491.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[q][p] = Hfmjk[p][q];

        for (unsigned j=0; j<nReagents; j++)
        {
          unsigned k = (j+1)%2;
          Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]] + Hfmjk[m][p]*concentrations[p][reagents[k]] + Hfmjk[m][q]*concentrations[q][reagents[k]]);
          Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]] + Hfmjk[n][p]*concentrations[p][reagents[k]] + Hfmjk[n][q]*concentrations[q][reagents[k]]);
          Hfmj[p] = 0.5*(Hfmjk[p][m]*concentrations[m][reagents[k]] + Hfmjk[p][n]*concentrations[n][reagents[k]] + Hfmjk[p][p]*concentrations[p][reagents[k]] + Hfmjk[p][q]*concentrations[q][reagents[k]]);
          Hfmj[q] = 0.5*(Hfmjk[q][m]*concentrations[m][reagents[k]] + Hfmjk[q][n]*concentrations[n][reagents[k]] + Hfmjk[q][p]*concentrations[p][reagents[k]] + Hfmjk[q][q]*concentrations[q][reagents[k]]);

          elementMat[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
          elementMat[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
          elementMat[eq(m,reagents[0])][var(p,reagents[j])] -= Hfmj[p];
          elementMat[eq(m,reagents[0])][var(q,reagents[j])] -= Hfmj[q];
          elementMat[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
          elementMat[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
          elementMat[eq(m,reagents[1])][var(p,reagents[j])] -= Hfmj[p];
          elementMat[eq(m,reagents[1])][var(q,reagents[j])] -= Hfmj[q];
          for (unsigned h = 0; h < nProducts; h++)
          {
            elementMat[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
            elementMat[eq(m,products[h])][var(n,reagents[j])] += Hfmj[n];
            elementMat[eq(m,products[h])][var(p,reagents[j])] += Hfmj[p];
            elementMat[eq(m,products[h])][var(q,reagents[j])] += Hfmj[q];
          }

        }
      }

      // Backward reaction
      if (nProducts == 1)
      {
        Hbmj[m] = (1245.*kb[m][r] + 335.*kb[n][r] + 335.*kb[p][r] + 335.*kb[q][r])*elementSize17280;
        Hbmj[n] = ( 335.*kb[m][r] + 161.*kb[n][r] +  97.*kb[p][r] +  97.*kb[q][r])*elementSize17280;
        Hbmj[p] = ( 335.*kb[m][r] +  97.*kb[n][r] + 161.*kb[p][r] +  97.*kb[q][r])*elementSize17280;
        Hbmj[q] = ( 335.*kb[m][r] +  97.*kb[n][r] +  97.*kb[p][r] + 161.*kb[q][r])*elementSize17280;

        elementMat[eq(m,products[0])][var(m,products[0])] -= Hbmj[m];
        elementMat[eq(m,products[0])][var(n,products[0])] -= Hbmj[n];
        elementMat[eq(m,products[0])][var(p,products[0])] -= Hbmj[p];
        elementMat[eq(m,products[0])][var(q,products[0])] -= Hbmj[q];
        for (unsigned j=0; j<nReagents; j++)
        {
          elementMat[eq(m,reagents[j])][var(m,products[0])] += Hbmj[m];
          elementMat[eq(m,reagents[j])][var(n,products[0])] += Hbmj[n];
          elementMat[eq(m,reagents[j])][var(p,products[0])] += Hbmj[p];
          elementMat[eq(m,reagents[j])][var(q,products[0])] += Hbmj[q];
        }
      }
      else if (nProducts == 2)
      {
        Hbmjk[m][m] = (17535.*kb[m][r] + 4115.*kb[n][r] + 4115.*kb[p][r] + 4115.*kb[q][r])*elementSize414720;
        Hbmjk[n][n] = ( 1815.*kb[m][r] + 1067.*kb[n][r] +  491.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[p][p] = ( 1815.*kb[m][r] +  491.*kb[n][r] + 1067.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[q][q] = ( 1815.*kb[m][r] +  491.*kb[n][r] +  491.*kb[p][r] + 1067.*kb[q][r])*elementSize414720;
        Hbmjk[m][n] = ( 4115.*kb[m][r] + 1815.*kb[n][r] + 1055.*kb[p][r] + 1055.*kb[q][r])*elementSize414720;
        Hbmjk[n][m] = Hbmjk[m][n];
        Hbmjk[m][p] = ( 4115.*kb[m][r] + 1055.*kb[n][r] + 1815.*kb[p][r] + 1055.*kb[q][r])*elementSize414720;
        Hbmjk[p][m] = Hbmjk[m][p];
        Hbmjk[m][q] = ( 4115.*kb[m][r] + 1055.*kb[n][r] + 1055.*kb[p][r] + 1815.*kb[q][r])*elementSize414720;
        Hbmjk[q][m] = Hbmjk[m][p];
        Hbmjk[n][p] = ( 1055.*kb[m][r] +  491.*kb[n][r] +  491.*kb[p][r] +  291.*kb[q][r])*elementSize414720;
        Hbmjk[p][n] = Hbmjk[n][p];
        Hbmjk[n][q] = ( 1055.*kb[m][r] +  491.*kb[n][r] +  291.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[q][n] = Hbmjk[n][q];
        Hbmjk[p][q] = ( 1055.*kb[m][r] +  291.*kb[n][r] +  491.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[q][p] = Hbmjk[p][q];

        for (unsigned j=0; j<nProducts; j++)
        {
          unsigned k = (j+1)%2;
          Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]] + Hbmjk[m][p]*concentrations[p][products[k]] + Hbmjk[m][q]*concentrations[q][products[k]]);
          Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]] + Hbmjk[n][p]*concentrations[p][products[k]] + Hbmjk[n][q]*concentrations[q][products[k]]);
          Hbmj[p] = 0.5*(Hbmjk[p][m]*concentrations[m][products[k]] + Hbmjk[p][n]*concentrations[n][products[k]] + Hbmjk[p][p]*concentrations[p][products[k]] + Hbmjk[p][q]*concentrations[q][products[k]]);
          Hbmj[q] = 0.5*(Hbmjk[q][m]*concentrations[m][products[k]] + Hbmjk[q][n]*concentrations[n][products[k]] + Hbmjk[q][p]*concentrations[p][products[k]] + Hbmjk[q][q]*concentrations[q][products[k]]);

          elementMat[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
          elementMat[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
          elementMat[eq(m,products[0])][var(p,products[j])] -= Hbmj[p];
          elementMat[eq(m,products[0])][var(q,products[j])] -= Hbmj[q];
          elementMat[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
          elementMat[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
          elementMat[eq(m,products[1])][var(p,products[j])] -= Hbmj[p];
          elementMat[eq(m,products[1])][var(q,products[j])] -= Hbmj[q];
          for (unsigned h = 0; h < nReagents; h++)
          {
            elementMat[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
            elementMat[eq(m,reagents[h])][var(n,products[j])] += Hbmj[n];
            elementMat[eq(m,reagents[h])][var(p,products[j])] += Hbmj[p];
            elementMat[eq(m,reagents[h])][var(q,products[j])] += Hbmj[q];
          }

        }
      }
    }
    delete[] products;
    delete[] reagents;
  }
}
//---------------------------------------------------------------------------
void HomReactionTerm_3D_MDC::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
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
  double elementSize414720 = elementSize/414720.;

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
        Hfmjk[m][m] = (17535.*kf[m][r] + 4115.*kf[n][r] + 4115.*kf[p][r] + 4115.*kf[q][r])*elementSize414720;
        Hfmjk[n][n] = ( 1815.*kf[m][r] + 1067.*kf[n][r] +  491.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[p][p] = ( 1815.*kf[m][r] +  491.*kf[n][r] + 1067.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[q][q] = ( 1815.*kf[m][r] +  491.*kf[n][r] +  491.*kf[p][r] + 1067.*kf[q][r])*elementSize414720;
        Hfmjk[m][n] = ( 4115.*kf[m][r] + 1815.*kf[n][r] + 1055.*kf[p][r] + 1055.*kf[q][r])*elementSize414720;
        Hfmjk[n][m] = Hfmjk[m][n];
        Hfmjk[m][p] = ( 4115.*kf[m][r] + 1055.*kf[n][r] + 1815.*kf[p][r] + 1055.*kf[q][r])*elementSize414720;
        Hfmjk[p][m] = Hfmjk[m][p];
        Hfmjk[m][q] = ( 4115.*kf[m][r] + 1055.*kf[n][r] + 1055.*kf[p][r] + 1815.*kf[q][r])*elementSize414720;
        Hfmjk[q][m] = Hfmjk[m][p];
        Hfmjk[n][p] = ( 1055.*kf[m][r] +  491.*kf[n][r] +  491.*kf[p][r] +  291.*kf[q][r])*elementSize414720;
        Hfmjk[p][n] = Hfmjk[n][p];
        Hfmjk[n][q] = ( 1055.*kf[m][r] +  491.*kf[n][r] +  291.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[q][n] = Hfmjk[n][q];
        Hfmjk[p][q] = ( 1055.*kf[m][r] +  291.*kf[n][r] +  491.*kf[p][r] +  491.*kf[q][r])*elementSize414720;
        Hfmjk[q][p] = Hfmjk[p][q];

        for (unsigned j=0; j<nReagents; j++)
        {
          unsigned k = (j+1)%2;
          Hfmj[m] = 0.5*(Hfmjk[m][m]*concentrations[m][reagents[k]] + Hfmjk[m][n]*concentrations[n][reagents[k]] + Hfmjk[m][p]*concentrations[p][reagents[k]] + Hfmjk[m][q]*concentrations[q][reagents[k]]);
          Hfmj[n] = 0.5*(Hfmjk[n][m]*concentrations[m][reagents[k]] + Hfmjk[n][n]*concentrations[n][reagents[k]] + Hfmjk[n][p]*concentrations[p][reagents[k]] + Hfmjk[n][q]*concentrations[q][reagents[k]]);
          Hfmj[p] = 0.5*(Hfmjk[p][m]*concentrations[m][reagents[k]] + Hfmjk[p][n]*concentrations[n][reagents[k]] + Hfmjk[p][p]*concentrations[p][reagents[k]] + Hfmjk[p][q]*concentrations[q][reagents[k]]);
          Hfmj[q] = 0.5*(Hfmjk[q][m]*concentrations[m][reagents[k]] + Hfmjk[q][n]*concentrations[n][reagents[k]] + Hfmjk[q][p]*concentrations[p][reagents[k]] + Hfmjk[q][q]*concentrations[q][reagents[k]]);

          elementJac[eq(m,reagents[0])][var(m,reagents[j])] -= Hfmj[m];
          elementJac[eq(m,reagents[0])][var(n,reagents[j])] -= Hfmj[n];
          elementJac[eq(m,reagents[0])][var(p,reagents[j])] -= Hfmj[p];
          elementJac[eq(m,reagents[0])][var(q,reagents[j])] -= Hfmj[q];
          elementJac[eq(m,reagents[1])][var(m,reagents[j])] -= Hfmj[m];
          elementJac[eq(m,reagents[1])][var(n,reagents[j])] -= Hfmj[n];
          elementJac[eq(m,reagents[1])][var(p,reagents[j])] -= Hfmj[p];
          elementJac[eq(m,reagents[1])][var(q,reagents[j])] -= Hfmj[q];
          for (unsigned h = 0; h < nProducts; h++)
          {
            elementJac[eq(m,products[h])][var(m,reagents[j])] += Hfmj[m];
            elementJac[eq(m,products[h])][var(n,reagents[j])] += Hfmj[n];
            elementJac[eq(m,products[h])][var(p,reagents[j])] += Hfmj[p];
            elementJac[eq(m,products[h])][var(q,reagents[j])] += Hfmj[q];
          }

        }
      }

      // Backward reaction
      if (nProducts == 2)
      {
        Hbmjk[m][m] = (17535.*kb[m][r] + 4115.*kb[n][r] + 4115.*kb[p][r] + 4115.*kb[q][r])*elementSize414720;
        Hbmjk[n][n] = ( 1815.*kb[m][r] + 1067.*kb[n][r] +  491.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[p][p] = ( 1815.*kb[m][r] +  491.*kb[n][r] + 1067.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[q][q] = ( 1815.*kb[m][r] +  491.*kb[n][r] +  491.*kb[p][r] + 1067.*kb[q][r])*elementSize414720;
        Hbmjk[m][n] = ( 4115.*kb[m][r] + 1815.*kb[n][r] + 1055.*kb[p][r] + 1055.*kb[q][r])*elementSize414720;
        Hbmjk[n][m] = Hbmjk[m][n];
        Hbmjk[m][p] = ( 4115.*kb[m][r] + 1055.*kb[n][r] + 1815.*kb[p][r] + 1055.*kb[q][r])*elementSize414720;
        Hbmjk[p][m] = Hbmjk[m][p];
        Hbmjk[m][q] = ( 4115.*kb[m][r] + 1055.*kb[n][r] + 1055.*kb[p][r] + 1815.*kb[q][r])*elementSize414720;
        Hbmjk[q][m] = Hbmjk[m][p];
        Hbmjk[n][p] = ( 1055.*kb[m][r] +  491.*kb[n][r] +  491.*kb[p][r] +  291.*kb[q][r])*elementSize414720;
        Hbmjk[p][n] = Hbmjk[n][p];
        Hbmjk[n][q] = ( 1055.*kb[m][r] +  491.*kb[n][r] +  291.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[q][n] = Hbmjk[n][q];
        Hbmjk[p][q] = ( 1055.*kb[m][r] +  291.*kb[n][r] +  491.*kb[p][r] +  491.*kb[q][r])*elementSize414720;
        Hbmjk[q][p] = Hbmjk[p][q];

        for (unsigned j=0; j<nProducts; j++)
        {
          unsigned k = (j+1)%2;
          Hbmj[m] = 0.5*(Hbmjk[m][m]*concentrations[m][products[k]] + Hbmjk[m][n]*concentrations[n][products[k]] + Hbmjk[m][p]*concentrations[p][products[k]] + Hbmjk[m][q]*concentrations[q][products[k]]);
          Hbmj[n] = 0.5*(Hbmjk[n][m]*concentrations[m][products[k]] + Hbmjk[n][n]*concentrations[n][products[k]] + Hbmjk[n][p]*concentrations[p][products[k]] + Hbmjk[n][q]*concentrations[q][products[k]]);
          Hbmj[p] = 0.5*(Hbmjk[p][m]*concentrations[m][products[k]] + Hbmjk[p][n]*concentrations[n][products[k]] + Hbmjk[p][p]*concentrations[p][products[k]] + Hbmjk[p][q]*concentrations[q][products[k]]);
          Hbmj[q] = 0.5*(Hbmjk[q][m]*concentrations[m][products[k]] + Hbmjk[q][n]*concentrations[n][products[k]] + Hbmjk[q][p]*concentrations[p][products[k]] + Hbmjk[q][q]*concentrations[q][products[k]]);

          elementJac[eq(m,products[0])][var(m,products[j])] -= Hbmj[m];
          elementJac[eq(m,products[0])][var(n,products[j])] -= Hbmj[n];
          elementJac[eq(m,products[0])][var(p,products[j])] -= Hbmj[p];
          elementJac[eq(m,products[0])][var(q,products[j])] -= Hbmj[q];
          elementJac[eq(m,products[1])][var(m,products[j])] -= Hbmj[m];
          elementJac[eq(m,products[1])][var(n,products[j])] -= Hbmj[n];
          elementJac[eq(m,products[1])][var(p,products[j])] -= Hbmj[p];
          elementJac[eq(m,products[1])][var(q,products[j])] -= Hbmj[q];
          for (unsigned h = 0; h < nReagents; h++)
          {
            elementJac[eq(m,reagents[h])][var(m,products[j])] += Hbmj[m];
            elementJac[eq(m,reagents[h])][var(n,products[j])] += Hbmj[n];
            elementJac[eq(m,reagents[h])][var(p,products[j])] += Hbmj[p];
            elementJac[eq(m,reagents[h])][var(q,products[j])] += Hbmj[q];
          }

        }
      }
    }
    delete[] products;
    delete[] reagents;
  }
}
//---------------------------------------------------------------------------

