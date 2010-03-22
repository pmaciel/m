//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------

#include <cstdlib>

#include "ElectrolyteModel_DH.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrolyteModel_DH::ElectrolyteModel_DH(ElectrolyteSolution* electrolyteSolution_)
  : ElectrolyteModel_Ideal (electrolyteSolution_)
{
  conductivityCorrectionFactor = 1.;
  f = new double[nIons];
  Alphas = new double[nIons];
  NormFactors = new double[nIons];
  Eigenvalues = new double[nIons];
  Eigenvectors = new double*[nIons];
  ls = new double*[nIons];
  lm = new double*[nIons];
  dlnf = new double*[nIons];
  for (unsigned i=0; i<nIons; i++) {
    ls[i] = new double[nIons];
    lm[i] = new double[nIons];
    dlnf[i] = new double[nIons];
    Eigenvectors[i] = new double[nIons];
  }
  M = new double[nIons];
}
//---------------------------------------------------------------------------
ElectrolyteModel_DH::~ElectrolyteModel_DH()
{
  for (unsigned i=0; i<nIons; i++) {
    delete[] ls[i];
    delete[] lm[i];
    delete[] dlnf[i];
    delete[] Eigenvectors[i];
  }
  delete[] f;
  delete[] Alphas;
  delete[] NormFactors;
  delete[] Eigenvalues;
  delete[] Eigenvectors;
  delete[] ls;
  delete[] lm;
  delete[] dlnf;
  delete[] M;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_DH::init(bool verbose)
{
  this->verbose = verbose;
  T = electrolyteSolution->getSolutionTemperature();
  viscosity = electrolyteSolution->getSolventDynamicViscosity();
  dielectricConstant = electrolyteSolution->getSolventDielectricConstant();
  density = electrolyteSolution->getSolutionDensity();
  for (unsigned i=0; i<nIons; i++)
  {
    z[i] = electrolyteSolution->getIonChargeNumber(i);
    c[i] = electrolyteSolution->getIonConcentration(i);
    D[i] = electrolyteSolution->getIonDiffusionConstant(i);
    M[i] = electrolyteSolution->getIonMolarMass(i);
  }

  LB = e_CONST*e_CONST/(4.*M_PI*kB_CONST*T*dielectricConstant);

  // calculate parameters
  calcConcentrationParameters();
  calcEigenproblem();

  std::ofstream output;
  if (verbose)
  {
    output.setf(std::ios::scientific);
    output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
    double totalConcentration = c[0]+c[1];
    if (nIons == 3) totalConcentration += 2*c[2];
    output << totalConcentration << '\t';
    for (unsigned i=0; i<nIons; i++)
      output << c[i] << '\t';
    output << HIDL << '\t';
    double osmoticPressureC = -HIDL*HIDL*HIDL/(3.*M_PI);
    output << osmoticPressureC/(totalConcentration*NA_CONST) + 1. << '\t';
    output.close();
  }

  for (unsigned i=0; i<nIons; i++) {
    // Activity coefficients in the McMillan-Mayer reference frame
    // THEY SHOULD BE CONVERTED TO THE LEWIS-RANDALL REFERENCE FRAME!!!
    double lnf = -LB*HIDL*z[i]*z[i];
    f[i] = exp(lnf);

    if (verbose)
    {
      output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
      output << f[i] << '\t';
      output.close();
    }

    for (unsigned j=i; j<nIons; j++) {
      // Partial derivatives of the activity coefficients --> d(lnfi)/d(cj) = d(lnfj)/d(ci)
      dlnf[i][j] = -M_PI*NA_CONST*LB*LB*z[i]*z[i]*z[j]*z[j]/(2.*HIDL);
      dlnf[j][i] = dlnf[i][j];

      // Lij in the solvent-fixed reference frame --> Lij=Lji
      double EC = calcElectrophoreticCorrection(i,j);
      double RC = calcRelaxationCorrection(i,j);
      ls[i][j] = (kroneck(i,j)*c[i]*D[i]/(R_CONST*T) + EC + RC)*conductivityCorrectionFactor;
      ls[j][i] = ls[i][j];
    }
  }

  double contrib3 = 0;
  for (unsigned k=0; k<nIons; k++)
  {
    for (unsigned l=0; l<nIons; l++)
    {
      contrib3 += M[k]*M[l]*ls[k][l];
    }
  }
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned j=i; j<nIons; j++)
    {
      // Lij in the mass-fixed reference frame
      lm[i][j] = ls[i][j];
      double contrib1 = 0;
      double contrib2 = 0;
      for (unsigned k=0; k<nIons; k++)
      {
        contrib1 += M[k]*ls[k][j];
        contrib2 += M[k]*ls[i][k];
      }
      lm[i][j] = ls[i][j]
                 - c[i]/density*contrib1
                 - c[j]/density*contrib2
                 + c[i]*c[j]/(density*density)*contrib3;
      lm[j][i] = lm[i][j];
    }
  }

  if (verbose)
  {
    output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
    output << std::endl;
    output.close();
  }
}
//---------------------------------------------------------------------------
// Concentration parameters (2x ionic strength (I), half inverse Debye length (HIDL) and limiting conductivity (D))
void ElectrolyteModel_DH::calcConcentrationParameters()
{
  I = 0.;
  cond0 = 0.;
  Mscs = density;
  for (unsigned i=0; i<nIons; i++) {
    I += z[i]*z[i]*c[i];
    cond0 += z[i]*z[i]*c[i]*D[i];
    Mscs -= M[i]*c[i];
  }
  if (I <= 0) errorZero("ElectrolyteModel_DH.cpp","calcConcentrationParameters");
  if (Mscs <= 0) errorZero("ElectrolyteModel_DH.cpp","calcConcentrationParameters");
  HIDL = sqrt(M_PI*LB*I*NA_CONST);
}
//---------------------------------------------------------------------------
// Electrophoretic correction for the Onsager coefficients
double ElectrolyteModel_DH::calcElectrophoreticCorrection(unsigned i, unsigned j) const
{
  if (verbose)
  {
    std::ofstream output;
    output.setf(std::ios::scientific);
    output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
    output << -LB/(3.*viscosity*HIDL)*z[i]*z[j]*c[j] << '\t';
    output.close();
  }

  return -LB/(3.*viscosity*HIDL)*z[i]*z[j]*c[i]*c[j];
}
//---------------------------------------------------------------------------
// Electrophoretic correction for the conductivity
/*double ElectrolyteModel_DH::calcElectrophoreticCorrection(unsigned i) const
{
  return -z[i]*HIDL/(3.*M_PI*NA_CONST*viscosity);
}*/
//---------------------------------------------------------------------------
// Eigenproblem
void ElectrolyteModel_DH::calcEigenproblem()
{
  std::vector< std::vector<unsigned> > DList;

  // Ordered list of D
  std::vector<unsigned> DIndex;
  std::vector<unsigned>::iterator it;
  DIndex.push_back(0);
  for (unsigned i=1; i<nIons; i++)
  {
    // Scroll through already ordered list of D...
    for (it = DIndex.begin(); it != DIndex.end(); it++)
    {
      //... until you find a D that is greater. Insert the new D before this one.
      if (D[*it] >= D[i])
      {
        DIndex.insert(it,i);
        break;
      }
      //... or until you are at the end of the list. Append the new D at the end of the list.
      else if (it+1 == DIndex.end())
      {
        DIndex.push_back(i);
        break; // NECESSARY!!!
      }
    }
  }

  // Ordered list of confluences of D
  std::vector<unsigned> DValue;
  std::vector<unsigned>::iterator begin,end;
  for (begin = DIndex.begin(); begin != DIndex.end(); begin++)
  {
    for (end = begin+1; end != DIndex.end(); end++)
    {
      if (D[*end] > D[*begin]) break;
    }
    DValue.assign(begin,end);
    DList.push_back(DValue);
    begin = end-1;
  }

  for (unsigned i=0; i<nIons; i++)
  {
    Alphas[i] = 0.;
    NormFactors[i] = 0.;
    Eigenvalues[i] = 0.;
    for (unsigned j=0; j<nIons; j++)
    {
      Eigenvectors[i][j] = 0.;
    }
  }

  unsigned p = 0;
  for (unsigned m=0; m<DList.size(); m++) {
    // If it is the first ion
    if (p == 0) {
      //Alphas[p] = 0.;
      for (unsigned i=0; i<nIons; i++) {
        NormFactors[p] += z[i]*z[i]*c[i]/D[i];
      }
      Eigenvalues[p] = 1.;
      NormFactors[p] = sqrt(cond0/NormFactors[p]);
      for (unsigned i=0; i<nIons; i++) {
        Eigenvectors[p][i] = z[i]*z[i]*c[i]*NormFactors[p];
        //std::cout << "chi(" << p << "," << i << ") = " << Eigenvectors[p][i] << endl;
      }
    }

    // Else
    else {
      double AlphaOld,AlphaMin,AlphaMax;
      AlphaMin = D[DList[m-1][0]]*D[DList[m-1][0]];
      AlphaMax = D[DList[m][0]]*D[DList[m][0]];

      // Starting value is the middle of the interval [D(p-1)^2,D(p)^2]
      Alphas[p] = (AlphaMin + AlphaMax)/2.;

      // Newton iterations
      unsigned iter = 0;
      const unsigned iterMax = 20;
      bool TracerIon;
      unsigned TracerIndex = 0;
      do {
        TracerIon = false;
        AlphaOld = Alphas[p];
        Alphas[p] = AlphaOld - calcAlphaFunction(AlphaOld)/calcAlphaFunctionDerivative(AlphaOld);
        iter++;
        // ----------------------------------------------------------------------------------------	//
        // Avoid leaving the interval [D(p-1)^2,D(p)^2] or you will not find the right solution!	//
        // Also avoid getting too close to a boundary, or you risk division by zero (can happen if	//
        // c(p) is low or if the interval is small)!												//
        // ----------------------------------------------------------------------------------------	//
        if (Alphas[p]/AlphaMin-1. < 1e-9) {
          TracerIon = true;
          TracerIndex = DList[m-1][0];
          Alphas[p] = (AlphaMin + AlphaOld)/2.;
        }
        else if (1.-Alphas[p]/AlphaMax < 1e-9) {
          TracerIon = true;
          TracerIndex = DList[m][0];
          Alphas[p] = (AlphaMax + AlphaOld)/2.;
        }
        // ----------------------------------------------------------------------------------------	//
        // The D are not known with more than 4 significant digits.									//
        // The relative difference between two different D is therefore minimally 10^-3.			//
        // The relative difference between two different D^2 is therefore minimally 10^-6.			//
        // The criterion for "too close to a boundary" must be smaller than this (e.g. 10^-9),		//
        // otherwise you can be too close to both boundaries at the same time!						//
        // ----------------------------------------------------------------------------------------	//
      } while ((fabs(1.-Alphas[p]/AlphaOld) > 1e-9) && (iter < iterMax));
      // Stop when relative change is small enough.

      // If it is a tracer ion
      if (TracerIon) {
        Alphas[p] = D[TracerIndex]*D[TracerIndex];
        for (unsigned i=0; i<nIons; i++) {
          Eigenvalues[p] += z[i]*z[i]*c[i]*D[i]/(D[i] + D[TracerIndex]);
        }
        Eigenvalues[p] /= I;
        NormFactors[p] = 0.;
        //Eigenvectors[p].assign(nIons,0.);
        for (unsigned i=0; i<nIons; i++)
          Eigenvectors[p][i] = 0.;
        Eigenvectors[p][TracerIndex] = sqrt(cond0*z[TracerIndex]*z[TracerIndex]*c[TracerIndex]*D[TracerIndex]);
        //std::cout << "chi(" << p << "," << TracerIndex << ") = " << Eigenvectors[p][TracerIndex] << endl;
      }

      // Else
      else {
        for (unsigned i=0; i<nIons; i++) {
          Eigenvalues[p] += z[i]*z[i]*c[i]*D[i]/(D[i] + sqrt(Alphas[p]));
          double factor = z[i]*D[i]/(D[i]*D[i] - Alphas[p]); // Risk of division by zero eliminated
          NormFactors[p] += c[i]*D[i]*factor*factor;
        }
        Eigenvalues[p] /= I;
        NormFactors[p] = sqrt(cond0/NormFactors[p]);
        for (unsigned i=0; i<nIons; i++) {
          Eigenvectors[p][i] = z[i]*z[i]*c[i]*NormFactors[p]*D[i]*D[i]/(D[i]*D[i] - Alphas[p]);
          //std::cout << "chi(" << p << "," << i << ") = " << Eigenvectors[p][i] << endl;
        }
      }
    }
    p++;

    // If it is a confluence
    if (DList[m].size() > 1) {
      double AlphaValue = D[DList[m][0]]*D[DList[m][0]];
      double Eigenvalue = 0.;
      for (unsigned i=0; i<nIons; i++) {
        Eigenvalue += z[i]*z[i]*c[i]*D[i]/(D[i] + D[DList[m][0]]);
      }
      Eigenvalue /= I;
      for (unsigned n=1; n<DList[m].size(); n++) {
        Alphas[p] = AlphaValue;
        Eigenvalues[p] = Eigenvalue;
        //NormFactors[p] = 0.;
        for (unsigned i=0; i<nIons; i++)
          Eigenvectors[p][i] = 0.;

        double Tp = 0.;
        for (unsigned r=0; r<n; r++)
          Tp += z[DList[m][r]]*z[DList[m][r]]*c[DList[m][r]]*D[DList[m][r]];
        double tp = z[DList[m][n]]*z[DList[m][n]]*c[DList[m][n]]*D[DList[m][n]];
        double Eigenvector = -sqrt(cond0*tp/(Tp*(Tp+tp)));
        for (unsigned r=0; r<n; r++) {
          if (z[DList[m][r]] != 0)
            Eigenvectors[p][DList[m][r]] = z[DList[m][r]]*z[DList[m][r]]*c[DList[m][r]]*D[DList[m][r]]*Eigenvector;
          //std::cout << "chi(" << p << "," << DList[m][r] << ") = " << Eigenvectors[p][DList[m][r]] << endl;
        }
        if (z[DList[m][n]] != 0)
          Eigenvectors[p][DList[m][n]] = sqrt(cond0*Tp*tp/(Tp+tp));
        //std::cout << "chi(" << p << "," << DList[m][n] << ") = " << Eigenvectors[p][DList[m][n]] << endl;
        p++;
      }
    }
  }

  /*ofstream CHECK("Eigenvalues.out");
  CHECK.precision(16);
  double ti;
  CHECK << "\t\t";
  for (unsigned i=0; i<nIons; i++) {
    ti = z[i]*z[i]*c[i]*D[i]/cond0;
    CHECK << '\t' << ti;
  }
  CHECK << endl;
  for (unsigned r=0; r<nIons; r++) {
    CHECK << Alphas[r] << '\t' << Eigenvalues[r] << '\t' << NormFactors[r];
    for (unsigned i=0; i<nIons; i++) {
      CHECK << '\t' << Eigenvectors[r][i];
    }
    CHECK << endl;
  }*/
  //double ti;
  //std::cout << "\t\t";
  //for (unsigned i=0; i<nIons; i++) {
  //	ti = z[i]*z[i]*c[i]*D[i]/cond0;
  //	std::cout << '\t' << ti;
  //}
  //std::cout << std::endl;
  //for (unsigned r=0; r<nIons; r++) {
  //	std::cout << Alphas[r] << '\t' << Eigenvalues[r] << '\t' << NormFactors[r];
  //	for (unsigned i=0; i<nIons; i++) {
  //		std::cout << '\t' << Eigenvectors[r][i];
  //	}
  //	std::cout << std::endl;
  //}
}
//---------------------------------------------------------------------------
// Alpha function
double ElectrolyteModel_DH::calcAlphaFunction(double a) const
{
  double aFunct = 0.;
  for (unsigned i=0; i<nIons; i++) {
    aFunct += z[i]*z[i]*c[i]*D[i]/(D[i]*D[i] - a);
  }
  return aFunct;
}
//---------------------------------------------------------------------------
// Derivative Alpha function
double ElectrolyteModel_DH::calcAlphaFunctionDerivative(double a) const
{
  double daFunct = 0.;
  for (unsigned i=0; i<nIons; i++) {
    daFunct += z[i]*z[i]*c[i]*D[i]/((D[i]*D[i] - a)*(D[i]*D[i] - a));
  }
  if (daFunct <= 0) errorZero("ElectrolyteModel_DH.cpp","calcAlphaFunctionDerivative");
  return daFunct;
}
//---------------------------------------------------------------------------
// Relaxation correction for the Onsager coefficients
double ElectrolyteModel_DH::calcRelaxationCorrection(unsigned i, unsigned j) const
{
  double RCorr = 0.;
  // no need to include p=0, because it produces a zero term
  for (unsigned p=1; p<nIons; p++) {
    RCorr += (1.-sqrt(Eigenvalues[p]))*Eigenvectors[p][i]*Eigenvectors[p][j];
  }

  if (verbose)
  {
    std::ofstream output;
    output.setf(std::ios::scientific);
    output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
    output << -2.*HIDL*LB/(3*cond0*R_CONST*T)*RCorr/c[i] << '\t';
    output.close();
  }

  return -2.*HIDL*LB/(3*cond0*R_CONST*T)*RCorr;
}
//---------------------------------------------------------------------------
// Relaxation correction for the conductivity
/*double ElectrolyteModel_DH::calcRelaxationCorrection(unsigned i) const
{
  double RCorr = 0.;
  // no need to include p=0, because it produces a zero term
  for (unsigned j=0; j<nIons; j++) {
    double contribj = 0.;
    for (unsigned p=1; p<nIons; p++) {
      contribj += (1.-sqrt(Eigenvalues[p]))*Eigenvectors[p][i]*Eigenvectors[p][j];
    }
    RCorr += z[j]*contribj;
  }
  return -2.*LB*HIDL/(3*cond0*R_CONST*T)*RCorr;
}*/
//---------------------------------------------------------------------------
// Dij
double ElectrolyteModel_DH::calcDiffusionFactor(unsigned i, unsigned j) const
{
  // POSSIBLE ERROR IN THIS FUNCTION (mistake in the indices?)
  double contrib = 0.;
  for (unsigned l=0; l<nIons; l++)
  {
    contrib += c[l]*dlnf[l][j];
  }
  double D = lm[i][j]/c[j];
  for (unsigned k=0; k<nIons; k++)
  {
    D += lm[i][k]*(dlnf[k][j] + M[k]/Mscs*(1. + contrib));
  }
  D *= R_CONST*T;
  return D;
}
//---------------------------------------------------------------------------
// Wi = zi*F*ui
double ElectrolyteModel_DH::calcMigrationFactor(unsigned i) const
{
  double W = 0.;
  for (unsigned j=0; j<nIons; j++)
  {
    W += z[j]*lm[i][j];
  }
  W *= F_CONST/c[i];
  return W;
}
//---------------------------------------------------------------------------
// Conductivity
double ElectrolyteModel_DH::calcConductivity() const
{
  double kappa = 0.;
  for (unsigned i=0; i<nIons; i++)
    for (unsigned j=0; j<nIons; j++)
      kappa += z[i]*z[j]*ls[i][j];
  return kappa*F_CONST*F_CONST;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_DH::calcOnsagerCoefficient_SF(unsigned i, unsigned j) const
{
  return ls[i][j];
}
//---------------------------------------------------------------------------
double ElectrolyteModel_DH::calcActivityCoefficient_MM(unsigned i) const
{
  return f[i];
}
//---------------------------------------------------------------------------
double ElectrolyteModel_DH::calcActivityCoefficientDerivative_MM(unsigned i, unsigned j) const
{
  return dlnf[i][j];
}
//---------------------------------------------------------------------------
double ElectrolyteModel_DH::calcOsmoticCoefficient_MM(double totalConcentration) const
{
  // Coulomb
  double osmoticPressureC = -HIDL*HIDL*HIDL/(3.*M_PI);

  // Ideal
  double osmoticPressureId = 0.;
  for (unsigned i=0; i<nIons; i++) {
    osmoticPressureId += c[i];
  }
  osmoticPressureId *= NA_CONST;

  return (osmoticPressureId + osmoticPressureC)/(totalConcentration*NA_CONST);
}
//---------------------------------------------------------------------------
double ElectrolyteModel_DH::calcMeanActivityCoefficient_MM(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const
{
  double fMean = pow(pow(f[0]*c[0]/sCation,(int)sCation)*pow(f[1]*c[1]/sAnion,(int)sAnion),1./(sCation+sAnion));
  return fMean/electrolyteConcentration;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_DH::calcBinaryOnsagerCoefficients_SF(int* stoichCation, int* stoichAnion)
{
  Ls00 = ls[0][0];
  Ls11 = ls[1][1];
  Ls01 = ls[0][1];
  unsigned nAssociatedIons = nIons-2;
  for (unsigned i=0; i<nAssociatedIons; i++)
  {
    unsigned i2 = i+2;
    double sumCat = 0.;
    double sumAn = 0.;
    for (unsigned j=0; j<nAssociatedIons; j++)
    {
      unsigned j2 = j+2;
      sumCat += stoichCation[j]*ls[j2][i2];
      sumAn += stoichAnion[j]*ls[j2][i2];
    }
    Ls00 -= stoichCation[i]*(ls[i2][0]+ls[0][i2]-sumCat);
    Ls11 -= stoichAnion[i]*(ls[i2][1]+ls[1][i2]-sumAn);
    Ls01 -= stoichCation[i]*ls[i2][1]+stoichAnion[i]*(ls[0][i2]-sumCat);
  }
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void ElectrolyteModel_DH::errorZero(const std::string &cppfile, const std::string &function) const
{
  std::cout << "ERROR IN " << cppfile << ".\nA PARAMETER IN " << function << " IS <= 0.\
  \nPROBABLY SOME ION CONCENTRATIONS ARE NEGATIVE OR UNREALISTICALLY HIGH." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------

