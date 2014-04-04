//---------------------------------------------------------------------------

#include "Fitter.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "DatFileReader.h"
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Fitter::Fitter(const std::string &name_)
: name(name_)
{
  electrolytesFile = name + ".electrolytes";
  ionsFile = name + ".ions";
  homReactionsFile = name + ".homreactions";

  readElectrolytes();
  readIons();
  readHomReactions();
}
//---------------------------------------------------------------------------
Fitter::~Fitter()
{
}
//---------------------------------------------------------------------------


//--- READ METHODS ----------------------------------------------------------
void Fitter::readElectrolytes()
{
  std::string label;

  DatFileReader datFile(electrolytesFile);

  datFile.readScalar("[maxIter]",maxIter);
  nExperimentalValues = 0;
  nElectrolytes = datFile.readMultipleVector_nVectors("[nElectrolytes]");
  electrolytes = new Electrolyte*[nElectrolytes];
  for (unsigned e=0; e<nElectrolytes; e++) {
    datFile.readMultipleVector("[nElectrolytes]", e, "<label>", label);
    electrolytes[e] = new Electrolyte(label);
    unsigned sCation, sAnion;
    datFile.readMultipleVector("[nElectrolytes]", e, "<sCation>", sCation);
    datFile.readMultipleVector("[nElectrolytes]", e, "<sAnion>", sAnion);
    electrolytes[e]->setStoichCoefficients(sCation,sAnion);
    double cMax;
    datFile.readMultipleVector("[nElectrolytes]", e, "<cMax>", cMax);
    electrolytes[e]->setMaximumConcentration(cMax);
    nExperimentalValues += electrolytes[e]->getNOsmoticCoefficients()
      + electrolytes[e]->getNActivityCoefficients()
      + electrolytes[e]->getNConductivities()
      + electrolytes[e]->getNTransportNumbers()
      + electrolytes[e]->getNDiffusionCoefficients();
  }
}
//---------------------------------------------------------------------------
void Fitter::readIons()
{
  double doubleValue;
  std::string stringValue;
  bool boolValue;

  DatFileReader datFile(ionsFile);

  nIons = datFile.readMultipleVector_nVectors("[nIons]");
  ions = new IonFit*[nIons];
  for (unsigned i=0; i<nIons; i++) {
    datFile.readMultipleVector("[nIons]", i, "<label>", stringValue);
    ions[i] = new IonFit(stringValue,electrolytes, nElectrolytes);
    datFile.readMultipleVector("[nIons]", i, "<d>", doubleValue);
    ions[i]->setDiameter(doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<dFit>", boolValue);
    if (boolValue==true) {
      datFile.readMultipleVector("[nIons]", i, "<dMin>", doubleValue);
      ions[i]->setMinimumDiameter(doubleValue);
      datFile.readMultipleVector("[nIons]", i, "<dMax>", doubleValue);
      ions[i]->setMaximumDiameter(doubleValue);
      diametersToFit.push_back(i);
    }

    datFile.readMultipleVector("[nIons]", i, "<D>", doubleValue);
    ions[i]->setDiffusionConstant(doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<DFit>", boolValue);
    if (boolValue==true) {
      datFile.readMultipleVector("[nIons]", i, "<DMin>", doubleValue);
      ions[i]->setMinimumDiffusionConstant(doubleValue);
      datFile.readMultipleVector("[nIons]", i, "<DMax>", doubleValue);
      ions[i]->setMaximumDiffusionConstant(doubleValue);
      diffusionConstantsToFit.push_back(i);
    }
    datFile.readMultipleVector("[nIons]", i, "<M>", doubleValue);
    ions[i]->setMolarMass(doubleValue);
  }
  nDiametersToFit = unsigned(diametersToFit.size());
  nDiffusionConstantsToFit = unsigned(diffusionConstantsToFit.size());
}
//---------------------------------------------------------------------------
void Fitter::readHomReactions()
{
  double doubleValue;
  std::string stringValue;
  bool boolValue;

  DatFileReader datFile(homReactionsFile);

  nHomReactions = datFile.readMultipleVector_nVectors("[nHomReactions]");
  homReactions = new HomReactionFit*[nHomReactions];
  for (unsigned r=0; r<nHomReactions; r++) {
    datFile.readMultipleVector("[nHomReactions]", r, "<label>", stringValue);
    homReactions[r] = new HomReactionFit(stringValue,electrolytes,nElectrolytes);
    datFile.readMultipleVector("[nHomReactions]", r, "<K>", doubleValue);
    homReactions[r]->setEquilibriumConstant(doubleValue);
    datFile.readMultipleVector("[nHomReactions]", r, "<KFit>", boolValue);
    if (boolValue==true) {
      datFile.readMultipleVector("[nHomReactions]", r, "<KMin>", doubleValue);
      homReactions[r]->setMinimumEquilibriumConstant(doubleValue);
      datFile.readMultipleVector("[nHomReactions]", r, "<KMax>", doubleValue);
      homReactions[r]->setMaximumEquilibriumConstant(doubleValue);
      equilibriumConstantsToFit.push_back(r);
    }
  }
  nEquilibriumConstantsToFit = unsigned(equilibriumConstantsToFit.size());
}
//---------------------------------------------------------------------------


//--- SOLVE METHODS ---------------------------------------------------------
void Fitter::solve ()
{
  nParametersToFit = nDiametersToFit + nDiffusionConstantsToFit + nEquilibriumConstantsToFit;

  std::cout << "Number of experimental values = " << nExperimentalValues << std::endl;
  std::cout << "Number of parameters to fit = " << nParametersToFit << std::endl;

  // Initialization
  R = new double[nExperimentalValues];
  for (unsigned m=0; m<nExperimentalValues; m++)
    R[m] = 0.;

  JtR = new double[nParametersToFit];
  for (unsigned n=0; n<nParametersToFit; n++)
    JtR[n] = 0.;

  J = new double*[nExperimentalValues];
  for (unsigned m=0; m<nExperimentalValues; m++)
  {
    J[m] = new double[nParametersToFit];
    for (unsigned n=0; n<nParametersToFit; n++)
      J[m][n] = 0.;
  }

  JtJ = new double*[nParametersToFit];
  for (unsigned n=0; n<nParametersToFit; n++)
  {
    JtJ[n] = new double[nParametersToFit];
    for (unsigned p=0; p<nParametersToFit; p++)
      JtJ[n][p] = 0.;
  }


  //std::ofstream CHECK("Controle.xls");
  //CHECK.precision(12);

  // Newton iterations
  unsigned iter = 0;
  //bool converged = false;
  //double concentration;
  //double T;
  do {
    unsigned row = 0;
    residu = 0.;
    // Fill in residual vector
    for (unsigned e=0; e<nElectrolytes; e++)
    {
      double scaleOsmoticCoefficient = sqrt((double)electrolytes[e]->getNOsmoticCoefficients())*calcLimitingOsmoticCoefficient(e);
      double scaleActivityCoefficient = sqrt((double)electrolytes[e]->getNActivityCoefficients())*calcLimitingActivityCoefficient(e);
      double scaleConductivity = sqrt((double)electrolytes[e]->getNConductivities())*calcLimitingConductivity(e);
      double scaleTransportNumber = sqrt((double)electrolytes[e]->getNTransportNumbers())*calcLimitingTransportNumber(e);
      double scaleDiffusionCoefficient = sqrt((double)electrolytes[e]->getNDiffusionCoefficients())*calcLimitingDiffusionCoefficient(e);

      /*std:: cout << scaleOsmoticCoefficient
        << '\t' << scaleActivityCoefficient
        << '\t' << scaleConductivity
        << '\t' << scaleTransportNumber
        << '\t' << scaleDiffusionCoefficient;
      char ch;
      std::cin >> ch;*/

      //std:: cout << "\nOsmotic coefficients\n";
      // Osmotic coefficients
      for (unsigned t=0; t<electrolytes[e]->getNOsmoticCoefficients(); t++)
      {
        R[row] = electrolytes[e]->calcOsmoticCoefficient(t) - electrolytes[e]->getOsmoticCoefficient(t);
        R[row] /= scaleOsmoticCoefficient;
        //std::cout << R[row] << '\n';
        residu += R[row]*R[row];
        row++;
      }
      /*std::ofstream output;
      output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
      output << std::endl;
      output.close();
      exit(1);*/

      //std:: cout << "\nActivity coefficients\n";
      // Activity coefficients
      for (unsigned t=0; t<electrolytes[e]->getNActivityCoefficients(); t++)
      {
        R[row] = electrolytes[e]->calcActivityCoefficient(t) - electrolytes[e]->getActivityCoefficient(t);
        R[row] /= scaleActivityCoefficient;
        //std::cout << R[row] << '\n';
        residu += R[row]*R[row];
        row++;
      }

      //std:: cout << "\nConductivities\n";
      // Conductivities
      for (unsigned t=0; t<electrolytes[e]->getNConductivities(); t++)
      {
        R[row] = electrolytes[e]->calcConductivity(t) - electrolytes[e]->getConductivity(t);
        R[row] /= scaleConductivity;
        //std::cout << R[row] << '\n';
        residu += R[row]*R[row];
        row++;
      }

      //std:: cout << "\nTransport numbers\n";
      // ElectrolyteModel numbers
      for (unsigned t=0; t<electrolytes[e]->getNTransportNumbers(); t++)
      {
        R[row] = electrolytes[e]->calcTransportNumber(t) - electrolytes[e]->getTransportNumber(t);
        R[row] /= scaleTransportNumber;
        //std::cout << R[row] << '\n';
        residu += R[row]*R[row];
        row++;
      }

      //std:: cout << "\nDiffusion coefficients\n";
      // Diffusion coefficients
      for (unsigned t=0; t<electrolytes[e]->getNDiffusionCoefficients(); t++)
      {
        R[row] = electrolytes[e]->calcDiffusionCoefficient(t) - electrolytes[e]->getDiffusionCoefficient(t);
        R[row] /= scaleDiffusionCoefficient;
        //std::cout << R[row] << '\n';
        residu += R[row]*R[row];
        row++;
      }
    }
    residu = sqrt(residu);
    std::cout << "Iteration " << iter << "\nResidu = " << residu << std::endl;

    //CHECK << "R :\n";
    //for (unsigned m=0; m<nExperimentalValues; m++) {
    //  CHECK << R[m] << std::endl;
    //  std::cout << R[m] << std::endl;
    //}
    //CHECK << std::endl;
    //std::cout << std::endl;
    //std::cin.get();

    unsigned col;
    // Fill in Jacobian matrix
    for (unsigned i=0; i<nDiametersToFit; i++)
    {
      col = i;
      for (unsigned e=0; e<ions[diametersToFit[i]]->getNElectrolytes(); e++)
      {
        unsigned f = ions[diametersToFit[i]]->getElectrolytes(e);
        double scaleOsmoticCoefficient = sqrt((double)electrolytes[f]->getNOsmoticCoefficients())*calcLimitingOsmoticCoefficient(f);
        double scaleActivityCoefficient = sqrt((double)electrolytes[f]->getNActivityCoefficients())*calcLimitingActivityCoefficient(f);
        double scaleConductivity = sqrt((double)electrolytes[f]->getNConductivities())*calcLimitingConductivity(f);
        double scaleTransportNumber = sqrt((double)electrolytes[f]->getNTransportNumbers())*calcLimitingTransportNumber(f);
        double scaleDiffusionCoefficient = sqrt((double)electrolytes[f]->getNDiffusionCoefficients())*calcLimitingDiffusionCoefficient(f);

        for (unsigned t=0; t<electrolytes[f]->getNOsmoticCoefficients(); t++)
        {
          row = calcRowOsmoticCoefficient(f,t);
          J[row][col] = ions[diametersToFit[i]]->calcOsmoticCoefficientDerivativeDiameter(e,t);
          J[row][col] /= scaleOsmoticCoefficient;
        }
        for (unsigned t=0; t<electrolytes[f]->getNActivityCoefficients(); t++)
        {
          row = calcRowActivityCoefficient(f,t);
          J[row][col] = ions[diametersToFit[i]]->calcActivityCoefficientDerivativeDiameter(e,t);
          J[row][col] /= scaleActivityCoefficient;
        }
        for (unsigned t=0; t<electrolytes[f]->getNConductivities(); t++)
        {
          row = calcRowConductivity(f,t);
          J[row][col] = ions[diametersToFit[i]]->calcConductivityDerivativeDiameter(e,t);
          J[row][col] /= scaleConductivity;
        }
        for (unsigned t=0; t<electrolytes[f]->getNTransportNumbers(); t++)
        {
          row = calcRowTransportNumber(f,t);
          J[row][col] = ions[diametersToFit[i]]->calcTransportNumberDerivativeDiameter(e,t);
          J[row][col] /= scaleTransportNumber;
        }
        for (unsigned t=0; t<electrolytes[f]->getNDiffusionCoefficients(); t++)
        {
          row = calcRowDiffusionCoefficient(f,t);
          J[row][col] = ions[diametersToFit[i]]->calcDiffusionCoefficientDerivativeDiameter(e,t);
          J[row][col] /= scaleDiffusionCoefficient;
        }
      }
    }
    //for (unsigned m=0; m<nExperimentalValues; m++) {
    //  for (unsigned n=0; n<nParametersToFit; n++) {
    //    std::cout << J[m][n] << '\t';
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //std::cin.get();

    for (unsigned i=0; i<nDiffusionConstantsToFit; i++)
    {
      col = nDiametersToFit + i;
      for (unsigned e=0; e<ions[diffusionConstantsToFit[i]]->getNElectrolytes(); e++)
      {
        unsigned f = ions[diffusionConstantsToFit[i]]->getElectrolytes(e);
        double scaleOsmoticCoefficient = sqrt((double)electrolytes[f]->getNOsmoticCoefficients())*calcLimitingOsmoticCoefficient(f);
        double scaleActivityCoefficient = sqrt((double)electrolytes[f]->getNActivityCoefficients())*calcLimitingActivityCoefficient(f);
        double scaleConductivity = sqrt((double)electrolytes[f]->getNConductivities())*calcLimitingConductivity(f);
        double scaleTransportNumber = sqrt((double)electrolytes[f]->getNTransportNumbers())*calcLimitingTransportNumber(f);
        double scaleDiffusionCoefficient = sqrt((double)electrolytes[f]->getNDiffusionCoefficients())*calcLimitingDiffusionCoefficient(f);

        for (unsigned t=0; t<electrolytes[f]->getNOsmoticCoefficients(); t++)
        {
          row = calcRowOsmoticCoefficient(f,t);
          J[row][col] = ions[diffusionConstantsToFit[i]]->calcOsmoticCoefficientDerivativeDiffusionConstant(e,t);
          J[row][col] /= scaleOsmoticCoefficient;
        }
        for (unsigned t=0; t<electrolytes[f]->getNActivityCoefficients(); t++)
        {
          row = calcRowActivityCoefficient(f,t);
          J[row][col] = ions[diffusionConstantsToFit[i]]->calcActivityCoefficientDerivativeDiffusionConstant(e,t);
          J[row][col] /= scaleActivityCoefficient;
        }
        for (unsigned t=0; t<electrolytes[f]->getNConductivities(); t++)
        {
          row = calcRowConductivity(f,t);
          J[row][col] = ions[diffusionConstantsToFit[i]]->calcConductivityDerivativeDiffusionConstant(e,t);
          J[row][col] /= scaleConductivity;
        }
        for (unsigned t=0; t<electrolytes[f]->getNTransportNumbers(); t++)
        {
          row = calcRowTransportNumber(f,t);
          J[row][col] = ions[diffusionConstantsToFit[i]]->calcTransportNumberDerivativeDiffusionConstant(e,t);
          J[row][col] /= scaleTransportNumber;
        }
        for (unsigned t=0; t<electrolytes[f]->getNDiffusionCoefficients(); t++)
        {
          row = calcRowDiffusionCoefficient(f,t);
          J[row][col] = ions[diffusionConstantsToFit[i]]->calcDiffusionCoefficientDerivativeDiffusionConstant(e,t);
          J[row][col] /= scaleDiffusionCoefficient;
        }
      }
    }
    //for (unsigned m=0; m<nExperimentalValues; m++) {
    //  for (unsigned n=0; n<nParametersToFit; n++) {
    //    std::cout << J[m][n] << '\t';
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //std::cin.get();

    for (unsigned i=0; i<nEquilibriumConstantsToFit; i++)
    {
      col = nDiametersToFit + nDiffusionConstantsToFit + i;
      unsigned f = homReactions[equilibriumConstantsToFit[i]]->getElectrolyte();
      double scaleOsmoticCoefficient = sqrt((double)electrolytes[f]->getNOsmoticCoefficients())*calcLimitingOsmoticCoefficient(f);
      double scaleActivityCoefficient = sqrt((double)electrolytes[f]->getNActivityCoefficients())*calcLimitingActivityCoefficient(f);
      double scaleConductivity = sqrt((double)electrolytes[f]->getNConductivities())*calcLimitingConductivity(f);
      double scaleTransportNumber = sqrt((double)electrolytes[f]->getNTransportNumbers())*calcLimitingTransportNumber(f);
      double scaleDiffusionCoefficient = sqrt((double)electrolytes[f]->getNDiffusionCoefficients())*calcLimitingDiffusionCoefficient(f);

      for (unsigned t=0; t<electrolytes[f]->getNOsmoticCoefficients(); t++)
      {
        row = calcRowOsmoticCoefficient(f,t);
        J[row][col] = homReactions[equilibriumConstantsToFit[i]]->calcOsmoticCoefficientDerivativeEquilibriumConstant(t);
        J[row][col] /= scaleOsmoticCoefficient;
      }
      for (unsigned t=0; t<electrolytes[f]->getNActivityCoefficients(); t++)
      {
        row = calcRowActivityCoefficient(f,t);
        J[row][col] = homReactions[equilibriumConstantsToFit[i]]->calcActivityCoefficientDerivativeEquilibriumConstant(t);
        J[row][col] /= scaleActivityCoefficient;
      }
      for (unsigned t=0; t<electrolytes[f]->getNConductivities(); t++)
      {
        row = calcRowConductivity(f,t);
        J[row][col] = homReactions[equilibriumConstantsToFit[i]]->calcConductivityDerivativeEquilibriumConstant(t);
        J[row][col] /= scaleConductivity;
      }
      for (unsigned t=0; t<electrolytes[f]->getNTransportNumbers(); t++)
      {
        row = calcRowTransportNumber(f,t);
        J[row][col] = homReactions[equilibriumConstantsToFit[i]]->calcTransportNumberDerivativeEquilibriumConstant(t);
        J[row][col] /= scaleTransportNumber;
      }
      for (unsigned t=0; t<electrolytes[f]->getNDiffusionCoefficients(); t++)
      {
        row = calcRowDiffusionCoefficient(f,t);
        J[row][col] = homReactions[equilibriumConstantsToFit[i]]->calcDiffusionCoefficientDerivativeEquilibriumConstant(t);
        J[row][col] /= scaleDiffusionCoefficient;
      }
    }
    //CHECK << "J :\n";
    //for (unsigned m=0; m<nExperimentalValues; m++) {
    //  for (unsigned n=0; n<nParametersToFit; n++) {
    //    CHECK << J[m][n] << '\t';
    //    std::cout << J[m][n] << '\t';
    //  }
    //  CHECK << std::endl;
    //  std::cout << std::endl;
    //}
    //CHECK << std::endl;
    //std::cout << std::endl;
    //std::cin.get();

    // Calculate system matrix Jt*J and right hand side vector Jt*R
    //std::ofstream CHECK("Jacobian.xls");
    //CHECK.precision(16);
    for (unsigned n=0; n<nParametersToFit; n++) {
      JtR[n] = 0.;
      for (unsigned m=0; m<nExperimentalValues; m++) {
        JtR[n] -= J[m][n]*R[m];
      }
      for (unsigned p=0; p<nParametersToFit; p++) {
        JtJ[n][p] = 0.;
        for (unsigned m=0; m<nExperimentalValues; m++) {
          JtJ[n][p] += J[m][n]*J[m][p];
        }
        //std::cout << JtJ[n][p] << '\t';
        //CHECK << JtJ[n][p] << '\t';
      }
      //std::cout << '\t' << JtR[n] << std::endl;
      //CHECK << '\t' << JtR[n] << std::endl;
    }
    //std::cin.get();

    solveGauss(JtJ,JtR,nParametersToFit);

    //CHECK << "delta d :\n";
    for (unsigned i=0; i<nDiametersToFit; i++)
    {
      double d = ions[diametersToFit[i]]->getDiameter();
      double dMin = ions[diametersToFit[i]]->getMinimumDiameter();
      double dMax = ions[diametersToFit[i]]->getMaximumDiameter();
      if(d + JtR[i] > dMax) ions[diametersToFit[i]]->setDiameter(dMax);
      else if(d + JtR[i] < dMin) ions[diametersToFit[i]]->setDiameter(dMin);
      else ions[diametersToFit[i]]->setDiameter(d+JtR[i]);
      //CHECK << JtR[i] << std::endl;
      //std::cout << d+JtR[i] << std::endl;
    }
    for (unsigned i=0; i<nDiffusionConstantsToFit; i++)
    {
      double D = ions[diffusionConstantsToFit[i]]->getDiffusionConstant();
      double DMin = ions[diffusionConstantsToFit[i]]->getMinimumDiffusionConstant();
      double DMax = ions[diffusionConstantsToFit[i]]->getMaximumDiffusionConstant();
      unsigned j = i + nDiametersToFit;
      if(D + JtR[j] > DMax) ions[diffusionConstantsToFit[i]]->setDiffusionConstant(DMax);
      else if(D + JtR[j] < DMin) ions[diffusionConstantsToFit[i]]->setDiffusionConstant(DMin);
      else ions[diffusionConstantsToFit[i]]->setDiffusionConstant(D+JtR[j]);
      //std::cout << D << '\t' << JtR[j] << std::endl;
    }
    for (unsigned i=0; i<nEquilibriumConstantsToFit; i++)
    {
      double K = homReactions[equilibriumConstantsToFit[i]]->getEquilibriumConstant();
      double KMin = homReactions[equilibriumConstantsToFit[i]]->getMinimumEquilibriumConstant();
      double KMax = homReactions[equilibriumConstantsToFit[i]]->getMaximumEquilibriumConstant();
      unsigned j = i + nDiametersToFit + nDiffusionConstantsToFit;
      if(K + JtR[j] > KMax) homReactions[equilibriumConstantsToFit[i]]->setEquilibriumConstant(KMax);
      else if(K + JtR[j] < KMin) homReactions[equilibriumConstantsToFit[i]]->setEquilibriumConstant(KMin);
      else homReactions[equilibriumConstantsToFit[i]]->setEquilibriumConstant(K+JtR[j]);
      //std::cout << K << std::endl;
    }
    //std::cin.get();
    iter++;
  } while(iter<maxIter);
  std::cin.get();

  for (unsigned e=0; e<nElectrolytes; e++)
  {
    double c,sqrtc,y,yCalc,yDev;
    std::cout << '\n' << electrolytes[e]->getName();
    std::string outputFile = electrolytes[e]->getName()/*+"_"+name*/+".xls";
    std::ofstream output(outputFile.c_str());
    std::cout << "Label\tD\td" << std::endl;
    output << "Label\tD\td" << std::endl;
    for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
    {
      std::cout << electrolytes[e]->getIonLabel(i) << '\t'
        << electrolytes[e]->getIonDiffusionConstant(i) << '\t'
        << electrolytes[e]->getIonDiameter(i) << std::endl;
      output << electrolytes[e]->getIonLabel(i) << '\t'
        << electrolytes[e]->getIonDiffusionConstant(i) << '\t'
        << electrolytes[e]->getIonDiameter(i) << std::endl;
    }
    std::cout << "Label\tK" << std::endl;
    output << "Label\tK" << std::endl;
    for (unsigned r=0; r<electrolytes[e]->getNHomReactions(); r++)
    {
      std::cout << electrolytes[e]->getHomReactionLabel(r) << '\t'
        << electrolytes[e]->getHomReactionEquilibriumConstant(r) << std::endl;
      output << electrolytes[e]->getHomReactionLabel(r) << '\t'
        << electrolytes[e]->getHomReactionEquilibriumConstant(r) << std::endl;
    }

    // Osmotic coefficients
    std::cout << "\nOsmotic coefficients\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nOsmotic coefficients\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNOsmoticCoefficients(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getOsmoticCoefficientConcentration(t));
      y = electrolytes[e]->getOsmoticCoefficient(t);
      yCalc = electrolytes[e]->calcOsmoticCoefficient(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // Activity coefficients
    std::cout << "\nActivity coefficients\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nActivity coefficients\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNActivityCoefficients(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getActivityCoefficientConcentration(t));
      y = electrolytes[e]->getActivityCoefficient(t);
      yCalc = electrolytes[e]->calcActivityCoefficient(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // Conductivities
    std::cout << "\nConductivities\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nConductivities\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNConductivities(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getConductivityConcentration(t));
      y = 1000*electrolytes[e]->getConductivity(t);
      yCalc = 1000*electrolytes[e]->calcConductivity(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // ElectrolyteModel numbers
    std::cout << "\nTransport numbers\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nTransport numbers\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNTransportNumbers(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getTransportNumberConcentration(t));
      y = electrolytes[e]->getTransportNumber(t);
      yCalc = electrolytes[e]->calcTransportNumber(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // Diffusion coefficients
    std:: cout << "\nDiffusion coefficients\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nDiffusion coefficients\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNDiffusionCoefficients(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getDiffusionCoefficientConcentration(t));
      y = 1e9*electrolytes[e]->getDiffusionCoefficient(t);
      yCalc = 1e9*electrolytes[e]->calcDiffusionCoefficient(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // L++/N
    std:: cout << "\nL++\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nL++\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNConductivities(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getConductivityConcentration(t));
      y = 1e13*electrolytes[e]->getL00N(t);
      yCalc = 1e13*electrolytes[e]->calcL00N(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // L+-/N
    std:: cout << "\nL+-\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nL+-\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNConductivities(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getConductivityConcentration(t));
      y = 1e13*electrolytes[e]->getL01N(t);
      yCalc = 1e13*electrolytes[e]->calcL01N(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // L--/N
    std:: cout << "\nL--\nsqrtc\ty\tyCalc\tyDev\n";
    output << "\nL--\nsqrtc\ty\tyCalc\tyDev\n";
    for (unsigned t=0; t<electrolytes[e]->getNConductivities(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getConductivityConcentration(t));
      y = 1e13*electrolytes[e]->getL11N(t);
      yCalc = 1e13*electrolytes[e]->calcL11N(t);
      yDev = 100.*(yCalc/y-1.);
      std::cout << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      output << sqrtc << '\t' << y << '\t' << yCalc << '\t' << yDev << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        c = 0.001*electrolytes[e]->getIonConcentration(i);
        std::cout << c << '\t';
        output << c << '\t';
      }
      std::cout << '\n';
      output << '\n';
    }

    // lij/N
    std:: cout << "\nlij\nsqrtc";
    output << "\nlij\nsqrtc";
    for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
    {
      std::string iLabel = electrolytes[e]->getIonLabel(i);
      for (unsigned j=0; j<electrolytes[e]->getNIons(); j++)
      {
        std::string jLabel = electrolytes[e]->getIonLabel(j);
        std:: cout << "\tls(" << iLabel << "," << jLabel << ")";
        output << "\tls(" << iLabel << "," << jLabel << ")";
      }
    }
    std:: cout << "\n";
    output << "\n";
    for (unsigned t=0; t<electrolytes[e]->getNConductivities(); t++)
    {
      sqrtc = sqrt(0.001*electrolytes[e]->getConductivityConcentration(t));
      std::cout << sqrtc << '\t';
      output << sqrtc << '\t';
      for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
      {
        for (unsigned j=0; j<electrolytes[e]->getNIons(); j++)
        {
          yCalc = 1e13*electrolytes[e]->calcLijN(t,i,j);
          std::cout << yCalc << '\t';
          output << yCalc << '\t';
        }
      }
      std::cout << '\n';
      output << '\n';
    }
    std:: cout << std::endl;
    output << std::endl;
  }
  std::cin.get();
}
//---------------------------------------------------------------------------
unsigned Fitter::calcRowOsmoticCoefficient(unsigned e, unsigned t) const
{
  unsigned row = 0;
  for (unsigned f=0; f<e; f++)
  {
    row += electrolytes[f]->getNOsmoticCoefficients();
    row += electrolytes[f]->getNActivityCoefficients();
    row += electrolytes[f]->getNConductivities();
    row += electrolytes[f]->getNTransportNumbers();
    row += electrolytes[f]->getNDiffusionCoefficients();
  }
  row += t;
  return row;
}
//---------------------------------------------------------------------------
unsigned Fitter::calcRowActivityCoefficient(unsigned e, unsigned t) const
{
  unsigned row = calcRowOsmoticCoefficient(e,t) + electrolytes[e]->getNOsmoticCoefficients();
  return row;
}
//---------------------------------------------------------------------------
unsigned Fitter::calcRowConductivity(unsigned e, unsigned t) const
{
  unsigned row = calcRowActivityCoefficient(e,t) + electrolytes[e]->getNActivityCoefficients();
  return row;
}
//---------------------------------------------------------------------------
unsigned Fitter::calcRowTransportNumber(unsigned e, unsigned t) const
{
  unsigned row = calcRowConductivity(e,t) + electrolytes[e]->getNConductivities();
  return row;
}
//---------------------------------------------------------------------------
unsigned Fitter::calcRowDiffusionCoefficient(unsigned e, unsigned t) const
{
  unsigned row = calcRowTransportNumber(e,t) + electrolytes[e]->getNTransportNumbers();
  return row;
}
//---------------------------------------------------------------------------
double Fitter::calcLimitingOsmoticCoefficient(unsigned e) const
{
  return 1.;
}
//---------------------------------------------------------------------------
double Fitter::calcLimitingActivityCoefficient(unsigned e) const
{
  return 1.;
}
//---------------------------------------------------------------------------
double Fitter::calcLimitingConductivity(unsigned e) const
{
  int z[2];
  z[0] = electrolytes[e]->getIonChargeNumber(0);
  z[1] = electrolytes[e]->getIonChargeNumber(1);
  double D[2];
  D[0] = electrolytes[e]->getIonDiffusionConstant(0);
  D[1] = electrolytes[e]->getIonDiffusionConstant(1);
  return F_CONST*F_CONST/(R_CONST*electrolytes[e]->getSolutionTemperature())
    *(z[0]*D[0]-z[1]*D[1]);
}
//---------------------------------------------------------------------------
double Fitter::calcLimitingTransportNumber(unsigned e) const
{
  int z[2];
  z[0] = electrolytes[e]->getIonChargeNumber(0);
  z[1] = electrolytes[e]->getIonChargeNumber(1);
  double D[2];
  D[0] = electrolytes[e]->getIonDiffusionConstant(0);
  D[1] = electrolytes[e]->getIonDiffusionConstant(1);
  return z[0]*D[0]/(z[0]*D[0] - z[1]*D[1]);
}
//---------------------------------------------------------------------------
double Fitter::calcLimitingDiffusionCoefficient(unsigned e) const
{
  int z[2];
  z[0] = electrolytes[e]->getIonChargeNumber(0);
  z[1] = electrolytes[e]->getIonChargeNumber(1);
  double D[2];
  D[0] = electrolytes[e]->getIonDiffusionConstant(0);
  D[1] = electrolytes[e]->getIonDiffusionConstant(1);
  return (z[0]-z[1])*D[0]*D[1]/(z[0]*D[0] - z[1]*D[1]);
}
//---------------------------------------------------------------------------
double Fitter::calcStepLength() const
{
  double stepLength = 1.;
  double stepLengthTemp = 0.;
  for (unsigned i=0; i<nDiametersToFit; i++)
  {
    double d = ions[diametersToFit[i]]->getDiameter();
    double dMin = ions[diametersToFit[i]]->getMinimumDiameter();
    double dMax = ions[diametersToFit[i]]->getMaximumDiameter();
    if(d + JtR[i] > dMax)
    {
      stepLengthTemp = (dMax - d)/JtR[i];
    }
    else if(d + JtR[i] < dMin)
    {
      stepLengthTemp = (dMin - d)/JtR[i];
    }
    if (stepLengthTemp < stepLength) stepLength = stepLengthTemp;
  }
  for (unsigned i=0; i<nDiffusionConstantsToFit; i++)
  {
    double D = ions[diffusionConstantsToFit[i]]->getDiffusionConstant();
    double DMin = ions[diffusionConstantsToFit[i]]->getMinimumDiffusionConstant();
    double DMax = ions[diffusionConstantsToFit[i]]->getMaximumDiffusionConstant();
    unsigned j = i + nDiametersToFit;
    if(D + JtR[j] > DMax)
    {
      stepLengthTemp = (DMax - D)/JtR[j];
    }
    else if(D + JtR[j] < DMin)
    {
      stepLengthTemp = (DMin - D)/JtR[j];
    }
    if (stepLengthTemp < stepLength) stepLength = stepLengthTemp;
  }
  for (unsigned i=0; i<nEquilibriumConstantsToFit; i++)
  {
    double K = homReactions[equilibriumConstantsToFit[i]]->getEquilibriumConstant();
    double KMin = homReactions[equilibriumConstantsToFit[i]]->getMinimumEquilibriumConstant();
    double KMax = homReactions[equilibriumConstantsToFit[i]]->getMaximumEquilibriumConstant();
    unsigned j = i + nDiametersToFit + nDiffusionConstantsToFit;
    if(K + JtR[j] > KMax)
    {
      stepLengthTemp = (KMax - K)/JtR[j];
    }
    else if(K + JtR[j] < KMin)
    {
      stepLengthTemp = (KMin - K)/JtR[j];
    }
    if (stepLengthTemp < stepLength) stepLength = stepLengthTemp;
  }
  return stepLength;
}
//---------------------------------------------------------------------------

