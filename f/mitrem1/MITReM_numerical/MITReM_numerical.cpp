// MITReM_numerical_application.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <time.h>
#include <string>
//#include <fstream>
#include <iostream>
#include "DatFileReader.h"
#include "Reactor.h"

int _tmain(int argc, _TCHAR* argv[])
{
  time_t t1,t2;

  time(&t1);

  DatFileReader datFile("MITReM_numerical.dat");
  std::string reactorFile,database,ecLabel,outputFile,outputFormat,initialConditionsFile;
  bool DoE;
  datFile.readScalar("[database]",database);
  datFile.readScalar("[ecLabel]",ecLabel);
  datFile.readScalar("[reactorFile]",reactorFile);
  datFile.readScalar("[outputFile]",outputFile);
  datFile.readScalar("[outputFormat]",outputFormat);
  datFile.readScalar("[initialConditionsFile]",initialConditionsFile);
  datFile.readScalar("[DoE]",DoE);

  Reactor reactor(reactorFile,database,ecLabel,outputFile,outputFormat,DoE,initialConditionsFile);

  if (DoE == false)
  {
    reactor.solve();
  }
  else
  {
    unsigned nParameters;
    datFile.readScalar("[nParameters]",nParameters);
    std::cout << nParameters << std::endl;
    unsigned nExperiments = 2;
    for (unsigned power=1; power<nParameters; power++)
    {
      nExperiments *= 2;
    }
    std::cout << nExperiments << std::endl;
    unsigned** parameterValues = new unsigned*[nParameters];
    double** surfaceValues = new double*[nExperiments];
    unsigned pow2p = 1;
    for (unsigned p=0; p<nParameters; p++)
    {
      parameterValues[p] = new unsigned[nExperiments];
      for (unsigned e=0; e<nExperiments; e++)
      {
        parameterValues[p][e] = (e/pow2p)%2;

        //0 = lower value
        //1 = higher value

      }
      pow2p *= 2;
    }
    unsigned nIons = (nParameters+1)/4;

    for (unsigned e=0; e<nExperiments; e++)
    {
      std::cout << e << std::endl;
      char buf[10];
      std::string experimentName = "";
      for (unsigned p=0; p<nParameters; p++)
      {
        experimentName += itoa(parameterValues[p][e],buf,10);
      }
      surfaceValues[e] = new double[2*nIons+8];

      //--- INITIALIZE IONS ---//
      double c0,c1,c2,D0,D1,D2,d0,d1,d2,M0,M1,M2;
      if (nParameters == 7)
      {
        //Binary electrolyte
        int z0 = reactor.getIonChargeNumber(0);
        int z1 = reactor.getIonChargeNumber(1);

        unsigned totalConcentration = parameterValues[0][e];
        unsigned DiffConst0 = parameterValues[1][e];
        unsigned DiffConst1 = parameterValues[2][e];
        unsigned Diameter0 = parameterValues[3][e];
        unsigned Diameter1 = parameterValues[4][e];
        unsigned MolarMass0 = parameterValues[5][e];
        unsigned MolarMass1 = parameterValues[6][e];
        //std::cout << totalConcentration << '\t'
        //<< DiffConst0 << '\t'
        //<< DiffConst1 << '\t'
        //<< Diameter0 << '\t'
        //<< Diameter1 << '\t'
        //<< MolarMass0 << '\t'
        //<< MolarMass1 << std::endl;

        if (totalConcentration == 0)
        {
          c0 = -10.*z1 / (z0 - z1);
          c1 = 10.*z0 / (z0 - z1);
        }
        else if (totalConcentration == 1)
        {
          c0 = -1000.*z1 / (z0 - z1);
          c1 = 1000.*z0 / (z0 - z1);
        }
        if (abs(z0) == 1)
        {
          if (DiffConst0 == 0)
          {
            D0 = 1e-9;
          }
          else if (DiffConst0 == 1)
          {
            D0 = 2e-9;
          }
          if (Diameter0 == 0)
          {
            d0 = 2.5e-10;
          }
          else if (Diameter0 == 1)
          {
            d0 = 5.5e-10;
          }
        }
        else if (abs(z0) == 2)
        {
          if (DiffConst0 == 0)
          {
            D0 = 0.7e-9;
          }
          else if (DiffConst0 == 1)
          {
            D0 = 1.1e-9;
          }
          if (Diameter0 == 0)
          {
            d0 = 3.5e-10;
          }
          else if (Diameter0 == 1)
          {
            d0 = 6.5e-10;
          }
        }
        if (abs(z1) == 1)
        {
          if (DiffConst1 == 0)
          {
            D1 = 1e-9;
          }
          else if (DiffConst1 == 1)
          {
            D1 = 2e-9;
          }
          if (Diameter1 == 0)
          {
            d1 = 2.5e-10;
          }
          else if (Diameter1 == 1)
          {
            d1 = 5.5e-10;
          }
        }
        else if (abs(z1) == 2)
        {
          if (DiffConst1 == 0)
          {
            D1 = 0.7e-9;
          }
          else if (DiffConst1 == 1)
          {
            D1 = 1.1e-9;
          }
          if (Diameter1 == 0)
          {
            d1 = 3.5e-10;
          }
          else if (Diameter1 == 1)
          {
            d1 = 6.5e-10;
          }
        }
        if (MolarMass0 == 0)
        {
          M0 = 0.040;
        }
        else if (MolarMass0 == 1)
        {
          M0 = 0.200;
        }
        if (MolarMass1 == 0)
        {
          M1 = 0.040;
        }
        else if (MolarMass1 == 1)
        {
          M1 = 0.200;
        }
        //std::cout << c0 << '\t'
        //<< c1 << '\t'
        //<< D0 << '\t'
        //<< D1 << '\t'
        //<< d0 << '\t'
        //<< d1 << '\t'
        //<< M0 << '\t'
        //<< M1 << std::endl;
        reactor.setIonConcentration(0,c0);
        reactor.setIonConcentration(1,c1);
        reactor.setIonDiffusionConstant(0,D0);
        reactor.setIonDiffusionConstant(1,D1);
        reactor.setIonDiameter(0,d0);
        reactor.setIonDiameter(1,d1);
        reactor.setIonMolarMass(0,M0);
        reactor.setIonMolarMass(1,M1);
      }
      else if (nParameters == 11)
      {
        // Ternary electrolyte
        int z0 = reactor.getIonChargeNumber(0);
        int z1 = reactor.getIonChargeNumber(1);
        int z2 = reactor.getIonChargeNumber(2);

        unsigned totalConcentration = parameterValues[0][e];
        unsigned ratioConcentration = parameterValues[1][e];
        unsigned DiffConst0 = parameterValues[2][e];
        unsigned DiffConst1 = parameterValues[3][e];
        unsigned DiffConst2 = parameterValues[4][e];
        unsigned Diameter0 = parameterValues[5][e];
        unsigned Diameter1 = parameterValues[6][e];
        unsigned Diameter2 = parameterValues[7][e];
        unsigned MolarMass0 = parameterValues[8][e];
        unsigned MolarMass1 = parameterValues[9][e];
        unsigned MolarMass2 = parameterValues[10][e];
        //std::cout << totalConcentration << '\t'
        //<< ratioConcentration << '\t'
        //<< DiffConst0 << '\t'
        //<< DiffConst1 << '\t'
        //<< DiffConst2 << '\t'
        //<< Diameter0 << '\t'
        //<< Diameter1 << '\t'
        //<< Diameter2 << '\t'
        //<< MolarMass0 << '\t'
        //<< MolarMass1 << '\t'
        //<< MolarMass1 << std::endl;

        if (totalConcentration == 0)
        {
          if (ratioConcentration == 0)
          {
            c0 = -10.*z2 / (z0 - z2 +(z1-z2)/0.01);
            c1 = c0 / 0.01;
            c2 = -(c0*z0 + c1*z1) / z2;
          }
          else if (ratioConcentration == 1)
          {
            c0 = -10.*z2 / (z0 - z2 +(z1-z2)/10.);
            c1 = c0 / 10.;
            c2 = -(c0*z0 + c1*z1) / z2;
          }
        }
        else if (totalConcentration == 1)
        {
          if (ratioConcentration == 0)
          {
            c0 = -1000.*z2 / (z0 - z2 +(z1-z2)/0.01);
            c1 = c0 / 0.01;
            c2 = -(c0*z0 + c1*z1) / z2;
          }
          else if (ratioConcentration == 1)
          {
            c0 = -1000.*z2 / (z0 - z2 +(z1-z2)/10.);
            c1 = c0 / 10.;
            c2 = -(c0*z0 + c1*z1) / z2;
          }
        }
        if (abs(z0) == 1)
        {
          if (DiffConst0 == 0)
          {
            D0 = 1e-9;
          }
          else if (DiffConst0 == 1)
          {
            D0 = 2e-9;
          }
          if (Diameter0 == 0)
          {
            d0 = 2.5e-10;
          }
          else if (Diameter0 == 1)
          {
            d0 = 5.5e-10;
          }
        }
        else if (abs(z0) == 2)
        {
          if (DiffConst0 == 0)
          {
            D0 = 0.7e-9;
          }
          else if (DiffConst0 == 1)
          {
            D0 = 1.1e-9;
          }
          if (Diameter0 == 0)
          {
            d0 = 3.5e-10;
          }
          else if (Diameter0 == 1)
          {
            d0 = 6.5e-10;
          }
        }
        if (abs(z1) == 1)
        {
          if (DiffConst1 == 0)
          {
            D1 = 1e-9;
          }
          else if (DiffConst1 == 1)
          {
            D1 = 2e-9;
          }
          if (Diameter1 == 0)
          {
            d1 = 2.5e-10;
          }
          else if (Diameter1 == 1)
          {
            d1 = 5.5e-10;
          }
        }
        else if (abs(z1) == 2)
        {
          if (DiffConst1 == 0)
          {
            D1 = 0.7e-9;
          }
          else if (DiffConst1 == 1)
          {
            D1 = 1.1e-9;
          }
          if (Diameter1 == 0)
          {
            d1 = 3.5e-10;
          }
          else if (Diameter1 == 1)
          {
            d1 = 6.5e-10;
          }
        }
        if (abs(z2) == 1)
        {
          if (DiffConst2 == 0)
          {
            D2 = 1e-9;
          }
          else if (DiffConst2 == 1)
          {
            D2 = 2e-9;
          }
          if (Diameter2 == 0)
          {
            d2 = 2.5e-10;
          }
          else if (Diameter2 == 1)
          {
            d2 = 5.5e-10;
          }
        }
        else if (abs(z2) == 2)
        {
          if (DiffConst2 == 0)
          {
            D2 = 0.7e-9;
          }
          else if (DiffConst2 == 1)
          {
            D2 = 1.1e-9;
          }
          if (Diameter2 == 0)
          {
            d2 = 3.5e-10;
          }
          else if (Diameter2 == 1)
          {
            d2 = 6.5e-10;
          }
        }
        if (MolarMass0 == 0)
        {
          M0 = 0.040;
        }
        else if (MolarMass0 == 1)
        {
          M0 = 0.200;
        }
        if (MolarMass1 == 0)
        {
          M1 = 0.040;
        }
        else if (MolarMass1 == 1)
        {
          M1 = 0.200;
        }
        if (MolarMass2 == 0)
        {
          M2 = 0.040;
        }
        else if (MolarMass2 == 1)
        {
          M2 = 0.200;
        }
        //std::cout << c0 << '\t'
        //<< c1 << '\t'
        //<< c2 << '\t'
        //<< D0 << '\t'
        //<< D1 << '\t'
        //<< D2 << '\t'
        //<< d0 << '\t'
        //<< d1 << '\t'
        //<< d2 << '\t'
        //<< M0 << '\t'
        //<< M1 << '\t'
        //<< M2 << std::endl;
        reactor.setIonConcentration(0,c0);
        reactor.setIonConcentration(1,c1);
        reactor.setIonConcentration(2,c2);
        reactor.setIonDiffusionConstant(0,D0);
        reactor.setIonDiffusionConstant(1,D1);
        reactor.setIonDiffusionConstant(2,D2);
        reactor.setIonDiameter(0,d0);
        reactor.setIonDiameter(1,d1);
        reactor.setIonDiameter(2,d2);
        reactor.setIonMolarMass(0,M0);
        reactor.setIonMolarMass(1,M1);
        reactor.setIonMolarMass(2,M2);
      }
      else
      {
        std::cout << "Not allowed" << std::endl;
        system("pause");
        exit(1);
      }

      reactor.setElectrolyteModel("MSA");
      reactor.setConductivity(0.);
      double conductivity = reactor.getConductivity();
      //std::cout << "conductivity = " << conductivity << std::endl;
      reactor.init();
      reactor.solve();
      for (unsigned i=0; i<nIons; i++)
      {
        surfaceValues[e][2*i] = reactor.xVec[i]; //surface concentrations
      }
      surfaceValues[e][2*nIons] = reactor.xVec[nIons]; //surface potential
      surfaceValues[e][2*nIons+2] = reactor.ionCurrentDensities[0][nIons][0];
      surfaceValues[e][2*nIons+4] = reactor.elecReactionCurrentDensities[0][0][0][0][0];
      surfaceValues[e][2*nIons+6] = reactor.residu;

      if (nParameters == 7)
      {
        //Binary electrolyte
        reactor.setIonConcentration(0,c0);
        reactor.setIonConcentration(1,c1);
        reactor.setIonDiffusionConstant(0,D0);
        reactor.setIonDiffusionConstant(1,D1);
        reactor.setIonDiameter(0,d0);
        reactor.setIonDiameter(1,d1);
        reactor.setIonMolarMass(0,M0);
        reactor.setIonMolarMass(1,M1);
      }
      else if (nParameters == 11)
      {
        // Ternary electrolyte
        reactor.setIonConcentration(0,c0);
        reactor.setIonConcentration(1,c1);
        reactor.setIonConcentration(2,c2);
        reactor.setIonDiffusionConstant(0,D0);
        reactor.setIonDiffusionConstant(1,D1);
        reactor.setIonDiffusionConstant(2,D2);
        reactor.setIonDiameter(0,d0);
        reactor.setIonDiameter(1,d1);
        reactor.setIonDiameter(2,d2);
        reactor.setIonMolarMass(0,M0);
        reactor.setIonMolarMass(1,M1);
        reactor.setIonMolarMass(2,M2);
      }
      //reactor.outputFile = experimentName + "_Ideal";
      reactor.setElectrolyteModel("Ideal");
      reactor.setConductivity(conductivity);
      reactor.init();
      reactor.solve();
      for (unsigned i=0; i<nIons; i++)
      {
        surfaceValues[e][2*i+1] = reactor.xVec[i]; //surface concentrations
      }
      surfaceValues[e][2*nIons+1] = reactor.xVec[nIons]; //surface potential
      surfaceValues[e][2*nIons+3] = reactor.ionCurrentDensities[0][nIons][0];
      surfaceValues[e][2*nIons+5] = reactor.elecReactionCurrentDensities[0][0][0][0][0];
      surfaceValues[e][2*nIons+7] = reactor.residu;

      //for (unsigned q=0; q<2*nIons+8; q++)
      //{
      //  std::cout << surfaceValues[e][q] << std::endl;
      //}
      //system("pause");
      //break;
    }

    outputFile += ".xls";
    std::ofstream output;
    output.setf(std::ios::scientific);
    output.precision(12);
    output.open(outputFile.c_str(), std::ofstream::out | std::ofstream::app);
    output << "Experiment\t";
    for (unsigned i=0; i<nIons; i++)
    {
      output << "c(" << i << ")<Ideal/MSA-1> %\t";
    }
    output << "Jt<Ideal/MSA-1> %\tJr<Ideal/MSA-1> %\t";
    for (unsigned i=0; i<nIons; i++)
    {
      output << "c(" << i << ")<MSA>\tc(" << i << ")<Ideal>\t";
    }
    output << "U<MSA>\tU<Ideal>\tJt<MSA>\tJt<Ideal>\tJr<MSA>\tJr<Ideal>\tResidu<MSA>\tResidu<Ideal>" << std::endl;
    for (unsigned e=0; e<nExperiments; e++)
    {
      output << e << '\t';
      for (unsigned i=0; i<nIons; i++)
      {
        output << 100.*(surfaceValues[e][2*i+1]/surfaceValues[e][2*i]-1.) << '\t';
      }
      output << 100.*(surfaceValues[e][2*nIons+3]/surfaceValues[e][2*nIons+2]-1.) << '\t'
           << 100.*(surfaceValues[e][2*nIons+5]/surfaceValues[e][2*nIons+4]-1.) << '\t';
      for (unsigned i=0; i<nIons; i++)
      {
        output << surfaceValues[e][2*i] << '\t' << surfaceValues[e][2*i+1] << '\t';
      }
      output << surfaceValues[e][2*nIons] << '\t' << surfaceValues[e][2*nIons+1] << '\t';
      output << surfaceValues[e][2*nIons+2] << '\t' << surfaceValues[e][2*nIons+3] << '\t';
      output << surfaceValues[e][2*nIons+4] << '\t' << surfaceValues[e][2*nIons+5] << '\t';
      output << surfaceValues[e][2*nIons+6] << '\t' << surfaceValues[e][2*nIons+7] << std::endl;
      //break;
    }
    output.close();
  }

  time(&t2);

  std::cout << "Total Time for the calculations = " <<  difftime (t2,t1) << " seconds" << std::endl;
  char ch;
  std::cin >> ch;

  return 0;
}

