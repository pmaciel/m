//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------

#include "Reactor.h"
#include <math.h>
#include "DatFileReader.h"
#include "TypeDefs.h"


#include <iostream>
#include <sstream>
#include <string>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Reactor::Reactor (const std::string &name, const std::string &database, const std::string &ecLabel, const std::string &outputFile, const std::string &outputFormat, bool DoE, const std::string &initialConditionsFile)
{
  this->DoE = DoE;
  this->initialConditionsFile = initialConditionsFile;

  mitrem = new MITReM(database, ecLabel);

  nIons = mitrem->getNIons();
  nVariables = nIons+1;

  reactorFile = name + ".reactor";
  nodesFile = name + ".nodes";
  elementsFile = name + ".elements";
  boundaryElementsFile = name + ".boundaryelements";
  boundariesFile = name + ".boundaries";
  electrodesFile = name + ".electrodes";
  flowFieldFile = name + ".flowfield";
  magneticFieldFile = name + ".magneticfield";
  solverFile = name + ".solver";
  electrodePotentialsFile = name + ".electrodepotentials";
  miotrasDatFile = name + ".dat";
  miotrasFlowFile = name + ".flow";
  this->outputFile = outputFile;
  this->outputFormat = outputFormat;

  std::cout << "Read reactor" << std::endl;
  readReactor();
  std::cout << "Read nodes" << std::endl;
  readNodes();
  std::cout << "Read elements" << std::endl;
  readElements();
  std::cout << "Read boundary elements" << std::endl;
  readBoundaryElements();
  std::cout << "Read boundaries" << std::endl;
  readBoundaries();
  std::cout << "Read electrodes" << std::endl;
  readElectrodes();
  std::cout << "Read flow field" << std::endl;
  readFlowField();
  std::cout << "Read magnetic field" << std::endl;
  readMagneticField();
  std::cout << "Read solver" << std::endl;
  readSolver();
  std::cout << "Read electrode potentials" << std::endl;
  readElectrodePotentials();

  coordinates = new DoubleVector[nElementNodes]; // 1 too large for boundary elements...
  velocities = new DoubleVector[nElementNodes]; // 1 too large for boundary elements...
  magneticFieldVectors = new DoubleVector[nElementNodes]; // 1 too large for boundary elements...
  concentrations = new DoubleVector[nElementNodes]; // 1 too large for boundary elements...
  potentials = new double[nElementNodes];
  temperatures = new double[nElementNodes];
  densities = new double[nElementNodes];
  voidFractions = new double[nElementNodes];
  double T = mitrem->getSolutionTemperature();
  double rho = mitrem->getSolutionDensity();
  double alpha = 0.;
  for (unsigned m=0; m<nElementNodes; m++)
  {
    temperatures[m] = T;
    densities[m] = rho;
    voidFractions[m] = alpha;
  }

  nIons = mitrem->getNIons();
  nVariables = nIons+1;
  inletVec = new double[nVariables];
  for (unsigned i=0; i<nIons; i++)
  {
    inletVec[i] = mitrem->getIonInletConcentration(i);
  }
  inletVec[nIons] = 0; // to be adjusted: 1) average of all electrode potentials, or 2) solve set per electrode and take average

  size = nNodes*nVariables;
  xVec = new double[size];
  xVecOld = new double[size];
  bVec = new double[size];
  bVecOld = new double[size];

  ionCurrentDensities = new double**[nElements];
  for (unsigned e=0; e<nElements; e++)
  {
    ionCurrentDensities[e] = new double*[nVariables];
    for (unsigned i=0; i<nVariables; i++)
    {
      ionCurrentDensities[e][i] = new double[nDimensions];
    }
  }

  gasReactionRates = new double****[nElectrodes];
  //std::cout << "nElectrodes = " << nElectrodes << std::endl;
  elecReactionCurrentDensities = new double****[nElectrodes];
  unsigned nGasReactions = mitrem->getNGasReactions();
  unsigned nElecReactions = mitrem->getNElecReactions();
  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    gasReactionRates[e] = new double***[nBoundariesOfElectrode];
    //std::cout << "nBoundariesOfElectrode = " << nBoundariesOfElectrode << std::endl;
    elecReactionCurrentDensities[e] = new double***[nBoundariesOfElectrode];
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      gasReactionRates[e][b] = new double**[nBoundaryElementsOfBoundary];
      //std::cout << "nBoundaryElementsOfBoundary = " << nBoundaryElementsOfBoundary << std::endl;
      elecReactionCurrentDensities[e][b] = new double**[nBoundaryElementsOfBoundary];
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        gasReactionRates[e][b][be] = new double*[nBoundaryElementNodes];
        //std::cout << "nBoundaryElementNodes = " << nBoundaryElementNodes << std::endl;
        elecReactionCurrentDensities[e][b][be] = new double*[nBoundaryElementNodes];
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          gasReactionRates[e][b][be][m] = new double[nGasReactions];
          //std::cout << "nElecReactions = " << nElecReactions << std::endl;
          elecReactionCurrentDensities[e][b][be][m] = new double[nElecReactions];
          //system("pause");
        }
      }
    }
  }
  //system("pause");
  currents = new double[nElectrodes];

  /*output.open("Arguments.txt", std::ofstream::out);
  output << "node\tx\ty\tc(A+)\tc(B-)\tU\n";
  for (unsigned m=0; m<nBoundaryElementNodes; m++)
  {
    output << m << '\t';
    coordinates[m] = nodes[m]->getComponents();
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << coordinates[m][c] << '\t';
    }
    for (unsigned i=0; i<nIons; i++)
    {
      xVec[var(m,i)] = mitrem->getIonConcentration(i);
    }
    concentrations[m] = &xVec[var(m,0)];
    for (unsigned i=0; i<nIons; i++)
    {
      output << concentrations[m][i] << '\t';
    }
    potentials[m] = 0.;
    output << potentials[m] << '\n';
  }
  output << std::endl;
  //unsigned nElecReactions = electrode->getNElecReactions();
  Electrode* electrode = electrodes[0];
  IndexList elecReactions = electrode->getElecReactions();
  double electrodePotential = 0.5;
  output << "V = " << electrodePotential << std::endl;
  output.close();
  output.open("boundaryElementVec.txt", std::ofstream::out);
  DoubleVector boundaryElementVec = elementMatrixAssembler->calcBoundaryElementVec(coordinates, concentrations, potentials, temperatures, densities, voidFractions, elecReactions, nElecReactions, electrodePotential);
  for (unsigned m=0; m<nBoundaryElementNodes; m++)
  {
    for (unsigned i=0; i<nVariables; i++)
    {
      unsigned rowLocal = var(m,i);
      output << boundaryElementVec[rowLocal] << std::endl;
    }
  }
  output.close();
  output.open("boundaryElementJac.txt", std::ofstream::out);
  DoubleMatrix boundaryElementJac = elementMatrixAssembler->calcBoundaryElementJac(coordinates, concentrations, potentials, temperatures, densities, voidFractions, elecReactions, nElecReactions, electrodePotential);
  for (unsigned m=0; m<nBoundaryElementNodes; m++)
  {
    for (unsigned i=0; i<nVariables; i++)
    {
      unsigned rowLocal = var(m,i);
      for (unsigned n=0; n<nBoundaryElementNodes; n++)
      {
        for (unsigned j=0; j<nVariables; j++)
        {
          unsigned colLocal = var(n,j);
          output << boundaryElementJac[rowLocal][colLocal] << '\t';
        }
      }
      output << std::endl;
    }
  }
  output.close();
  exit(1);*/
}
//---------------------------------------------------------------------------
Reactor::~Reactor ()
{
  for (unsigned m=0; m<nNodes; m++)
  {
    delete nodes[m];
    delete flowField[m];
    delete magneticField[m];
  }
  for (unsigned e=0; e<nElements; e++)
  {
    delete elements[e];
  }
  for (unsigned b=0; b<nBoundaryElements; b++)
  {
    delete boundaryElements[b];
  }
  for (unsigned b=0; b<nBoundaries; b++)
  {
    delete boundaries[b];
  }
  for (unsigned e=0; e<nElectrodes; e++)
  {
    delete electrodes[e];
  }
  delete[] nodes;
  delete[] elements;
  delete[] boundaryElements;
  delete[] boundaries;
  delete[] electrodes;
  delete[] flowField;
  delete[] magneticField;

  delete[] coordinates;
  delete[] velocities;
  delete[] magneticFieldVectors;
  delete[] concentrations;
  delete[] potentials;

  delete[] inletVec;
  delete[] xVec;
  delete[] xVecOld;
  delete[] bVec;
  delete[] bVecOld;


  for (unsigned e=0; e<nElements; e++)
  {
    for (unsigned i=0; i<nVariables; i++)
    {
      delete[] ionCurrentDensities[e][i];
    }
    delete[] ionCurrentDensities[e];
  }
  delete[] ionCurrentDensities;

  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          delete[] gasReactionRates[e][b][be][m];
          delete[] elecReactionCurrentDensities[e][b][be][m];
        }
        delete[] gasReactionRates[e][b][be];
        delete[] elecReactionCurrentDensities[e][b][be];
      }
      delete[] gasReactionRates[e][b];
      delete[] elecReactionCurrentDensities[e][b];
    }
    delete[] gasReactionRates[e];
    delete[] elecReactionCurrentDensities[e];
  }
  delete[] gasReactionRates;
  delete[] elecReactionCurrentDensities;

  delete[] currents;

}
//---------------------------------------------------------------------------
void Reactor::readReactor ()
{
  DatFileReader datFile(reactorFile);

  //datFile.readScalar("[nDimensions]",nDimensions);
  datFile.readScalar("[dimensions]",dimensions);

  if (dimensions == "1D")
  {
    nDimensions = 1;
  }

  else if (dimensions == "2D")
  {
    nDimensions = 2;
  }

  else if (dimensions == "AX")
  {
    nDimensions = 2;
  }

  else if (dimensions == "3D")
  {
    nDimensions = 3;
  }

  else if (dimensions == "RDE_1D")
  {
    dimensions = "1D";
    nDimensions = 1;

    double length,refFac,vRot;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[vRot]",vRot);

    // Write elements file
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e++)
    {
      elementsOutput << e << '\t' << e+1 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+1;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    vRot *= M_PI/30.; // convert rpm to rad/s
    double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
    double Quot = vRot/kinematicViscosity;
    double sqrtQuot = sqrt(Quot);
    double sqrtProd = sqrt(vRot*kinematicViscosity);
    if (refFac > 1.)
    {
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << std::endl;
        flowFieldOutput << v << std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      for (unsigned m=0; m<nNodes; m++)
      {
        double x = length*m/nElements;
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << std::endl;
        flowFieldOutput << v << std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 0 << std::endl;
    boundaryElementsOutput << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Electrode\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 1;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
    electrodesOutput << "[nBoundaries0] = 1" << std::endl;
    electrodesOutput << "\t<boundary> = 0\n" << std::endl;
    unsigned nElecReactions = mitrem->getNElecReactions();
    electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
    for (unsigned r=0; r<nElecReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
    }
    unsigned nGasReactions = mitrem->getNGasReactions();
    electrodesOutput << "\n[nGasReactions0] = " << nGasReactions << std::endl;
    for (unsigned r=0; r<nGasReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getGasReactionLabel(r) << std::endl;
    }
    electrodesOutput.close();
  }

  else if (dimensions == "RDE_2D")
  {
    dimensions = "2D";
    nDimensions = 2;

    double length,refFac,vRot;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[vRot]",vRot);

    // Write elements file
    nElements *= 2;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=2)
    {
      elementsOutput << e << '\t' << e+2 << '\t' << e+1 << std::endl;
      elementsOutput << e+1 << '\t' << e+2 << '\t' << e+3 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+2;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    vRot *= M_PI/30.; // convert rpm to rad/s
    double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
    double Quot = vRot/kinematicViscosity;
    double sqrtQuot = sqrt(Quot);
    double sqrtProd = sqrt(vRot*kinematicViscosity);
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/2; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 2.*length/nElements;
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = 2.*length*m/nElements;
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 1 << '\t' << 0 << std::endl;
    boundaryElementsOutput << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Electrode\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 1;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
    electrodesOutput << "[nBoundaries0] = 1" << std::endl;
    electrodesOutput << "\t<boundary> = 0\n" << std::endl;
    unsigned nElecReactions = mitrem->getNElecReactions();
    electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
    for (unsigned r=0; r<nElecReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
    }
    unsigned nGasReactions = mitrem->getNGasReactions();
    electrodesOutput << "\n[nGasReactions0] = " << nGasReactions << std::endl;
    for (unsigned r=0; r<nGasReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getGasReactionLabel(r) << std::endl;
    }
    electrodesOutput.close();
  }

  //else if (dimensions == "RDE_2Dbis")
  //{
  //  dimensions = "2D";
  //  nDimensions = 2;
  //
  //  double length,refFac,vRot;
  //
  //  datFile.readScalar("[nElements]",nElements);
  //  datFile.readScalar("[length]",length);
  //  datFile.readScalar("[refFac]",refFac);
  //  datFile.readScalar("[vRot]",vRot);

  //  // Write elements file
  //  nElements *= 2;
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  for (unsigned e=0; e<nElements; e+=4)
  //  {
  //    unsigned e_ = e*3/4;
  //    elementsOutput << e_ << '\t' << e_+2 << '\t' << e_+1 << std::endl;
  //    elementsOutput << e_ << '\t' << e_+3 << '\t' << e_+2 << std::endl;
  //    elementsOutput << e_+1 << '\t' << e_+2 << '\t' << e_+4 << std::endl;
  //    elementsOutput << e_+2 << '\t' << e_+3 << '\t' << e_+4 << std::endl;
  //  }
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = nElements*3/4+2;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  vRot *= M_PI/30.; // convert rpm to rad/s
  //  double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
  //  double Quot = vRot/kinematicViscosity;
  //  double sqrtQuot = sqrt(Quot);
  //  double sqrtProd = sqrt(vRot*kinematicViscosity);
  //  if (refFac > 1.)
  //  {
  //    double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
  //    double refFacN = refFac;
  //    double refFacm = 1.;
  //    for (unsigned p=1; p<nElements/2; p++)
  //    {
  //      refFacN *= refFac;
  //    }
  //    for (unsigned m=0; m<nElements/2+1; m++)
  //    {
  //      double x = length*(1-refFacm)/(1-refFacN);
  //      double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
  //      if (m%2 == 0)
  //      {
  //        nodesOutput << x << '\t' << -H/2 << std::endl;
  //        nodesOutput << x << '\t' << H/2 << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
  //      }
  //      else
  //      {
  //        nodesOutput << x << '\t' << 0.0 << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
  //      }
  //      refFacm *= refFac;
  //    }
  //  }
  //  else if (refFac == 1.)
  //  {
  //    double H = 2.*length/nElements;
  //    for (unsigned m=0; m<nElements/2+1; m++)
  //    {
  //      double x = 2.*length*m/nElements;
  //      double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
  //      if (m%2 == 0)
  //      {
  //        nodesOutput << x << '\t' << -H/2 << std::endl;
  //        nodesOutput << x << '\t' << H/2 << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
  //      }
  //      else
  //      {
  //        nodesOutput << x << '\t' << 0.0 << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
  //      }
  //    }
  //  }
  //  else
  //  {
  //    // error message
  //  }
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 2;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 1 << '\t' << 0 << std::endl;
  //  boundaryElementsOutput << nNodes-2 << '\t' << nNodes-1 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 2;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<dimensions> = Electrode\n" << std::endl;
  //  boundariesOutput << "\t<dimensions> = VirtualElectrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 1;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}

  else if (dimensions == "RDE_AX")
  {
    dimensions = "AX";
    nDimensions = 2;

    double length,refFac,vRot;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[vRot]",vRot);

    // Write elements file
    nElements *= 2;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=2)
    {
      elementsOutput << e << '\t' << e+2 << '\t' << e+1 << std::endl;
      elementsOutput << e+1 << '\t' << e+2 << '\t' << e+3 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+2;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    vRot *= M_PI/30.; // convert rpm to rad/s
    double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
    double Quot = vRot/kinematicViscosity;
    double sqrtQuot = sqrt(Quot);
    double sqrtProd = sqrt(vRot*kinematicViscosity);
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/2; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 2.*length/nElements;
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = 2.*length*m/nElements;
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 1 << '\t' << 0 << std::endl;
    boundaryElementsOutput << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Electrode\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 1;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
    electrodesOutput << "[nBoundaries0] = 1" << std::endl;
    electrodesOutput << "\t<boundary> = 0\n" << std::endl;
    unsigned nElecReactions = mitrem->getNElecReactions();
    electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
    for (unsigned r=0; r<nElecReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
    }
    unsigned nGasReactions = mitrem->getNGasReactions();
    electrodesOutput << "\n[nGasReactions0] = " << nGasReactions << std::endl;
    for (unsigned r=0; r<nGasReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getGasReactionLabel(r) << std::endl;
    }
    electrodesOutput.close();
  }

  else if (dimensions == "RDE_3D")
  {
    dimensions = "3D";
    nDimensions = 3;

    double length,refFac,vRot;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[vRot]",vRot);

    // Write elements file
    nElements *= 3;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=3)
    {
      elementsOutput << e << '\t' << e+1 << '\t' << e+2 << '\t' << e+3 << std::endl;
      elementsOutput << e+1 << '\t' << e+3 << '\t' << e+4 << '\t' << e+5 << std::endl;
      elementsOutput << e+1 << '\t' << e+2 << '\t' << e+3 << '\t' << e+5 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+3;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    vRot *= M_PI/30.; // convert rpm to rad/s
    double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
    double Quot = vRot/kinematicViscosity;
    double sqrtQuot = sqrt(Quot);
    double sqrtProd = sqrt(vRot*kinematicViscosity);
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/3; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/3; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 0 << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 3.*length/nElements;
      for (unsigned m=0; m<nNodes/3; m++)
      {
        double x = 3.*length*m/nElements;
        double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
        nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 0 << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 0 << '\t' << 2 << '\t' << 1 << std::endl;
    boundaryElementsOutput << nNodes-3 << '\t' << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Electrode\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 1;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
    electrodesOutput << "[nBoundaries0] = 1" << std::endl;
    electrodesOutput << "\t<boundary> = 0\n" << std::endl;
    unsigned nElecReactions = mitrem->getNElecReactions();
    electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
    for (unsigned r=0; r<nElecReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
    }
    unsigned nGasReactions = mitrem->getNGasReactions();
    electrodesOutput << "\n[nGasReactions0] = " << nGasReactions << std::endl;
    for (unsigned r=0; r<nGasReactions; r++)
    {
      electrodesOutput << "\t<label> = " << mitrem->getGasReactionLabel(r) << std::endl;
    }
    electrodesOutput.close();
  }

  else if ((dimensions == "MIoTraS_2D") || (dimensions == "MIoTraS_AX"))
  {
    if (dimensions == "MIoTraS_2D") dimensions = "2D";
    else if (dimensions == "MIoTraS_AX") dimensions = "AX";
    nDimensions = 2;

    std::cout << "Reading Miotras dat file..." << std::endl;
    input.open(miotrasDatFile.c_str());
    if (input == NULL) std::cout << "Error: Miotras dat file does not exist" << std::endl;
    input.setf(std::ios::scientific);
    char temp[256];

    skip(input,5,'\n');

    // Read number of parameters
    input.getline(temp,256,'\n');
    unsigned nParameters = atoi(temp);

    skip(input,nParameters,'\n');


    input.getline(temp,256,'\n');
    const unsigned n = atoi(temp);
    skip(input,n,'\n');

    // Read number of boundaries
    input.getline(temp,256,'\n');
    nBoundaries = atoi(temp);
    std::ofstream boundariesOutput(boundariesFile.c_str());
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;

    // Read boundary types
    for (unsigned b=0; b<nBoundaries; b++)
    {
      skip(input,8,'\t');
      input.getline(temp,256,'\t');
      unsigned boundaryTypeIndex = atoi(temp);
      switch (boundaryTypeIndex)
      {
        case 1 :
          boundariesOutput << "\t<type> = " << "Electrode" << std::endl;
          break;
        case 2 :
          boundariesOutput << "\t<type> = " << "VirtualElectrode" << std::endl;
          break;
        case 3 :
          boundariesOutput << "\t<type> = " << "Inlet" << std::endl;
          break;
        case 4 :
          boundariesOutput << "\t<type> = " << "Outlet" << std::endl;
          break;
        default :
          boundariesOutput << "\t<type> = " << "Insulator" << std::endl;
      }
      skip(input,1,'\n');
    }
    boundariesOutput << std::endl;
    input.close();

    std::cout << "Reading Miotras flow file..." << std::endl;

    input.open(miotrasFlowFile.c_str());
    std::cout << miotrasFlowFile.c_str() << std::endl;
    if (input == NULL) std::cout << "Error: Miotras flow file does not exist" << std::endl;
    input.setf(std::ios::scientific);
    input.getline(temp,256,'\n');

    std::cout << "Writing nodes and flowfield file..." << std::endl;

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());

    input.getline(temp,256,'\n');
    nNodes = atoi(temp);
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    for (unsigned m=0; m<nNodes; m++)
    {
      skip(input,1,'\t');
      input.getline(temp,256,'\t');
      nodesOutput << temp << '\t';
      input.getline(temp,256,'\t');
      nodesOutput << temp << std::endl;
      input.getline(temp,256,'\t');
      flowFieldOutput << temp << '\t';
      input.getline(temp,256,'\t');
      flowFieldOutput << temp << std::endl;
      input.getline(temp,256,'\n');
    }
    nodesOutput.close();
    flowFieldOutput.close();

    std::cout << "Writing elements and boundaryelements file..." << std::endl;

    // Write elements and boundary elements file
    std::ofstream elementsOutput(elementsFile.c_str());
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());

    skip(input,1,' ');
    input.getline(temp,256,' ');
    nElements = atoi(temp);
    input.getline(temp,256,'\n');
    nBoundaryElements = atoi(temp);
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    for (unsigned e=0; e<nElements; e++)
    {
      skip(input,1,' ');
      skip(input,1,'\t');
      for (unsigned m=0; m<3; m++)
      {
        input.getline(temp,256,' ');
        unsigned index = atoi(temp) - 1;
        elementsOutput << index << "\t";
      }
      elementsOutput << std::endl;
      input.getline(temp,256,'\n');
    }
    elementsOutput.close();
    std::vector< std::vector<unsigned> > boundaryElementsTemp;
    boundaryElementsTemp.resize(nBoundaries);
    unsigned boundaryIndex;
    for (unsigned e=0; e<nBoundaryElements; e++)
    {
      input.getline(temp,256,' ');
      boundaryIndex = atoi(temp) - 1;
      boundaryElementsTemp[boundaryIndex].push_back(e);
      skip(input,1,'\t');
      for (unsigned m=0; m<2; m++)
      {
        input.getline(temp,256,' ');
        unsigned index = atoi(temp) - 1;
        boundaryElementsOutput << index << "\t";
      }
      boundaryElementsOutput << std::endl;
      input.getline(temp,256,'\n');
    }
    boundaryElementsOutput.close();

    std::cout << "Writing boundaries file..." << std::endl;

    // Write boundaries file
    unsigned nBoundaryElementsTemp;
    for (unsigned b=0; b<nBoundaries; b++)
    {
      nBoundaryElementsTemp = boundaryElementsTemp[b].size();
      boundariesOutput << "[nBoundaryElements" << b << "] = " << nBoundaryElementsTemp << std::endl;
      for (unsigned e=0; e<nBoundaryElementsTemp; e++)
      {
        boundariesOutput << "\t<boundaryElement> = " << boundaryElementsTemp[b][e] << std::endl;
      }
      boundariesOutput << std::endl;
    }
    boundariesOutput.close();
  }
  //else if (dimensions == "RDE_3Dbis")
  //{
  //  dimensions = "3D";
  //  nDimensions = 3;
  //
  //  double length,refFac,vRot;
  //
  //  datFile.readScalar("[nElements]",nElements);
  //  datFile.readScalar("[length]",length);
  //  datFile.readScalar("[refFac]",refFac);
  //  datFile.readScalar("[vRot]",vRot);

  //  // Write elements file
  //  nElements *= 4;
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  for (unsigned e=0; e<nElements; e+=8)
  //  {
  //    unsigned e_ = e/2;
  //    elementsOutput << e_ << '\t' << e_+2 << '\t' << e_+1 << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_ << '\t' << e_+1 << '\t' << e_+6 << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_+1 << '\t' << e_+2 << '\t' << e_+4 << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_+2 << '\t' << e_ << '\t' << e_+5 << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_+5 << '\t' << e_+4 << '\t' << e_+2 << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_+4 << '\t' << e_+6 << '\t' << e_+1 << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_+6 << '\t' << e_+5 << '\t' << e_ << '\t' << e_+3 << std::endl;
  //    elementsOutput << e_+4 << '\t' << e_+5 << '\t' << e_+6 << '\t' << e_+3 << std::endl;
  //  }
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = nElements/2+3;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  vRot *= M_PI/30.; // convert rpm to rad/s
  //  double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
  //  double Quot = vRot/kinematicViscosity;
  //  double sqrtQuot = sqrt(Quot);
  //  double sqrtProd = sqrt(vRot*kinematicViscosity);
  //  if (refFac > 1.)
  //  {
  //    double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
  //    double refFacN = refFac;
  //    double refFacm = 1.;
  //    for (unsigned p=1; p<nElements/4; p++)
  //    {
  //      refFacN *= refFac;
  //    }
  //    for (unsigned m=0; m<nElements/4+1; m++)
  //    {
  //      double x = length*(1-refFacm)/(1-refFacN);
  //      double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
  //      if (m%4 == 0)
  //      {
  //        nodesOutput << x << '\t' << 0.0 << '\t' << -H/sqrt(3.) << std::endl;
  //        nodesOutput << x << '\t' << H/2. << '\t' << 0.5*H/sqrt(3.) << std::endl;
  //        nodesOutput << x << '\t' << -H/2. << '\t' << 0.5*H/sqrt(3.) << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      }
  //      else if (m%4 == 2)
  //      {
  //        nodesOutput << x << '\t' << 0.0 << '\t' << H/sqrt(3.) << std::endl;
  //        nodesOutput << x << '\t' << -H/2. << '\t' << -0.5*H/sqrt(3.) << std::endl;
  //        nodesOutput << x << '\t' << H/2. << '\t' << -0.5*H/sqrt(3.) << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      }
  //      else
  //      {
  //        nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 << std::endl;
  //      }
  //      refFacm *= refFac;
  //    }
  //  }
  //  else if (refFac == 1.)
  //  {
  //    double H = 4.*length/nElements;
  //    for (unsigned m=0; m<nElements/4+1; m++)
  //    {
  //      double x = 4.*length*m/nElements;
  //      double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
  //      if (m%2 == 0)
  //      {
  //        nodesOutput << x << '\t' << 0.0 << '\t' << -H/sqrt(3.) << std::endl;
  //        nodesOutput << x << '\t' << H/2. << '\t' << 0.5*H/sqrt(3.) << std::endl;
  //        nodesOutput << x << '\t' << -H/2. << '\t' << 0.5*H/sqrt(3.) << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      }
  //      else
  //      {
  //        nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
  //        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 << std::endl;
  //      }
  //    }
  //  }
  //  else
  //  {
  //    // error message
  //  }
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 2;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 0 << '\t' << 2 << '\t' << 1 << std::endl;
  //  boundaryElementsOutput << nNodes-3 << '\t' << nNodes-2 << '\t' << nNodes-1 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 2;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<dimensions> = Electrode\n" << std::endl;
  //  boundariesOutput << "\t<dimensions> = VirtualElectrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 1;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}
  //else if (dimensions == "Debug_1D")
  //{
  //  nDimensions = 1;

  //  // Write elements file
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  nElements = 1;
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  elementsOutput << 0 << '\t' << 1 << std::endl;
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = 2;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  nodesOutput << 0.0 << std::endl;
  //  nodesOutput << 1.0 << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << 1.0 << std::endl;
  //  flowFieldOutput << 1.0 << std::endl;
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 1;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 0 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 1;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<dimensions> = Electrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 1;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}
  //else if (dimensions == "Debug_2D")
  //{
  //  nDimensions = 2;
  //
  //  // Write elements file
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  nElements = 1;
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  elementsOutput << 0 << '\t' << 1 << '\t' << 2 << std::endl;
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = 3;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  nodesOutput << 0.0 << '\t' << 0.0 << std::endl;
  //  nodesOutput << 1.0 << '\t' << 0.0 << std::endl;
  //  nodesOutput << 1./sqrt(3.) << '\t' << 1.0 << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << std::endl;
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 1;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 0 << '\t' << 1 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 1;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<dimensions> = Electrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 1;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}
  //else if (dimensions == "Debug_3D")
  //{
  //  nDimensions = 3;
  //
  //  // Write elements file
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  nElements = 1;
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  elementsOutput << 0 << '\t' << 1 << '\t' << 2 << '\t' << 3 << std::endl;
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = 4;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  nodesOutput << 0.0 << '\t' << 0.0 << '\t' << 0.0 << std::endl;
  //  nodesOutput << 1.0 << '\t' << 0.0 << '\t' << 0.0 << std::endl;
  //  nodesOutput << 1./sqrt(3.) << '\t' << 1.0 << '\t' << 0.0 << std::endl;
  //  nodesOutput << 0.5 << '\t' << 0.5 << '\t' << -1.0 << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << '\t' << 0.0 << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << '\t' << 0.0 << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << '\t' << 0.0 << std::endl;
  //  flowFieldOutput << cos(M_PI/6.) << '\t' << sin(M_PI/6.) << '\t' << 0.0 << std::endl;
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 1;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 0 << '\t' << 1 << '\t' << 2 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 1;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<dimensions> = Electrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 1;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}
  //else if (dimensions == "RDE_3Dbis")
  //{
  //  nDimensions = 3;
  //
  //  double length,refFac,vRot;
  //
  //  datFile.readScalar("[nElements]",nElements);
  //  datFile.readScalar("[length]",length);
  //  datFile.readScalar("[refFac]",refFac);
  //  datFile.readScalar("[vRot]",vRot);

  //  // Write elements file
  //  nElements *= 4;
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  for (unsigned e=0; e<nElements; e+=4)
  //  {
  //    elementsOutput << e+1 << '\t' << e+2 << '\t' << e+3 << '\t' << e+7 << std::endl;
  //    elementsOutput << e+2 << '\t' << e+4 << '\t' << e+6 << '\t' << e+7 << std::endl;
  //    elementsOutput << e+1 << '\t' << e+5 << '\t' << e+4 << '\t' << e+7 << std::endl;
  //    elementsOutput << e+1 << '\t' << e+4 << '\t' << e+2 << '\t' << e+7 << std::endl;
  //  }
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = nElements+4;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  vRot *= M_PI/30.; // convert rpm to rad/s
  //  double kinematicViscosity = mitrem->getSolutionKinematicViscosity();
  //  double Quot = vRot/kinematicViscosity;
  //  double sqrtQuot = sqrt(Quot);
  //  double sqrtProd = sqrt(vRot*kinematicViscosity);
  //  if (refFac > 1.)
  //  {
  //    double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
  //    double refFacN = refFac;
  //    double refFacm = 1.;
  //    for (unsigned p=1; p<nElements/3; p++)
  //    {
  //      refFacN *= refFac;
  //    }
  //    for (unsigned m=0; m<nNodes/3; m++)
  //    {
  //      double x = length*(1-refFacm)/(1-refFacN);
  //      double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
  //      nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
  //      nodesOutput << x << '\t' << H << '\t' << 0.0 << std::endl;
  //      nodesOutput << x << '\t' << 0 << '\t' << H << std::endl;
  //      nodesOutput << x << '\t' << H << '\t' << H << std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      refFacm *= refFac;
  //    }
  //  }
  //  else if (refFac == 1.)
  //  {
  //    double H = 3.*length/nElements;
  //    for (unsigned m=0; m<nNodes/3; m++)
  //    {
  //      double x = 3.*length*m/nElements;
  //      double v = sqrtProd*Quot*x*x*(-0.51023 + 0.333*sqrtQuot*x - 0.103*Quot*x*x);
  //      nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
  //      nodesOutput << x << '\t' << H << '\t' << 0.0 << std::endl;
  //      nodesOutput << x << '\t' << 0 << '\t' << H << std::endl;
  //      nodesOutput << x << '\t' << H << '\t' << H << std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //      flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
  //    }
  //  }
  //  else
  //  {
  //    // error message
  //  }
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 4;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 0 << '\t' << 2 << '\t' << 1 << std::endl;
  //  boundaryElementsOutput << 1 << '\t' << 2 << '\t' << 3 << std::endl;
  //  boundaryElementsOutput << nNodes-3 << '\t' << nNodes-2 << '\t' << nNodes-1 << std::endl;
  //  boundaryElementsOutput << nNodes-3 << '\t' << nNodes-2 << '\t' << nNodes-1 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 2;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<dimensions> = Electrode\n" << std::endl;
  //  boundariesOutput << "\t<dimensions> = VirtualElectrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 2" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements1] = 2" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 1;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}

  else if (dimensions == "ConvectionDiffusion_1D")
  {
    dimensions = "1D";
    nDimensions = 1;

    double length,refFac,v;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[v]",v);

    // Write elements file
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e++)
    {
      elementsOutput << e << '\t' << e+1 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+1;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    if (refFac > 1.)
    {
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        nodesOutput << x << std::endl;
        flowFieldOutput << v << std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      for (unsigned m=0; m<nNodes; m++)
      {
        double x = length*m/nElements;
        nodesOutput << x << std::endl;
        flowFieldOutput << v << std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 0 << std::endl;
    boundaryElementsOutput << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Diffusion\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 0;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput.close();
  }

  else if (dimensions == "ConvectionDiffusion_2D")
  {
    dimensions = "2D";
    nDimensions = 2;

    double length,refFac,v;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[v]",v);

    // Write elements file
    nElements *= 2;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=4)
    {
      elementsOutput << e << '\t' << e+2 << '\t' << e+1 << std::endl;
      elementsOutput << e+1 << '\t' << e+2 << '\t' << e+3 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+2;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/2; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 1e-5 << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 2.*length/nElements;
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = 2.*length*m/nElements;
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 1e-5 << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 1 << '\t' << 0 << std::endl;
    boundaryElementsOutput << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Diffusion\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 0;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput.close();
  }

  else if (dimensions == "ConvectionDiffusion_AX")
  {
    dimensions = "AX";
    nDimensions = 2;

    double length,refFac,v;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[v]",v);

    // Write elements file
    nElements *= 2;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=2)
    {
      elementsOutput << e << '\t' << e+2 << '\t' << e+1 << std::endl;
      elementsOutput << e+1 << '\t' << e+2 << '\t' << e+3 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+2;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/2; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        nodesOutput << x << '\t' << 0 << std::endl;
        nodesOutput << x << '\t' << 1 << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 2.*length/nElements;
      for (unsigned m=0; m<nNodes/2; m++)
      {
        double x = 2.*length*m/nElements;
        nodesOutput << x << '\t' << 0 << std::endl;
        nodesOutput << x << '\t' << 1 << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 1 << '\t' << 0 << std::endl;
    boundaryElementsOutput << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Diffusion\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 0;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput.close();
  }

  else if (dimensions == "ConvectionDiffusion_3D")
  {
    dimensions = "3D";
    nDimensions = 3;

    double length,refFac,v;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[v]",v);

    // Write elements file
    nElements *= 3;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=3)
    {
      elementsOutput << e << '\t' << e+2 << '\t' << e+1 << '\t' << e+3 << std::endl;
      elementsOutput << e+1 << '\t' << e+4 << '\t' << e+3 << '\t' << e+5 << std::endl;
      elementsOutput << e+1 << '\t' << e+3 << '\t' << e+2 << '\t' << e+5 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = nElements+3;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/4.)*(1.-refFac)/(1.-pow(refFac,nElements/2.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/3; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/3; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 0 << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 3.*length/nElements;
      for (unsigned m=0; m<nNodes/3; m++)
      {
        double x = 3.*length*m/nElements;
        nodesOutput << x << '\t' << 0.0 << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << H << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 0 << '\t' << H << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 2;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 0 << '\t' << 2 << '\t' << 1 << std::endl;
    boundaryElementsOutput << nNodes-3 << '\t' << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Diffusion\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 0;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput.close();
  }

  else if (dimensions == "ConvectionDiffusion_2Dbis")
  {
    dimensions = "2D";
    nDimensions = 2;

    double length,refFac,v;

    datFile.readScalar("[nElements]",nElements);
    datFile.readScalar("[length]",length);
    datFile.readScalar("[refFac]",refFac);
    datFile.readScalar("[v]",v);

    // Write elements file
    nElements *= 4;
    std::ofstream elementsOutput(elementsFile.c_str());
    elementsOutput << "Versie 1.0\n" << nElements << std::endl;
    for (unsigned e=0; e<nElements; e+=4)
    {
      unsigned e_ = e*3/4;
      elementsOutput << e_ << '\t' << e_+3 << '\t' << e_+1 << std::endl;
      elementsOutput << e_+1 << '\t' << e_+3 << '\t' << e_+4 << std::endl;
      elementsOutput << e_+1 << '\t' << e_+4 << '\t' << e_+2 << std::endl;
      elementsOutput << e_+2 << '\t' << e_+4 << '\t' << e_+5 << std::endl;
    }
    elementsOutput.close();

    // Write nodes and flow field file
    std::ofstream nodesOutput(nodesFile.c_str());
    std::ofstream flowFieldOutput(flowFieldFile.c_str());
    nNodes = (nElements/4+1)*3;
    nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
    flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
    if (refFac > 1.)
    {
      double H = length*pow(refFac,-0.5+nElements/8.)*(1.-refFac)/(1.-pow(refFac,nElements/4.));
      double refFacN = refFac;
      double refFacm = 1.;
      for (unsigned p=1; p<nElements/4; p++)
      {
        refFacN *= refFac;
      }
      for (unsigned m=0; m<nNodes/3; m++)
      {
        double x = length*(1-refFacm)/(1-refFacN);
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 1.0 << std::endl;
        nodesOutput << x << '\t' << 2.0 << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        refFacm *= refFac;
      }
    }
    else if (refFac == 1.)
    {
      double H = 4.*length/nElements;
      for (unsigned m=0; m<nNodes/3; m++)
      {
        double x = 4.*length*m/nElements;
        nodesOutput << x << '\t' << 0.0 << std::endl;
        nodesOutput << x << '\t' << 1.0 << std::endl;
        nodesOutput << x << '\t' << 2.0 << std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
        flowFieldOutput << v << '\t' << 0.0 <<  std::endl;
      }
    }
    else
    {
      // error message
    }
    nodesOutput.close();
    flowFieldOutput.close();

    // Write boundary elements file
    std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
    nBoundaryElements = 4;
    boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
    boundaryElementsOutput << 2 << '\t' << 1 << std::endl;
    boundaryElementsOutput << 1 << '\t' << 0 << std::endl;
    boundaryElementsOutput << nNodes-3 << '\t' << nNodes-2 << std::endl;
    boundaryElementsOutput << nNodes-2 << '\t' << nNodes-1 << std::endl;
    boundaryElementsOutput.close();

    // Write boundaries file
    std::ofstream boundariesOutput(boundariesFile.c_str());
    nBoundaries = 2;
    boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
    boundariesOutput << "\t<type> = Diffusion\n" << std::endl;
    boundariesOutput << "\t<type> = VirtualElectrode\n" << std::endl;
    boundariesOutput << "[nBoundaryElements0] = 2" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
    boundariesOutput << "[nBoundaryElements1] = 2" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 2\n" << std::endl;
    boundariesOutput << "\t<boundaryElement> = 3\n" << std::endl;
    boundariesOutput.close();

    // Write electrodes file
    std::ofstream electrodesOutput(electrodesFile.c_str());
    nElectrodes = 0;
    electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
    electrodesOutput.close();
  }

  //else if (dimensions == "ParallelPlate_1D")
  //{
  //  dimensions = "1D";
  //  nDimensions = 1;
  //
  //  double length,refFac;
  //
  //  datFile.readScalar("[nElements]",nElements);
  //  datFile.readScalar("[length]",length);
  //  datFile.readScalar("[refFac]",refFac);

  //  // Write elements file
  //  nElements *= 2;
  //  std::ofstream elementsOutput(elementsFile.c_str());
  //  elementsOutput << "Versie 1.0\n" << nElements << std::endl;
  //  for (unsigned e=0; e<nElements; e++)
  //  {
  //    elementsOutput << e << '\t' << e+1 << std::endl;
  //  }
  //  elementsOutput.close();

  //  // Write nodes and flow field file
  //  std::ofstream nodesOutput(nodesFile.c_str());
  //  std::ofstream flowFieldOutput(flowFieldFile.c_str());
  //  nNodes = nElements+1;
  //  nodesOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  flowFieldOutput << "Versie 1.0\n" << nNodes << std::endl;
  //  if (refFac > 1.)
  //  {
  //    // left side
  //    double refFacN = refFac;
  //    double refFacm = 1.;
  //    for (unsigned p=1; p<nElements/2; p++)
  //    {
  //      refFacN *= refFac;
  //    }
  //    for (unsigned m=0; m<nElements/2+1; m++)
  //    {
  //      double x = 0.5*length*((1.-refFacm)/(1.-refFacN) - 1.);
  //      nodesOutput << x << std::endl;
  //      flowFieldOutput << 0 << std::endl;
  //      refFacm *= refFac;
  //    }

  //    // right side
  //    refFacm = refFacN/refFac;
  //    for (unsigned m=0; m<nElements/2; m++)
  //    {
  //      double x = 0.5*length*(1. - (1.-refFacm)/(1.-refFacN));
  //      nodesOutput << x << std::endl;
  //      flowFieldOutput << 0 << std::endl;
  //      refFacm /= refFac;
  //    }
  //  }
  //  else if (refFac == 1.)
  //  {
  //    for (unsigned m=0; m<nNodes; m++)
  //    {
  //      double x = -length*(0.5 - m/nElements);
  //      nodesOutput << x << std::endl;
  //      flowFieldOutput << 0 << std::endl;
  //    }
  //  }
  //  else
  //  {
  //    // error message
  //  }
  //  nodesOutput.close();
  //  flowFieldOutput.close();

  //  // Write boundary elements file
  //  std::ofstream boundaryElementsOutput(boundaryElementsFile.c_str());
  //  nBoundaryElements = 2;
  //  boundaryElementsOutput << "Versie 1.0\n" << nBoundaryElements << std::endl;
  //  boundaryElementsOutput << 0 << std::endl;
  //  boundaryElementsOutput << nNodes-1 << std::endl;
  //  boundaryElementsOutput.close();

  //  // Write boundaries file
  //  std::ofstream boundariesOutput(boundariesFile.c_str());
  //  nBoundaries = 2;
  //  boundariesOutput << "Versie 1.0\n[nBoundaries] = " << nBoundaries << std::endl;
  //  boundariesOutput << "\t<type> = Electrode\n" << std::endl;
  //  boundariesOutput << "\t<type> = Electrode\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements0] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 0\n" << std::endl;
  //  boundariesOutput << "[nBoundaryElements1] = 1" << std::endl;
  //  boundariesOutput << "\t<boundaryElement> = 1\n" << std::endl;
  //  boundariesOutput.close();

  //  // Write electrodes file
  //  std::ofstream electrodesOutput(electrodesFile.c_str());
  //  nElectrodes = 2;
  //  electrodesOutput << "Versie 1.0\n[nElectrodes] = " << nElectrodes << std::endl;
  //  electrodesOutput << "\t<label> = WorkingElectrode\n" << std::endl;
  //  electrodesOutput << "\t<label> = CounterElectrode\n" << std::endl;
  //  electrodesOutput << "[nBoundaries0] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 0\n" << std::endl;
  //  unsigned nElecReactions = mitrem->getNElecReactions();
  //  electrodesOutput << "[nElecReactions0] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput << "[nBoundaries1] = 1" << std::endl;
  //  electrodesOutput << "\t<boundary> = 1\n" << std::endl;
  //  electrodesOutput << "[nElecReactions1] = " << nElecReactions << std::endl;
  //  for (unsigned r=0; r<nElecReactions; r++)
  //  {
  //    electrodesOutput << "\t<label> = " << mitrem->getElecReactionLabel(r) << std::endl;
  //  }
  //  electrodesOutput.close();
  //}

  else
  {
    std::cout << "Invalid [dimensions] in reactor file." << std::endl;
    system("pause");
    exit(1);
  }
}
//---------------------------------------------------------------------------
void Reactor::readNodes ()
{
  size = nNodes*nVariables;

  char temp[256];

  std::ifstream istrm(nodesFile.c_str());
  if (istrm == NULL) errorFileDoesNotExist(nodesFile);
  istrm.setf(std::ios::scientific);
  // Skip version info
  istrm.getline(temp,256,'\n');
  // Read nNodes
  istrm.getline(temp,256,'\n');
  nNodes = atoi(temp);
  nodes = new Vector*[nNodes];
  // Read coordinates
  for (unsigned m=0; m<nNodes; m++)
  {
    nodes[m] = new Vector(nDimensions);
    for (unsigned c=0; c<nDimensions-1; c++)
    {
      istrm.getline(temp,256,'\t');
      nodes[m]->setComponents(c,atof(temp));
    }
    istrm.getline(temp,256,'\n');
    nodes[m]->setComponents(nDimensions-1,atof(temp));
  }
  istrm.close();
}
//---------------------------------------------------------------------------
void Reactor::readElements ()
{
  nElementNodes = nDimensions+1;

  char temp[256];

  std::ifstream istrm(elementsFile.c_str());
  if (istrm == NULL) errorFileDoesNotExist(elementsFile);
  istrm.setf(std::ios::scientific);
  // Skip version info
  istrm.getline(temp,256,'\n');
  // Read nElements
  istrm.getline(temp,256,'\n');
  nElements = atoi(temp);
  elements = new Element*[nElements];
  // Read elements
  for (unsigned e=0; e<nElements; e++)
  {
    elements[e] = new Element(nElementNodes);
    for (unsigned m=0; m<nElementNodes-1; m++)
    {
      istrm.getline(temp,256,'\t');
      elements[e]->setNodes(m,atoi(temp));
    }
    istrm.getline(temp,256,'\n');
    elements[e]->setNodes(nElementNodes-1,atoi(temp));
  }
  istrm.close();
}
//---------------------------------------------------------------------------
void Reactor::readBoundaryElements ()
{
  nBoundaryElementNodes = nDimensions;

  char temp[256];

  std::ifstream istrm(boundaryElementsFile.c_str());
  if (istrm == NULL) errorFileDoesNotExist(boundaryElementsFile);
  istrm.setf(std::ios::scientific);
  // Skip version info
  istrm.getline(temp,256,'\n');
  // Read nBoundaryElements
  istrm.getline(temp,256,'\n');
  nBoundaryElements = atoi(temp);
  boundaryElements = new Element*[nBoundaryElements];
  // Read boundary elements
  for (unsigned e=0; e<nBoundaryElements; e++)
  {
    boundaryElements[e] = new Element(nBoundaryElementNodes);
    for (unsigned m=0; m<nBoundaryElementNodes-1; m++)
    {
      istrm.getline(temp,256,'\t');
      boundaryElements[e]->setNodes(m,atoi(temp));
    }
    istrm.getline(temp,256,'\n');
    boundaryElements[e]->setNodes(nBoundaryElementNodes-1,atoi(temp));
  }
  istrm.close();
}
//---------------------------------------------------------------------------
void Reactor::readFlowField ()
{
  char temp[256];

  std::ifstream istrm(flowFieldFile.c_str());
  if (istrm == NULL) errorFileDoesNotExist(flowFieldFile);
  istrm.setf(std::ios::scientific);
  // Skip version info
  istrm.getline(temp,256,'\n');
  // Read nNodes
  istrm.getline(temp,256,'\n');
  unsigned nNodes_ = atoi(temp);
  if (nNodes_ != nNodes) errorConflictingData(nodesFile,flowFieldFile,"NUMBERS OF NODES");
  flowField = new Vector*[nNodes];
  // Read coordinates
  for (unsigned m=0; m<nNodes; m++)
  {
    flowField[m] = new Vector(nDimensions);
    for (unsigned c=0; c<nDimensions-1; c++)
    {
      istrm.getline(temp,256,'\t');
      flowField[m]->setComponents(c,atof(temp));
    }
    istrm.getline(temp,256,'\n');
    flowField[m]->setComponents(nDimensions-1,atof(temp));
  }
  istrm.close();
}
//---------------------------------------------------------------------------
void Reactor::readMagneticField ()
{
  char temp[256];

  std::ifstream istrm(magneticFieldFile.c_str());
  if (istrm == NULL) errorFileDoesNotExist(magneticFieldFile);
  istrm.setf(std::ios::scientific);
  // Skip version info
  istrm.getline(temp,256,'\n');
  // Read nNodes
  istrm.getline(temp,256,'\n');
  unsigned nNodes_ = atoi(temp);
  if (nNodes_ != nNodes) errorConflictingData(nodesFile,magneticFieldFile,"NUMBERS OF NODES");
  magneticField = new Vector*[nNodes];
  // Read coordinates
  for (unsigned m=0; m<nNodes; m++)
  {
    magneticField[m] = new Vector(nDimensions);
    for (unsigned c=0; c<nDimensions-1; c++)
    {
      istrm.getline(temp,256,'\t');
      magneticField[m]->setComponents(c,atof(temp));
    }
    istrm.getline(temp,256,'\n');
    magneticField[m]->setComponents(nDimensions-1,atof(temp));
  }
  istrm.close();
}
//---------------------------------------------------------------------------
void Reactor::readBoundaries ()
{
  DatFileReader datFile(boundariesFile);

  nBoundaries = datFile.readMultipleVector_nVectors("[nBoundaries]");
  boundaries = new Boundary*[nBoundaries];
  for (unsigned b=0; b<nBoundaries; b++)
  {
    char buf[10];
    sprintf(buf,"%d",b);
    std::string bString(buf);
    std::string nBoundaryElementsb = "[nBoundaryElements" + bString + ']';
    unsigned nBoundaryElements = datFile.readMultipleVector_nVectors(nBoundaryElementsb);
    boundaries[b] = new Boundary(nBoundaryElements);
    std::string type;
    datFile.readMultipleVector("[nBoundaries]", b, "<type>", type);
    boundaries[b]->setType(type);
    //std::cout << "Boundary " << b << '\n';
    //std::cout << "Type " << type << '\n';
    for (unsigned e=0; e<nBoundaryElements; e++)
    {
      unsigned boundaryElement;
      datFile.readMultipleVector(nBoundaryElementsb, e, "<boundaryElement>", boundaryElement);
      boundaries[b]->setBoundaryElements(e,boundaryElement);
      //std::cout << boundaries[b]->getBoundaryElements(e) << " : ";
      //for (unsigned m=0; m<nBoundaryElementNodes; m++)
      //{
      //  std::cout << boundaryElements[boundaries[b]->getBoundaryElements(e)]->getNodes(m) << '\t';
      //}
      //std::cout << std::endl;
    }
  }
  //system("pause");
}
//---------------------------------------------------------------------------
void Reactor::readElectrodes ()
{
  DatFileReader datFile(electrodesFile);

  nElectrodes = datFile.readMultipleVector_nVectors("[nElectrodes]");
  electrodes = new Electrode*[nElectrodes];
  for (unsigned e=0; e<nElectrodes; e++)
  {
    char buf[10];
    sprintf(buf,"%d",e);
    std::string eString(buf);
    std::string nBoundariese = "[nBoundaries" + eString + ']';
    unsigned nBoundaries = datFile.readMultipleVector_nVectors(nBoundariese);

    std::string nElecReactionse = "[nElecReactions" + eString + ']';
    unsigned nElecReactions = datFile.readMultipleVector_nVectors(nElecReactionse);

    std::string nGasReactionse = "[nGasReactions" + eString + ']';
    unsigned nGasReactions = datFile.readMultipleVector_nVectors(nGasReactionse);

    electrodes[e] = new Electrode(nBoundaries,nElecReactions,nGasReactions);

    std::string label;
    datFile.readMultipleVector("[nElectrodes]", e, "<label>", label);
    electrodes[e]->setLabel(label);
    //double V;
    //datFile.readMultipleVector("[nElectrodes]", e, "<V>", V);
    //electrodes[e]->setPotential(V);

    for (unsigned b=0; b<nBoundaries; b++)
    {
      unsigned boundary;
      datFile.readMultipleVector(nBoundariese, b, "<boundary>", boundary);
      electrodes[e]->setBoundaries(b,boundary);
    }

    unsigned nElecReactionsMax = mitrem->getNElecReactions();
    for (unsigned r=0; r<nElecReactions; r++) {
      electrodes[e]->setElecReactions(r,nElecReactionsMax);
      std::string ElecReacrLabel;
      datFile.readMultipleVector(nElecReactionse, r, "<label>", ElecReacrLabel);
      for (unsigned s=0; s<nElecReactionsMax; s++)
      {
        if (mitrem->getElecReactionLabel(s) == ElecReacrLabel)
        {
          electrodes[e]->setElecReactions(r,s);
          break;
        }
      }
      if (electrodes[e]->getElecReactions(r) == nElecReactionsMax)
        errorConflictingData(electrodesFile,".elecreactions"/*mitrem->elecReactionsFile*/,"ELECTRODE REACTION LABELS");
    }

    unsigned nGasReactionsMax = mitrem->getNGasReactions();
    for (unsigned r=0; r<nGasReactions; r++) {
      electrodes[e]->setGasReactions(r,nGasReactionsMax);
      std::string GasReacrLabel;
      datFile.readMultipleVector(nGasReactionse, r, "<label>", GasReacrLabel);
      for (unsigned s=0; s<nGasReactionsMax; s++)
      {
        if (mitrem->getGasReactionLabel(s) == GasReacrLabel)
        {
          electrodes[e]->setGasReactions(r,s);
          break;
        }
      }
      if (electrodes[e]->getGasReactions(r) == nGasReactionsMax)
        errorConflictingData(electrodesFile,".elecreactions"/*mitrem->elecReactionsFile*/,"ELECTRODE REACTION LABELS");
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::readSolver ()
{
  DatFileReader datFile(solverFile);

  datFile.readScalar("[residuMax]",residuMax);
  datFile.readScalar("[iterMax]",iterMax);
  datFile.readScalar("[relaxFacConc]",relaxFacConc);
  datFile.readScalar("[relaxFacConcInc]",relaxFacConcInc);
  datFile.readScalar("[relaxFacPot]",relaxFacPot);
  datFile.readScalar("[relaxFacPotInc]",relaxFacPotInc);

  // Read numerical schemes
  std::string convectionScheme;
  std::string diffusionScheme;
  std::string migrationScheme;
  std::string magneticScheme;
  std::string homReactionScheme;
  std::string electrostaticsScheme;
  std::string timeScheme;
  std::string elecReactionScheme;
  std::string gasReactionScheme;

  datFile.readScalar("[convectionScheme]",convectionScheme);
  datFile.readScalar("[diffusionScheme]",diffusionScheme);
  datFile.readScalar("[migrationScheme]",migrationScheme);
  datFile.readScalar("[magneticScheme]",magneticScheme);
  datFile.readScalar("[homReactionScheme]",homReactionScheme);
  datFile.readScalar("[electrostaticsScheme]",electrostaticsScheme);
  datFile.readScalar("[timeScheme]",timeScheme);
  datFile.readScalar("[elecReactionScheme]",elecReactionScheme);
  datFile.readScalar("[gasReactionScheme]",gasReactionScheme);

  elementMatrixAssembler = new ElementMatrixAssembler(dimensions,mitrem,convectionScheme.c_str(),diffusionScheme.c_str(),migrationScheme.c_str(),magneticScheme.c_str(),homReactionScheme.c_str(),electrostaticsScheme.c_str(),timeScheme.c_str(),elecReactionScheme.c_str(),gasReactionScheme.c_str());
  // read linear solver
  std::string linearSolver;

  datFile.readScalar("[linearSolver]",linearSolver);

  if (linearSolver == "BandLU")
  {
    AMat = new SolverLinear_BandLU(elements,nElements,nNodes,nVariables);
  }
  else if (linearSolver == "GMRES")
  {
    AMat = new SolverLinear_GMRES(elements,nElements,nNodes,nVariables);
  }
  else
  {
    // error message
  }
}
//---------------------------------------------------------------------------
void Reactor::readElectrodePotentials()
{
  DatFileReader datFile(electrodePotentialsFile);

  datFile.readScalar("[timeAccurate]",timeAccurate);
  nElectrodePotentials = datFile.readMultipleVector_nVectors("[nElectrodePotentials]");
  electrodePotentials = new ElectrodePotential[nElectrodePotentials];
  for (unsigned p=0; p<nElectrodePotentials; p++)
  {
    datFile.readMultipleVector("[nElectrodePotentials]",p,"<t>",electrodePotentials[p].t);
    electrodePotentials[p].V = new double[nElectrodes];
    for (unsigned e=0; e<nElectrodes; e++)
    {
      std::string VElectrode = "<V" + electrodes[e]->getLabel() + ">";
      datFile.readMultipleVector("[nElectrodePotentials]",p,VElectrode,electrodePotentials[p].V[e]);
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::solve()
{
  //double kappa = mitrem->calcTransportConductivity();
  //std::cout << "c+ = " << inletVec[0] /*mitrem->getIonConcentration(0)*/ << std::endl;
  //std::cout << "c- = " << inletVec[1] /*mitrem->getIonConcentration(1)*/ << std::endl;
  //std::cout << "conductivity = " << kappa << std::endl;
  std::cout << "Solve..." << std::endl;
  if (timeAccurate)
  {
    for (unsigned m=0; m<nNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        xVec[var(m,i)] = inletVec[i];
        xVecOld[var(m,i)] = inletVec[i];
      }
    }
    dtOld = 1.;
    dt = 1.;
    solve_TimeAccurate(0);
    if (outputFormat == "Excel")
    {
      write_TimeAccurate_Excel(0);
    }
    else if (outputFormat == "TecPlot")
    {
      write_TecPlot(0);
    }
    for (unsigned p=1; p<nElectrodePotentials; p++)
    {
      dtOld = dt;
      dt = electrodePotentials[p].t - electrodePotentials[p-1].t;
      solve_TimeAccurate(p);
      if (outputFormat == "Excel")
      {
        write_TimeAccurate_Excel(p);
      }
      else if (outputFormat == "TecPlot")
      {
        write_TecPlot(p);
      }
    }
  }

  else
  {
    std::cout << "Init" << std::endl;
    if (DoE == false)
    {
      initLinearSystem_Stationary(); //remove here to start from bulk solution
    }
    for (unsigned p=0; p<nElectrodePotentials; p++)
    {
      std::cout << "Solve stationary" << std::endl;
      solve_Stationary(p);
      if (DoE == false)
      {
        if (outputFormat == "Excel")
        {
          std::cout << "Writing output files..." << std::endl;
          write_Stationary_Excel(p);
        }
        else if (outputFormat == "TecPlot")
        {
          std::cout << "Writing output files..." << std::endl;
          write_TecPlot(p);
        }
      }
      else
      {
        calcCurrentDensities();
      }
    }
  }
  //system("pause");
}
//---------------------------------------------------------------------------
void Reactor::solve_Stationary (unsigned p)
{
  std::cout << "Set potential" << std::endl;
  for (unsigned e=0; e<nElectrodes; e++)
  {
    electrodes[e]->setPotential(electrodePotentials[p].V[e]);
  }
  if (DoE == true)
  {
    initLinearSystem_Stationary(); //remove here to start from solution of previously applied potential
  }
  unsigned iter = 0;
  converged = false;
  do
  {
    solveLinearSystem_Stationary(iter);
    unsigned nNegativeConcentrations = 0;
    for (unsigned m=0; m<nNodes; m++)
    {
      for (unsigned i=0; i<nIons; i++)
      {
        if (xVec[var(m,i)] < 0)
        {
          nNegativeConcentrations++;
        }
      }
    }
    iter++;
    std::cout << "Iteration = " << iter << std::endl;
    std::cout << "Residual = " << residu << std::endl;
    std::cout << "Number of negative concentrations = " << nNegativeConcentrations << std::endl;
    std::cout << "Damping factor concentrations = " << relaxFacConc << std::endl;
    std::cout << "Damping factor potential = " << relaxFacPot << std::endl;
    relaxFacConc = std::min(1., relaxFacConc + relaxFacConcInc);
    relaxFacPot = std::min(1., relaxFacPot + relaxFacPotInc);
    //system("pause");

    std::ifstream killProgram("Kill.txt");
    std::ifstream eachIterationOutput;
    eachIterationOutput.open("eachIterationOutput.txt");

    if (eachIterationOutput != NULL)
    {
      std::string originalOutputFile = outputFile;
      char buf[10];
      outputFile = "iteration";// + itoa(iter,buf,10);// + "_" + originalOutputFile;
      outputFile = outputFile + itoa(iter,buf,10) + "_" + originalOutputFile;

      if (outputFormat == "Excel")
      {
        std::cout << "Writing output files..." << std::endl;
        write_Stationary_Excel(p);
      }
      else if (outputFormat == "TecPlot")
      {
        std::cout << "Writing output files..." << std::endl;
        write_TecPlot(p);
      }

      outputFile = originalOutputFile;
    }

    eachIterationOutput.close();

    if (killProgram != NULL)
    {
      converged = true;
    }
  } while ((!converged) && (iter < iterMax));
  std::cout << std::endl;
}
//---------------------------------------------------------------------------
void Reactor::solve_TimeAccurate (unsigned p)
{
  timeFactor1 = 1.5/dt;
  timeFactor2 = 0.5/dtOld;
  timeFactor3 = timeFactor1 + timeFactor2;
  for (unsigned e=0; e<nElectrodes; e++)
  {
    electrodes[e]->setPotential(electrodePotentials[p].V[e]);
  }
  initLinearSystem_TimeAccurate();
  unsigned iter = 0;
  converged = false;
  do
  {
    solveLinearSystem_TimeAccurate();
    unsigned nNegativeConcentrations = 0;
    for (unsigned m=0; m<nNodes; m++)
    {
      for (unsigned i=0; i<nIons; i++)
      {
        if (xVec[var(m,i)] < 0)
        {
          nNegativeConcentrations++;
        }
      }
    }
    iter++;
    std::cout << "Iteration = " << iter << std::endl;
    std::cout << "Residual = " << residu << std::endl;
    std::cout << "Number of negative concentrations = " << nNegativeConcentrations << std::endl;
    std::cout << "Damping factor concentrations = " << relaxFacConc << std::endl;
    std::cout << "Damping factor potential = " << relaxFacPot << std::endl;
    relaxFacConc = std::min(1., relaxFacConc + relaxFacConcInc);
    relaxFacPot = std::min(1., relaxFacPot + relaxFacPotInc);
    //system("pause");
    std::ifstream killProgram("Kill.txt");
    if (killProgram != NULL)
    {
      converged = true;
    }
  } while ((!converged) && (iter < iterMax));
  std::cout << std::endl;
}
//---------------------------------------------------------------------------
void Reactor::write_Stationary_Excel (unsigned p)
{
  calcCurrentDensities();

  char axis[3] = {'x','y','z'};
  unsigned nElecReactions = mitrem->getNElecReactions();
  unsigned nGasReactions = mitrem->getNGasReactions();
  unsigned nHomReactions = mitrem->getNHomReactions();
  char buf[10];

  std::string outputFileV = outputFile + ".xls";
  output.open(outputFileV.c_str(), std::ofstream::out | std::ofstream::app);
  output.precision(12);
  if (p==0)
  {
    output << "V\t";
    for (unsigned i=0; i<nIons; i++)
    {
      output << "c(" << mitrem->getIonLabel(i) << ")\t";
    }
    output << "U\t";
    for (unsigned i=0; i<nIons; i++)
    {
      output << "Jx(" << mitrem->getIonLabel(i) << ")\t";
    }
    output << "Jx\t";
    for (unsigned r=0; r<nElecReactions; r++)
    {
      output << "J(" << mitrem->getElecReactionLabel(r) << ")\t";
    }
    for (unsigned r=0; r<nGasReactions; r++)
    {
      output << "v(" << mitrem->getGasReactionLabel(r) << ")\t";
    }
    output << '\n';
  }
  output << electrodePotentials[p].V[0] << '\t';
  for (unsigned i=0; i<nVariables; i++)
  {
    output << xVec[var(0,i)] << '\t';
  }
  for (unsigned i=0; i<nVariables; i++)
  {
    output << ionCurrentDensities[0][i][0] << '\t';
  }
  for (unsigned r=0; r<nElecReactions; r++)
  {
    output << elecReactionCurrentDensities[0][0][0][0][r] << '\t';
  }
  for (unsigned r=0; r<nGasReactions; r++)
  {
    output << gasReactionRates[0][0][0][0][r] << '\t';
  }
  output << std::endl;
  output.close();


  sprintf(buf,"%f",electrodePotentials[p].V[0]);
  outputFileV = outputFile + "_vars_" + buf + "V.xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << axis[c] << '\t';
  }
  for (unsigned i=0; i<nIons; i++)
  {
    output << "c(" << mitrem->getIonLabel(i) << ")\t";
  }
  output << "U\t";
  for (unsigned r=0; r<nHomReactions; r++)
  {
    output << "v(" << mitrem->getHomReactionLabel(r) << ")\t";
  }
  output << '\n';
  for (unsigned m=0; m<nNodes; m++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << nodes[m]->getComponents(c) << '\t';
    }
    for (unsigned i=0; i<nVariables; i++)
    {
      output << xVec[var(m,i)] << '\t';
    }

    concentrations[0] = &xVec[var(m,0)];
    potentials[0] = xVec[var(m,nIons)];
    mitrem->init(concentrations[0],potentials[0],temperatures[0],densities[0]);
    for (unsigned r=0; r<nHomReactions; r++)
    {
      output << mitrem->calcHomReactionRate(r) << '\t';
    }
    output << std::endl;
  }
  output.close();


  outputFileV = outputFile + "_currentdensities_" + buf + "V.xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << axis[c] << '\t';
  }
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << 'J' << axis[c] << '(' << mitrem->getIonLabel(i) << ")\t";
    }
  }
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << 'J' << axis[c] << '\t';
  }
  output << '\n';
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned c=0; c<nDimensions; c++)
    {
      double coordinate = 0.;
      for (unsigned m=0; m<nElementNodes; m++)
      {
        coordinates[m] = nodes[elementNodes[m]]->getComponents();
        coordinate += coordinates[m][c];
      }
      coordinate /= nElementNodes;
      output << coordinate << '\t';
    }
    for (unsigned i=0; i<nVariables; i++)
    {
      for (unsigned c=0; c<nDimensions; c++)
      {
        output << ionCurrentDensities[e][i][c] << '\t';
      }
    }
    output << std::endl;
  }
  output.close();


  outputFileV = outputFile + "_electrodes_" + buf + "V.xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    output << electrode->getLabel() << "\tI = \t" << currents[e] << "A (3D) or A/m (2D), NONSENSE FOR AXISYMMETRIC!!!\n";
    std::cout << "I(" << electrode->getLabel() << ") = " << currents[e] << "A (3D) or A/m (2D), NONSENSE FOR AXISYMMETRIC!!!\n";
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << axis[c] << '\t';
    }
    for (unsigned r=0; r<nElecReactions; r++)
    {
      output << "J(" << mitrem->getElecReactionLabel(r) << "\t)";
    }
    for (unsigned r=0; r<nGasReactions; r++)
    {
      output << "v(" << mitrem->getGasReactionLabel(r) << "\t)";
    }
    output << '\n';
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          for (unsigned c=0; c<nDimensions; c++)
          {
            output << nodes[boundaryElementNodes[m]]->getComponents(c) << '\t';
          }
          for (unsigned r=0; r<nElecReactions; r++)
          {
            output << elecReactionCurrentDensities[e][b][be][m][r] << '\t';
          }
          for (unsigned r=0; r<nGasReactions; r++)
          {
            output << gasReactionRates[e][b][be][m][r] << '\t';
          }
          output << std::endl;
        }
      }
    }
    output << std::endl;
  }
  output.close();

  //system("pause");
}
//---------------------------------------------------------------------------
void Reactor::write_TimeAccurate_Excel (unsigned p)
{
  calcCurrentDensities();

  char axis[3] = {'x','y','z'};
  unsigned nElecReactions = mitrem->getNElecReactions();
  unsigned nGasReactions = mitrem->getNGasReactions();
  unsigned nHomReactions = mitrem->getNHomReactions();
  char buf[10];

  std::string outputFileV = outputFile + ".xls";
  output.open(outputFileV.c_str(), std::ofstream::out | std::ofstream::app);
  output.precision(12);
  if (p==0)
  {
    output << "t\tV\t";
    for (unsigned i=0; i<nIons; i++)
    {
      output << "c(" << mitrem->getIonLabel(i) << ")\t";
    }
    output << "U\t";
    for (unsigned i=0; i<nIons; i++)
    {
      output << "Jx(" << mitrem->getIonLabel(i) << ")\t";
    }
    output << "Jx\t";
    for (unsigned r=0; r<nElecReactions; r++)
    {
      output << "J(" << mitrem->getElecReactionLabel(r) << ")\t";
    }
    for (unsigned r=0; r<nGasReactions; r++)
    {
      output << "v(" << mitrem->getGasReactionLabel(r) << ")\t";
    }
    output << '\n';
  }
  output << electrodePotentials[p].t << '\t' << electrodePotentials[p].V[0] << '\t';
  for (unsigned i=0; i<nVariables; i++)
  {
    output << xVec[var(0,i)] << '\t';
  }
  for (unsigned i=0; i<nVariables; i++)
  {
    output << ionCurrentDensities[0][i][0] << '\t';
  }
  for (unsigned r=0; r<nElecReactions; r++)
  {
    output << elecReactionCurrentDensities[0][0][0][0][r] << '\t';
  }
  for (unsigned r=0; r<nGasReactions; r++)
  {
    output << gasReactionRates[0][0][0][0][r] << '\t';
  }
  output << std::endl;
  output.close();


  sprintf(buf,"%f",electrodePotentials[p].t);
  outputFileV = outputFile + "_vars_" + buf + "s.xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << axis[c] << '\t';
  }
  for (unsigned i=0; i<nIons; i++)
  {
    output << "c(" << mitrem->getIonLabel(i) << ")\t";
  }
  output << "U\t";
  for (unsigned r=0; r<nHomReactions; r++)
  {
    output << mitrem->getHomReactionLabel(r) << '\t';
  }
  output << '\n';
  for (unsigned m=0; m<nNodes; m++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << nodes[m]->getComponents(c) << '\t';
    }
    for (unsigned i=0; i<nVariables; i++)
    {
      output << xVec[var(m,i)] << '\t';
    }

    concentrations[0] = &xVec[var(m,0)];
    potentials[0] = xVec[var(m,nIons)];
    mitrem->init(concentrations[0],potentials[0],temperatures[0],densities[0]);
    for (unsigned r=0; r<nHomReactions; r++)
    {
      output << mitrem->calcHomReactionRelativeDeviationFromEquilibrium(r) << '\t';
    }
    output << std::endl;
  }
  output.close();


  outputFileV = outputFile + "_currentdensities_" + buf + "s.xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << axis[c] << '\t';
  }
  for (unsigned i=0; i<nIons; i++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << 'J' << axis[c] << '(' << mitrem->getIonLabel(i) << ")\t";
    }
  }
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << 'J' << axis[c] << '\t';
  }
  output << '\n';
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned c=0; c<nDimensions; c++)
    {
      double coordinate = 0.;
      for (unsigned m=0; m<nElementNodes; m++)
      {
        coordinates[m] = nodes[elementNodes[m]]->getComponents();
        coordinate += coordinates[m][c];
      }
      coordinate /= nElementNodes;
      output << coordinate << '\t';
    }
    for (unsigned i=0; i<nVariables; i++)
    {
      for (unsigned c=0; c<nDimensions; c++)
      {
        output << ionCurrentDensities[e][i][c] << '\t';
      }
    }
    output << std::endl;
  }
  output.close();


  outputFileV = outputFile + "_electrodes_" + buf + "s.xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    output << electrode->getLabel() << "\tI = \t" << currents[e] << '\n';
    std::cout << "I(" << electrode->getLabel() << ") = " << currents[e] << '\n';
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << axis[c] << '\t';
    }
    for (unsigned r=0; r<nElecReactions; r++)
    {
      output << "J(" << mitrem->getElecReactionLabel(r) << "\t)";
    }
    for (unsigned r=0; r<nGasReactions; r++)
    {
      output << "v(" << mitrem->getGasReactionLabel(r) << "\t)";
    }
    output << '\n';
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          for (unsigned c=0; c<nDimensions; c++)
          {
            output << nodes[boundaryElementNodes[m]]->getComponents(c) << '\t';
          }
          for (unsigned r=0; r<nElecReactions; r++)
          {
            output << elecReactionCurrentDensities[e][b][be][m][r] << '\t';
          }
          for (unsigned r=0; r<nGasReactions; r++)
          {
            output << gasReactionRates[e][b][be][m][r] << '\t';
          }
          output << std::endl;
        }
      }
    }
    output << std::endl;
  }
  output.close();
  //system("pause");
}
//---------------------------------------------------------------------------
void Reactor::write_TecPlot (unsigned p)
{
  calcCurrentDensities();
  char axis[3] = {'X','Y','Z'};
  std::string velo[3] = {"vx","vy","vz"};
  char buf[10];
  unsigned nElecReactions = mitrem->getNElecReactions();
  unsigned nGasReactions = mitrem->getNGasReactions();
  std::string outputFileV = outputFile + itoa(p,buf,10) + ".plt";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(20);

  // Title
  output << "Title = output" << std::endl;

  // Variable names
  output << "VARIABLES = ";
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << axis[c] << ' ';
  }
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << velo[c] << ' ';
  }
  for (unsigned i=0; i<nIons; i++)
  {
    output << "c(" << mitrem->getIonLabel(i) << ") ";
  }
  output << "U\n";

  // Variabe values
  output << "ZONE N=" << nNodes << "\t,E= " << nElements
    << ",    F=FEPOINT       ET=";
  if (nDimensions == 1)
  {
    output << "LINE" << std::endl;
  }
  else if (nDimensions == 2)
  {
    output << "TRIANGLE" << std::endl;
  }
  else if (nDimensions == 3)
  {
    output << "TETRAHEDRON" << std::endl;
  }

  for (unsigned m=0; m<nNodes; m++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << nodes[m]->getComponents(c) << '\t';
    }
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << flowField[m]->getComponents(c) << '\t';
    }
    for (unsigned i=0; i<nVariables; i++)
    {
      output << xVec[var(m,i)] << '\t';
    }
    output << std::endl;
  }
  for (unsigned e=0; e<nElements; e++)
  {
    for (unsigned m=0; m<nElementNodes; m++)
    {
      output << elements[e]->getNodes(m)+1 << '\t';
    }
    output << std::endl;
  }
  output.close();


  outputFileV = outputFile + "_electrodes_" + itoa(p,buf,10) + ".xls";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);
  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    output << electrode->getLabel() << "\tI = \t" << currents[e] << '\n';
    std::cout << "I(" << electrode->getLabel() << ") = " << currents[e] << '\n';
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << axis[c] << '\t';
    }
    for (unsigned i=0; i<nIons; i++)
    {
      output << "c(" << mitrem->getIonLabel(i) << ")\t";
    }
    output << "U\t";
    for (unsigned r=0; r<nElecReactions; r++)
    {
      output << "J(" << mitrem->getElecReactionLabel(r) << ")\t";
    }
    for (unsigned r=0; r<nGasReactions; r++)
    {
      output << "v(" << mitrem->getGasReactionLabel(r) << ")\t";
    }
    output << '\n';
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        //for (unsigned m=0; m<nBoundaryElementNodes; m++)
        for (unsigned m=nBoundaryElementNodes; m>0; m--)
        {
          for (unsigned c=0; c<nDimensions; c++)
          {
            output << nodes[boundaryElementNodes[m-1]]->getComponents(c) << '\t';
          }
          for (unsigned i=0; i<nVariables; i++)
          {
            output << xVec[var(boundaryElementNodes[m-1],i)] << '\t';
          }
          for (unsigned r=0; r<nElecReactions; r++)
          {
            output << elecReactionCurrentDensities[e][b][be][m-1][r] << '\t';
          }
          for (unsigned r=0; r<nGasReactions; r++)
          {
            output << gasReactionRates[e][b][be][m-1][r] << '\t';
          }
          output << std::endl;
        }
      }
    }
    output << std::endl;
  }
  output.close();
}
//---------------------------------------------------------------------------
void Reactor::write_TecPlot_Residual (unsigned p)
{
  char axis[3] = {'X','Y','Z'};
  char buf[10];
  std::string outputFileV = outputFile + "_Residual" + itoa(p,buf,10) + ".plt";
  output.open(outputFileV.c_str(), std::ofstream::out);
  output.precision(12);

  // Title
  output << "Title = output" << std::endl;

  // Variable names
  output << "VARIABLES = ";
  for (unsigned c=0; c<nDimensions; c++)
  {
    output << axis[c] << ' ';
  }
  for (unsigned i=0; i<nVariables; i++)
  {
    output << "equation" << i << " ";
  }
  output << "\n";

  // Variabe values
  output << "ZONE N=" << nNodes << "\t,E= " << nElements
    << ",    F=FEPOINT       ET=";
  if (nDimensions == 1)
  {
    output << "LINE" << std::endl;
  }
  else if (nDimensions == 2)
  {
    output << "TRIANGLE" << std::endl;
  }
  else if (nDimensions == 3)
  {
    output << "TETRAHEDRON" << std::endl;
  }

  for (unsigned m=0; m<nNodes; m++)
  {
    for (unsigned c=0; c<nDimensions; c++)
    {
      output << nodes[m]->getComponents(c) << '\t';
    }
    for (unsigned i=0; i<nVariables; i++)
    {
      output << AMat->getResidual(xVec, bVec, var(m,i)) << '\t';
    }
    output << std::endl;
  }
  for (unsigned e=0; e<nElements; e++)
  {
    for (unsigned m=0; m<nElementNodes; m++)
    {
      output << elements[e]->getNodes(m)+1 << '\t';
    }
    output << std::endl;
  }
  output.close();
}
//---------------------------------------------------------------------------
void Reactor::initLinearSystem_Stationary ()
{
  char temp[256];

  std::ifstream istrm(initialConditionsFile.c_str());
  if (istrm == NULL)
  {
    std::cout << "Applying BULK concentrations...";
    //std::cin.getline(temp,256,'\n');
    for (unsigned m=0; m<nNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        xVec[var(m,i)] = inletVec[i];
      }
    }
  }
  else
  {
    std::cout << "Applying concentrations from file...";
    //std::cin.getline(temp,256,'\n');
    istrm.setf(std::ios::scientific);
    istrm.getline(temp,256,'\n');
    istrm.getline(temp,256,'\n');
    istrm.getline(temp,256,'\n');

    for (unsigned m=0; m<nNodes; m++)
    {
      for (unsigned c=0; c<nDimensions; c++)
      {
        istrm.getline(temp,256,'\t');//pointcoordinates
      }
      for (unsigned c=0; c<nDimensions; c++)
      {
        istrm.getline(temp,256,'\t');//flowfield
      }
      for (unsigned i=0; i<nVariables; i++)
      {
        istrm.getline(temp,256,'\t');//concentrations
        xVec[var(m,i)] = atof(temp);
      }
      istrm.getline(temp,256,'\n');
    }
  }
  istrm.close();


/*  for (unsigned m=0; m<nNodes; m++) {
    for (unsigned i=0; i<nVariables; i++) {
      xVec[var(m,i)] = inletVec[i];
    }*/


    /*for (unsigned i=0; i<nIons; i++) {
      xVec[var(m,i)] = inletVec[i];
    }
    xVec[var(m,nIons)] = nodes[m]->getComponents(0)*1.;*/
//  }
  //system("pause");
}
//---------------------------------------------------------------------------
void Reactor::initLinearSystem_TimeAccurate ()
{
  // Compute epsilon/dt(-1) * [T(X(t-1))]*{X(t-1)}
  xVecPtr = xVec;
  xVec = xVecOld;
  // XVec points to variables at t-1. This is necessary if the time matrix is concentration dependent, so that addElementMatTime(e) is computed correctly!

  AMat->init();
  addElementTimeMatrices();
  AMat->multiplyWithVector(xVec,bVecOld);
  for (unsigned row=0; row<size; row++)
  {
    bVecOld[row] *= timeFactor2;
  }


  // Compute ((1+epsilon)/dt+epsilon/dt(-1)) * [T(X(t))]*{X(t)}
  xVec = xVecPtr;
  // xVec points to variables at t

  AMat->init();
  addElementTimeMatrices();
  AMat->multiplyWithVector(xVec,bVec);
  for (unsigned row=0; row<size; row++)
  {
    bVecOld[row] -= timeFactor3*bVec[row];
  }


  // The old part is now stored in bVecOld. You no longer need X(t-1).
  for (unsigned row=0; row<size; row++)
  {
    xVecOld[row] = xVec[row];
  }
}
//---------------------------------------------------------------------------
void Reactor::solveLinearSystem_Stationary (unsigned iter)
{
  std::cout << "Init system matrix" << std::endl;
  AMat->init();
  std::cout << "Init rhs vector" << std::endl;
  for (unsigned row=0; row<size; row++)
  {
    bVec[row] = 0.;
  }

  // A.X = B
  std::cout << "Add element matrices" << std::endl;
  addElementMatrices_Stationary();
  std::cout << "Add boundary element vectors" << std::endl;
  addBoundaryElementVectors();

  // R = B - A.X
  //std::string f = outputFile + ".xls";
  //AMat->print(f);
  //std::ofstream CHECK;
  //CHECK.open(f.c_str(), std::ofstream::out | std::ofstream::app);
  //CHECK.precision(12);
  //CHECK << "\nX :\t\tB :" << std::endl;
  //for (unsigned row=0; row<size; row++)
  //{
  //  CHECK << xVec[row] << "\t\t" << bVec[row] << std::endl;
  //}
  //CHECK << std::endl;
  //CHECK.close();
  //exit(1);
  std::cout << "Calc residu" << std::endl;
  AMat->residu(xVec,bVec);
//  AMat->print("matrix2.txt");

  std::ifstream writeResidual("writeResidual.txt");
  if (writeResidual != NULL)
  {
    write_TecPlot_Residual(iter);
  }
  //CHECK.open(f.c_str(), std::ofstream::out | std::ofstream::app);
  //for (unsigned row=0; row<size; row++)
  //  CHECK<< bVec[row] << std::endl;
  //CHECK << std::endl;
  //CHECK.close();
  //~ BVec is now the residual vector R. ~//

  // K = AMat + AJac - dB/dX
  //system("pause");
  std::cout << "Add element jacobians" << std::endl;
  addElementJacobians_Stationary();
  std::cout << "Add boundary element jacobians" << std::endl;
  addBoundaryElementJacobians();

  //~ AMat is now the linear system matrix K. ~//
  //AMat->print(f);

  // Impose essential boundary conditions
  std::cout << "Impose BC" << std::endl;
  imposeEssentialBoundaryConditions();

  // Check convergence
  residu = 0.;
  for (unsigned row=0; row<size; row++)
  {
    residu += bVec[row]*bVec[row];
  }
  residu = sqrt(residu);
  if (residu < residuMax)
  {
    converged = true;
    return;
  }

  //AMat->print(f);
  //exit(1);
  ////std::ofstream CHECK;
  //CHECK.open(f.c_str(), std::ofstream::out | std::ofstream::app);
  //CHECK.precision(12);
  //CHECK << "\nR :" << std::endl;
  //for (unsigned row=0; row<size; row++)
  //{
  //  CHECK << bVec[row] << std::endl;
  //}
  //CHECK << std::endl;
  //CHECK.close();
  //exit(1);

  // Solve the linear system
  std::cout << "Solve linear system" << std::endl;
  AMat->solve(bVec);
  //~ BVec is now the incremental vector DX. ~//

  //for (unsigned row=0; row<size; row++)
  //{
  //  CHECK << bVec[row] << std::endl;
  //}
  //CHECK << std::endl;
  //CHECK.close();
  //exit(1);

  // Update the solution vector
  std::cout << "Update" << std::endl;

  for (unsigned node=0; node<nNodes; node++)
  {
    for (unsigned ion=0; ion<nIons; ion++)
    {
      xVec[node*nVariables+ion] += relaxFacConc*bVec[node*nVariables+ion];
    }
    xVec[node*nVariables+nIons] += relaxFacPot*bVec[node*nVariables+nIons];
  }

/*  for (unsigned row=0; row<size; row++)
  {
    xVec[row] += relaxFac*bVec[row];
  }*/
}
//---------------------------------------------------------------------------
void Reactor::solveLinearSystem_TimeAccurate ()
{
  AMat->init();
  for (unsigned row=0; row<size; row++)
  {
    bVec[row] = bVecOld[row];
  }

  // A.X = B
  addElementMatrices_TimeAccurate();
  addBoundaryElementVectors();

  // R = B - A.X
  //std::string f = outputFile + ".xls";
  ////AMat->print(f);
  //std::ofstream CHECK;
  //CHECK.open(f.c_str(), std::ofstream::out | std::ofstream::app);
  //CHECK.precision(12);
  //CHECK << "\nX :\t\tB :" << std::endl;
  //for (unsigned row=0; row<size; row++)
  //{
  //  CHECK << xVec[row] << "\t\t" << bVec[row] << std::endl;
  //}
  //CHECK << std::endl;
  //exit(1);
  AMat->residu(xVec,bVec);
  //for (unsigned row=0; row<size; row++)
  //  CHECK<< bVec[row] << std::endl;
  //CHECK << std::endl;
  //CHECK.close();
  //~ BVec is now the residual vector R. ~//

  // K = AMat + AJac - dB/dX
  //system("pause");
  addElementJacobians_TimeAccurate();
  addBoundaryElementJacobians();

  //~ AMat is now the linear system matrix K. ~//
  //AMat->print(f);

  // Impose essential boundary conditions
  imposeEssentialBoundaryConditions();

  // Check convergence
  residu = 0.;
  for (unsigned row=0; row<size; row++)
  {
    residu += bVec[row]*bVec[row];
  }
  residu = sqrt(residu);
  if (residu < residuMax)
  {
    converged = true;
    return;
  }

  //AMat->print(f);
  //std::ofstream CHECK;
  //CHECK.open(f.c_str(), std::ofstream::out | std::ofstream::app);
  //CHECK.precision(12);
  //for (unsigned row=0; row<size; row++)
  //{
  //  CHECK << bVec[row] << std::endl;
  //}
  //CHECK << std::endl;

  // Solve the linear system
  AMat->solve(bVec);
  //~ BVec is now the incremental vector DX. ~//

  //for (unsigned row=0; row<size; row++)
  //{
  //  CHECK << bVec[row] << std::endl;
  //}
  //CHECK << std::endl;
  //CHECK.close();
  //exit(1);

  // Update the solution vector
  for (unsigned node=0; node<nNodes; node++)
  {
    for (unsigned ion=0; ion<nIons; ion++)
    {
      xVec[node*nVariables+ion] += relaxFacConc*bVec[node*nVariables+ion];
    }
    xVec[node*nVariables+nIons] += relaxFacPot*bVec[node*nVariables+nIons];
  }

/*  for (unsigned row=0; row<size; row++)
  {
    xVec[row] += relaxFac*bVec[row];
  }*/
}
//---------------------------------------------------------------------------
void Reactor::addElementMatrices_Stationary ()
{
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned m=0; m<nElementNodes; m++)
    {
      coordinates[m] = nodes[elementNodes[m]]->getComponents();
      velocities[m] = flowField[elementNodes[m]]->getComponents();
      magneticFieldVectors[m] = magneticField[elementNodes[m]]->getComponents();
      concentrations[m] = &xVec[var(elementNodes[m],0)];
      potentials[m] = xVec[var(elementNodes[m],nIons)];
    }
    DoubleMatrix elementMat = elementMatrixAssembler->calcElementMat(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    for (unsigned m=0; m<nElementNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        unsigned rowLocal = eq(m,i);//var(m,i);
        unsigned rowGlobal = var(elementNodes[m],i);


        for (unsigned n=0; n<nElementNodes; n++)
        {
          for (unsigned j=0; j<nVariables; j++)
          {
            unsigned colLocal = var(n,j);
            unsigned colGlobal = var(elementNodes[n],j);
            AMat->add(rowGlobal,colGlobal,elementMat[rowLocal][colLocal]);
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::addElementMatrices_TimeAccurate ()
{
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned m=0; m<nElementNodes; m++)
    {
      coordinates[m] = nodes[elementNodes[m]]->getComponents();
      velocities[m] = flowField[elementNodes[m]]->getComponents();
      magneticFieldVectors[m] = magneticField[elementNodes[m]]->getComponents();
      concentrations[m] = &xVec[var(elementNodes[m],0)];
      potentials[m] = xVec[var(elementNodes[m],nIons)];
    }
    DoubleMatrix elementMat = elementMatrixAssembler->calcElementMat(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    DoubleMatrix elementTimeMat = elementMatrixAssembler->calcElementTimeMat(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    for (unsigned m=0; m<nElementNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        unsigned rowLocal = eq(m,i);//var(m,i);
        unsigned rowGlobal = var(elementNodes[m],i);
        for (unsigned n=0; n<nElementNodes; n++)
        {
          for (unsigned j=0; j<nVariables; j++)
          {
            unsigned colLocal = var(n,j);
            unsigned colGlobal = var(elementNodes[n],j);
            double contribution = - timeFactor1*elementTimeMat[rowLocal][colLocal] + elementMat[rowLocal][colLocal];
            AMat->add(rowGlobal,colGlobal,contribution);
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::addElementJacobians_Stationary ()
{
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned m=0; m<nElementNodes; m++)
    {
      coordinates[m] = nodes[elementNodes[m]]->getComponents();
      velocities[m] = flowField[elementNodes[m]]->getComponents();
      magneticFieldVectors[m] = magneticField[elementNodes[m]]->getComponents();
      concentrations[m] = &xVec[var(elementNodes[m],0)];
      potentials[m] = xVec[var(elementNodes[m],nIons)];
    }
    DoubleMatrix elementJac = elementMatrixAssembler->calcElementJac(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    for (unsigned m=0; m<nElementNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        unsigned rowLocal = eq(m,i);//var(m,i);
        unsigned rowGlobal = var(elementNodes[m],i);
        for (unsigned n=0; n<nElementNodes; n++)
        {
          for (unsigned j=0; j<nVariables; j++)
          {
            unsigned colLocal = var(n,j);
            unsigned colGlobal = var(elementNodes[n],j);
            AMat->add(rowGlobal,colGlobal,elementJac[rowLocal][colLocal]);
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::addElementJacobians_TimeAccurate ()
{
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned m=0; m<nElementNodes; m++)
    {
      coordinates[m] = nodes[elementNodes[m]]->getComponents();
      velocities[m] = flowField[elementNodes[m]]->getComponents();
      magneticFieldVectors[m] = magneticField[elementNodes[m]]->getComponents();
      concentrations[m] = &xVec[var(elementNodes[m],0)];
      potentials[m] = xVec[var(elementNodes[m],nIons)];
    }
    DoubleMatrix elementJac = elementMatrixAssembler->calcElementJac(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    DoubleMatrix elementTimeJac = elementMatrixAssembler->calcElementTimeJac(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    for (unsigned m=0; m<nElementNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        unsigned rowLocal = eq(m,i);//var(m,i);
        unsigned rowGlobal = var(elementNodes[m],i);
        for (unsigned n=0; n<nElementNodes; n++)
        {
          for (unsigned j=0; j<nVariables; j++)
          {
            unsigned colLocal = var(n,j);
            unsigned colGlobal = var(elementNodes[n],j);
            double contribution =  - timeFactor1*elementTimeJac[rowLocal][colLocal] + elementJac[rowLocal][colLocal];
            AMat->add(rowGlobal,colGlobal,contribution);
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::addBoundaryElementVectors ()
{
  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          coordinates[m] = nodes[boundaryElementNodes[m]]->getComponents();
          concentrations[m] = &xVec[var(boundaryElementNodes[m],0)];
          potentials[m] = xVec[var(boundaryElementNodes[m],nIons)];
        }
        unsigned nElecReactions = electrode->getNElecReactions();
        IndexList elecReactions = electrode->getElecReactions();
        double electrodePotential = electrode->getPotential();
        unsigned nGasReactions = electrode->getNGasReactions();
        IndexList gasReactions = electrode->getGasReactions();
        DoubleVector boundaryElementVec = elementMatrixAssembler->calcBoundaryElementVec(coordinates, concentrations, potentials, temperatures, densities, voidFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          for (unsigned i=0; i<nVariables; i++)
          {
            unsigned rowLocal = eq(m,i);//var(m,i);
            unsigned rowGlobal = var(boundaryElementNodes[m],i);
            bVec[rowGlobal] += boundaryElementVec[rowLocal];
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::addBoundaryElementJacobians ()
{
  for (unsigned e=0; e<nElectrodes; e++)
  {
    Electrode* electrode = electrodes[e];
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {

          coordinates[m] = nodes[boundaryElementNodes[m]]->getComponents();
          concentrations[m] = &xVec[var(boundaryElementNodes[m],0)];
          potentials[m] = xVec[var(boundaryElementNodes[m],nIons)];
        }
        unsigned nElecReactions = electrode->getNElecReactions();
        IndexList elecReactions = electrode->getElecReactions();
        double electrodePotential = electrode->getPotential();
        unsigned nGasReactions = electrode->getNGasReactions();
        IndexList gasReactions = electrode->getGasReactions();
        DoubleMatrix boundaryElementJac = elementMatrixAssembler->calcBoundaryElementJac(coordinates, concentrations, potentials, temperatures, densities, voidFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          for (unsigned i=0; i<nVariables; i++)
          {
            unsigned rowLocal = eq(m,i);//var(m,i);
            unsigned rowGlobal = var(boundaryElementNodes[m],i);
            for (unsigned n=0; n<nBoundaryElementNodes; n++)
            {
              for (unsigned j=0; j<nVariables; j++)
              {
                unsigned colLocal = var(n,j);
                unsigned colGlobal = var(boundaryElementNodes[n],j);
                AMat->add(rowGlobal,colGlobal,boundaryElementJac[rowLocal][colLocal]);
              }
            }
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::addElementTimeMatrices ()
{
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned m=0; m<nElementNodes; m++)
    {
      coordinates[m] = nodes[elementNodes[m]]->getComponents();
      velocities[m] = flowField[elementNodes[m]]->getComponents();
      magneticFieldVectors[m] = magneticField[elementNodes[m]]->getComponents();
      concentrations[m] = &xVec[var(elementNodes[m],0)];
      potentials[m] = xVec[var(elementNodes[m],nIons)];
    }
    DoubleMatrix elementTimeMat = elementMatrixAssembler->calcElementTimeMat(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    for (unsigned m=0; m<nElementNodes; m++)
    {
      for (unsigned i=0; i<nVariables; i++)
      {
        unsigned rowLocal = eq(m,i);//var(m,i);
        unsigned rowGlobal = var(elementNodes[m],i);
        for (unsigned n=0; n<nElementNodes; n++)
        {
          for (unsigned j=0; j<nVariables; j++)
          {
            unsigned colLocal = var(n,j);
            unsigned colGlobal = var(elementNodes[n],j);
            AMat->add(rowGlobal,colGlobal,elementTimeMat[rowLocal][colLocal]);
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::imposeEssentialBoundaryConditions ()
{
  for (unsigned b=0; b<nBoundaries; b++)
  {
    Boundary* boundary = boundaries[b];
    unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
    IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
    std::string type = boundary->getType();

    if (type == "Inlet")
    {
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
        }
      }

      /*for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          AMat->clearRow(var(node,nIons));
          bVec[var(node,nIons)] = 0.;
          for (unsigned i=0; i<nIons; i++)
          {
            double zi = (double)mitrem->getIonChargeNumber(i);
            AMat->linearCombination(var(node,nIons),var(node,i),zi);
            bVec[var(node,nIons)] += zi*bVec[var(node,i)];
          }
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
        }
      }*/
    }

    else if (type == "VirtualElectrode")
    {
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nVariables; i++)
          {
            bVec[var(node,i)] = inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }

/*          bVec[var(node,nIons)] = -0.0314583921944 - xVec[var(node,nIons)];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = inletVec[i] - xVec[var(node,i)];
          }
          for (unsigned i=0; i<nVariables; i++)
          {
            AMat->imposeVar(var(node,i));
          }*/
        }
      }
    }

    else if (type == "Concentration_0.5")
    {
      std::cout << "Impose concentration to 0.5 times bulk at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = 0.5*inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
          //bVec[var(node,nIons)] = inletVec[nIons] - xVec[var(node,nIons)];
          //AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Concentration_1.5")
    {
      std::cout << "Impose concentration to 1.5 times bulk at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = 1.5*inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
          //bVec[var(node,nIons)] = inletVec[nIons] - xVec[var(node,nIons)];
          //AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Potential_0")
    {
      std::cout << "Impose potential to 0V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          //for (unsigned i=0; i<nIons; i++)
          //{
          //  bVec[var(node,i)] = inletVec[i] - xVec[var(node,i)];
          //  AMat->imposeVar(var(node,i));
          //}
          bVec[var(node,nIons)] = 0. - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Potential_1")
    {
      std::cout << "Impose potential to 1V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          //for (unsigned i=0; i<nIons; i++)
          //{
          //  bVec[var(node,i)] = inletVec[i] - xVec[var(node,i)];
          //  AMat->imposeVar(var(node,i));
          //}
          bVec[var(node,nIons)] = 1. - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    if (type == "ImposeConcentrationsEVERYWHERE")
    {
      std::cout << "Impose bulk concentrations everywhere." << std::endl;
      for (unsigned m=0; m<nNodes; m++)
      {
        /*AMat->clearRow(var(m,nIons));
        bVec[var(m,nIons)] = 0.;
        for (unsigned i=0; i<nIons; i++)
        {
          double zi = (double)mitrem->getIonChargeNumber(i);
          AMat->linearCombination(var(m,nIons),var(m,i),zi);
          bVec[var(m,nIons)] += zi*bVec[var(m,i)];
        }*/
        for (unsigned i=0; i<nIons; i++)
        {
          bVec[var(m,i)] = inletVec[i] - xVec[var(m,i)];
          AMat->imposeVar(var(m,i));
        }
      }
    }

    if (type == "ImposePotentialEVERYWHERE")
    {
      for (unsigned m=0; m<nNodes; m++)
      {
        bVec[var(m,nIons)] = 0. - xVec[var(m,nIons)];
        AMat->imposeVar(var(m,nIons));
      }
    }

    else if (type == "Dirichlet_0.5_0")
    {
      std::cout << "Impose concentration to 0.5 times bulk and potential to 0V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = 0.5*inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
          bVec[var(node,nIons)] = 0. - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Dirichlet_1.5_0")
    {
      std::cout << "Impose concentration to 1.5 times bulk and potential to 0V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = 1.5*inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
          bVec[var(node,nIons)] = 0. - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Dirichlet_1.5_1")
    {
      std::cout << "Impose concentration to 1.5 times bulk and potential to 1V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          for (unsigned i=0; i<nIons; i++)
          {
            bVec[var(node,i)] = 1.5*inletVec[i] - xVec[var(node,i)];
            AMat->imposeVar(var(node,i));
          }
          bVec[var(node,nIons)] = 1. - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Dirichlet_50_first_-0.10031")
    {
      std::cout << "Impose concentration of FIRST ion to 50 times bulk and potential -0.10031V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          bVec[var(node,0)] = 25*inletVec[0] - xVec[var(node,0)];
          AMat->imposeVar(var(node,0));
          bVec[var(node,nIons)] = -0.088376740286 - xVec[var(node,nIons)];//-0.10031 - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Dirichlet_Johan_en_Steven_cathode")
    {
      std::cout << "Impose concentration of FIRST ion to 50 times bulk and potential -0.10031V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          bVec[var(node,0)] = 1e-16 - xVec[var(node,0)];
          AMat->imposeVar(var(node,0));
          bVec[var(node,1)] = 499.4995 - xVec[var(node,1)];
          AMat->imposeVar(var(node,1));
          bVec[var(node,2)] = 499.4995 - xVec[var(node,2)];
          AMat->imposeVar(var(node,2));
          bVec[var(node,nIons)] = 0.94  - xVec[var(node,nIons)];//-0.10031 - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

    else if (type == "Dirichlet_Johan_en_Steven_anode")
    {
      std::cout << "Impose concentration of FIRST ion to 50 times bulk and potential -0.10031V at boundary " << b << std::endl;
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          unsigned node = boundaryElementNodes[m];
          bVec[var(node,0)] = 1 - xVec[var(node,0)];
          AMat->imposeVar(var(node,0));
          bVec[var(node,1)] = 499 - xVec[var(node,1)];
          AMat->imposeVar(var(node,1));
          bVec[var(node,2)] = 500 - xVec[var(node,2)];
          AMat->imposeVar(var(node,2));
          bVec[var(node,nIons)] = 0.94 - xVec[var(node,nIons)];//-0.10031 - xVec[var(node,nIons)];
          AMat->imposeVar(var(node,nIons));
        }
      }
    }

  }
}
//---------------------------------------------------------------------------
void Reactor::calcCurrentDensities ()
{
  std::cout << "Calculating current densities..." << std::endl;
  for (unsigned e=0; e<nElements; e++)
  {
    Element* element = elements[e];
    //unsigned nElementNodes = element->getNNodes();
    IndexList elementNodes = element->getNodes();
    for (unsigned m=0; m<nElementNodes; m++)
    {
      coordinates[m] = nodes[elementNodes[m]]->getComponents();
      velocities[m] = flowField[elementNodes[m]]->getComponents();
      concentrations[m] = &xVec[var(elementNodes[m],0)];
      potentials[m] = xVec[var(elementNodes[m],nIons)];
    }
    DoubleVectorList ionCurrentDensity = elementMatrixAssembler->calcIonCurrentDensities(coordinates, velocities, concentrations, potentials, temperatures, densities, voidFractions, magneticFieldVectors);
    for (unsigned i=0; i<nVariables; i++)
    {
      for (unsigned c=0; c<nDimensions; c++)
      {
        ionCurrentDensities[e][i][c] = ionCurrentDensity[i][c];
      }
    }
  }

  unsigned nElecReactionsMax = mitrem->getNElecReactions();
  unsigned nGasReactionsMax = mitrem->getNGasReactions();
  for (unsigned e=0; e<nElectrodes; e++)
  {
    currents[e] = 0.;
    Electrode* electrode = electrodes[e];
    unsigned nBoundariesOfElectrode = electrode->getNBoundaries();
    IndexList boundariesOfElectrode = electrode->getBoundaries();
    for (unsigned b=0; b<nBoundariesOfElectrode; b++)
    {
      Boundary* boundary = boundaries[boundariesOfElectrode[b]];
      unsigned nBoundaryElementsOfBoundary = boundary->getNBoundaryElements();
      IndexList boundaryElementsOfBoundary = boundary->getBoundaryElements();
      for (unsigned be=0; be<nBoundaryElementsOfBoundary; be++)
      {
        Element* boundaryElement = boundaryElements[boundaryElementsOfBoundary[be]];
        //unsigned nBoundaryElementNodes = boundaryElement->getNNodes();
        IndexList boundaryElementNodes = boundaryElement->getNodes();
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          coordinates[m] = nodes[boundaryElementNodes[m]]->getComponents();
          concentrations[m] = &xVec[var(boundaryElementNodes[m],0)];
          potentials[m] = xVec[var(boundaryElementNodes[m],nIons)];
        }
        unsigned nElecReactions = electrode->getNElecReactions();
        IndexList elecReactions = electrode->getElecReactions();
        double electrodePotential = electrode->getPotential();
        unsigned nGasReactions = electrode->getNGasReactions();
        IndexList gasReactions = electrode->getGasReactions();
        DoubleListList elecReactionCurrentDensity = elementMatrixAssembler->calcElecReactionCurrentDensities(coordinates, concentrations, potentials, temperatures, densities, voidFractions, nElecReactions, elecReactions, electrodePotential);
        DoubleListList gasReactionRate = elementMatrixAssembler->calcGasReactionRates(coordinates, concentrations, potentials, temperatures, densities, voidFractions, nGasReactions, gasReactions);
        currents[e] += elementMatrixAssembler->calcCurrent(coordinates, concentrations, potentials, temperatures, densities, voidFractions, nElecReactions, elecReactions, electrodePotential);
        for (unsigned m=0; m<nBoundaryElementNodes; m++)
        {
          for (unsigned r=0; r<nElecReactionsMax; r++)
          {
            elecReactionCurrentDensities[e][b][be][m][r] = elecReactionCurrentDensity[m][r];
          }
          for (unsigned r=0; r<nGasReactionsMax; r++)
          {
            gasReactionRates[e][b][be][m][r] = gasReactionRate[m][r];
          }
        }
      }
    }
  }
}
//---------------------------------------------------------------------------
void Reactor::skip(std::ifstream &istrm, unsigned nSkips, char delimiter)
{
  for (unsigned i=0; i<nSkips; i++)
  {
    istrm.ignore(256,delimiter);
  }
}
//---------------------------------------------------------------------------
void Reactor::init()
{
  for (unsigned i=0; i<nIons; i++)
  {
    inletVec[i] = mitrem->getIonInletConcentration(i);
  }
  inletVec[nIons] = 0.;
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void Reactor::errorConflictingData (const std::string &file1, const std::string &file2, const std::string &conflict)
{
  char *file1CHAR = new char[file1.size()+1];
  char *file2CHAR = new char[file2.size()+1];
  char *conflictCHAR = new char[conflict.size()+1];
  strcpy(file1CHAR, file1.c_str());
  strcpy(file2CHAR, file2.c_str());
  strcpy(conflictCHAR, conflict.c_str());
  std::cout << "CONFLICTING DATA IN " << file1CHAR << " AND " << file2CHAR << ".\
      \nTHE " << conflictCHAR <<" ARE INCOMPATIBLE." << std::endl;
  system("pause");
  exit(1);
  delete[] file1CHAR;
  delete[] file2CHAR;
  delete[] conflictCHAR;
}
//---------------------------------------------------------------------------
void Reactor::errorFileDoesNotExist(const std::string &file)
{
  char *fileCHAR = new char[file.size()+1];
  strcpy(fileCHAR, file.c_str());
  std::cout << "ERROR IN DatFileReader.cpp.\
      \nTHE FILE " << fileCHAR <<" DOES NOT EXIST." << std::endl;
  system("pause");
  exit(1);
  delete[] fileCHAR;
}
//---------------------------------------------------------------------------

