//---------------------------------------------------------------------------

#include "ElectrolyteSolution.h"
#include "DatFileReader.h"
#include <sstream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrolyteSolution::ElectrolyteSolution(unsigned nIons/*const std::string &name*/)
{
  // /*
  // Required formats:

  // *******************************
  // .electrolytesolution:
  // Versie 1.0
  // [solventDielectricConst] = ...
  // [solventDynamicViscosity] = ...
  // [solutionDensity] = ...
  // [solutionKinematicViscosity] = ...
  // [temperature] = ...
  // [nIons] = ...
  //   <label> = ...  <z> = ...  <D> = ...  <d> = ...  <M> = ...  <cInlet> = ...
  // *******************************
  // */

  /*electrolyteSolutionFile = name + ".electrolytesolution";

  readElectrolyteSolution();*/
  this->nIons = nIons;
  ions = new Ion[nIons];
}
//---------------------------------------------------------------------------
ElectrolyteSolution::~ElectrolyteSolution()
{
  delete[] ions;
}
//---------------------------------------------------------------------------

//--- READ METHODS ----------------------------------------------------------
/*void ElectrolyteSolution::readElectrolyteSolution()
{
  double doubleValue;
  std::string stringValue;
  int intValue;

  DatFileReader datFile(electrolyteSolutionFile);

  datFile.readScalar("[solutionKinematicViscosity]",doubleValue);
  setSolutionKinematicViscosity(doubleValue);
  datFile.readScalar("[solutionDensity]",doubleValue);
  setSolutionDensity(doubleValue);
  datFile.readScalar("[solutionTemperature]",doubleValue);
  setSolutionTemperature(doubleValue);
  setSolutionPotential(0.);
  datFile.readScalar("[solventDielectricConstant]",doubleValue);
  setSolventDielectricConstant(doubleValue);
  datFile.readScalar("[solventDynamicViscosity]",doubleValue);
  setSolventDynamicViscosity(doubleValue);

  nIons = datFile.readMultipleVector_nVectors("[nIons]");
  ions = new Ion[nIons];
  for (unsigned i=0; i<nIons; i++)
  {
    datFile.readMultipleVector("[nIons]", i, "<label>", stringValue);
    setIonLabel(i,stringValue);
    datFile.readMultipleVector("[nIons]", i, "<z>", intValue);
    setIonChargeNumber(i,intValue);
    datFile.readMultipleVector("[nIons]", i, "<D>", doubleValue);
    setIonDiffusionConstant(i,doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<d>", doubleValue);
    setIonDiameter(i,doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<M>", doubleValue);
    setIonMolarMass(i,doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<DLim>", doubleValue);
    setIonDiffusionConstantLimit(i,doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<DPow>", doubleValue);
    setIonDiffusionConstantPower(i,doubleValue);
    datFile.readMultipleVector("[nIons]", i, "<cInlet>", doubleValue);
    //setIonInletConcentration(i,doubleValue);
    setIonConcentration(i,doubleValue);
//    setIonTVExpansionCoefficient(i,s.getAttribute< double >("alpha",0.));
//    setIonCDensificationCoefficient(i,s.getAttribute< double >("beta",0.));
//    setIonMMagneticSusceptibility(i,s.getAttribute< double >("MMChi",0.));
  }
}*/
//---------------------------------------------------------------------------


//--- METHODS ---------------------------------------------------------------
/*void ElectrolyteSolution::swapIons(unsigned i, unsigned j)
{
  Ion ionTemp = ions[i];
  ions[i] = ions[j];
  ions[j] = ionTemp;
}*/
//---------------------------------------------------------------------------

