//---------------------------------------------------------------------------

#ifndef ElectrolyteSolutionH
#define ElectrolyteSolutionH

//---------------------------------------------------------------------------

#include <string>
#include "xmlParser.h"

//---------------------------------------------------------------------------

class ElectrolyteSolution
{
public :
  ElectrolyteSolution(unsigned nIons/*const std::string &name*/);
  ~ElectrolyteSolution();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solution
  void    setSolutionDensity(double solutionDensity);
  void    setSolutionKinematicViscosity(double solutionKinematicViscosity);
  void    setSolutionTemperature(double solutionTemperature);
  void    setSolutionPotential(double solutionPotential);
  double    getSolutionDensity() const;
  double    getSolutionKinematicViscosity() const;
  double    getSolutionTemperature() const;
  double    getSolutionPotential() const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solvent
  void    setSolventDielectricConstant(double solventDielectricConst);
  void    setSolventDynamicViscosity(double solventDynamicViscosity);
  double    getSolventDielectricConstant() const;
  double    getSolventDynamicViscosity() const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Ions
  void    setIonLabel(unsigned i, const std::string &label);
  void    setIonChargeNumber(unsigned i, int chargeNumber);
  void    setIonDiffusionConstant(unsigned i, double diffusionConstant);
  void    setIonDiffusionConstantLimit(unsigned i, double diffusionConstantLimit);
  void    setIonDiffusionConstantPower(unsigned i, double diffusionConstantPower);
  void    setIonDiameter(unsigned i, double diameter);
  void    setIonMolarMass(unsigned i, double molarMass);
  void    setIonConcentration(unsigned i, double concentration);
  void    setIonInletConcentration(unsigned i, double inletConcentration);
  inline void setIonTVExpansionCoefficient(unsigned i, double _alpha)    {  ions[i].alpha = _alpha;  }
  inline void setIonCDensificationCoefficient(unsigned i, double _beta)  {  ions[i].beta  = _beta;   }
  inline void setIonMMagneticSusceptibility(unsigned i, double _mmchi)   {  ions[i].mmchi = _mmchi;  }
  unsigned  getNIons() const;
  std::string  getIonLabel(unsigned i) const;
  int      getIonChargeNumber(unsigned i) const;
  double    getIonDiffusionConstant(unsigned i) const;
  double    getIonDiffusionConstantLimit(unsigned i) const;
  double    getIonDiffusionConstantPower(unsigned i) const;
  double    getIonDiameter(unsigned i) const;
  double    getIonMolarMass(unsigned i) const;
  double    getIonConcentration(unsigned i) const;
  double    getIonInletConcentration(unsigned i) const;
  inline double getIonTVExpansionCoefficient(unsigned i) const     {  return ions[i].alpha;  }
  inline double getIonCDensificationCoefficient(unsigned i) const  {  return ions[i].beta;   }
  inline double getIonMMagneticSusceptibility(unsigned i) const    {  return ions[i].mmchi;  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Methods
  //void    swapIons(unsigned i, unsigned j);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

private :
  class Ion
  {
  public :
    std::string  label;
    int      z;            // charge number
    double    M;            // molar mass [kg/mol]
    double    d;            // diameter [m]
    double    D;            // diffusion constant [m2/s]
    double    DLim;          // diffusion constant limit [m2/s]
    double    DPow;          // diffusion constant power [m3/mol]
    double    cInlet;          // inlet concentration [mol/m3]
    double    c;            // concentration [mol/m3]
    double    alpha;          // temperature densification coefficient (alpha) [1/K]
    double    beta;          // concentration-wise densification coefficient (beta) [1/mol]
    double    mmchi;          // molar magnetic susceptibility [m3/mol]
  };

  // Members
  double    solutionKinematicViscosity;  // solution kinematic viscosity [m2/s]
  double    solutionDensity;      // solution density [kg/m3]
  double    solutionTemperature;    // solution temperature [K]
  double    solutionPotential;      // solution potential [V]
  double    solventDielectricConst;    // solvent dielectric constant [C/Vm]  ( = 4*pi*eps0*epsr)
  double    solventDynamicViscosity;  // solvent dynamic viscosity [kg/ms]
  unsigned  nIons;
  Ion*    ions;
  //std::string electrolyteSolutionFile;

  // Read methods
  //void    readElectrolyteSolution();
};

//---------------------------------------------------------------------------


//--- SOLUTION --------------------------------------------------------------
inline void ElectrolyteSolution::setSolutionDensity(double solutionDensity)
{
  this->solutionDensity = solutionDensity;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setSolutionKinematicViscosity(double solutionKinematicViscosity)
{
  this->solutionKinematicViscosity = solutionKinematicViscosity;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setSolutionTemperature(double solutionTemperature)
{
  this->solutionTemperature = solutionTemperature;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setSolutionPotential(double solutionPotential)
{
  this->solutionPotential = solutionPotential;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getSolutionDensity() const
{
  return solutionDensity;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getSolutionKinematicViscosity() const
{
  return solutionKinematicViscosity;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getSolutionTemperature()const
{
  return solutionTemperature;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getSolutionPotential() const
{
  return solutionPotential;
}
//---------------------------------------------------------------------------


//--- SOLVENT ----------------------------------------------------------------
inline void ElectrolyteSolution::setSolventDielectricConstant(double solventDielectricConst)
{
  this->solventDielectricConst = solventDielectricConst;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setSolventDynamicViscosity(double solventDynamicViscosity)
{
  this->solventDynamicViscosity = solventDynamicViscosity;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getSolventDielectricConstant() const
{
  return solventDielectricConst;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getSolventDynamicViscosity() const
{
  return solventDynamicViscosity;
}
//---------------------------------------------------------------------------


//--- IONS ------------------------------------------------------------------
inline void ElectrolyteSolution::setIonLabel(unsigned i, const std::string &label)
{
  this->ions[i].label = label;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonChargeNumber(unsigned i, int chargeNumber)
{
  this->ions[i].z = chargeNumber;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonDiffusionConstant(unsigned i, double diffusionConstant)
{
  this->ions[i].D = diffusionConstant;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonDiffusionConstantLimit(unsigned i, double diffusionConstantLimit)
{
  this->ions[i].DLim = diffusionConstantLimit;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonDiffusionConstantPower(unsigned i, double diffusionConstantPower)
{
  this->ions[i].DPow = diffusionConstantPower;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonDiameter(unsigned i, double diameter)
{
  this->ions[i].d = diameter;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonMolarMass(unsigned i, double molarMass)
{
  this->ions[i].M = molarMass;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonConcentration(unsigned i, double concentration)
{
  this->ions[i].c = concentration;
}
//---------------------------------------------------------------------------
inline void ElectrolyteSolution::setIonInletConcentration(unsigned i, double inletConcentration)
{
  this->ions[i].cInlet = inletConcentration;
}
//---------------------------------------------------------------------------
inline unsigned ElectrolyteSolution::getNIons() const
{
  return nIons;
}
//---------------------------------------------------------------------------
inline std::string ElectrolyteSolution::getIonLabel(unsigned i) const
{
  return ions[i].label;
}
//---------------------------------------------------------------------------
inline int ElectrolyteSolution::getIonChargeNumber(unsigned i) const
{
  return ions[i].z;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonDiffusionConstant(unsigned i) const
{
  return ions[i].D;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonDiffusionConstantLimit(unsigned i) const
{
  return ions[i].DLim;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonDiffusionConstantPower(unsigned i) const
{
  return ions[i].DPow;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonDiameter(unsigned i) const
{
  return ions[i].d;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonMolarMass(unsigned i) const
{
  return ions[i].M;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonConcentration(unsigned i) const
{
  return ions[i].c;
}
//---------------------------------------------------------------------------
inline double ElectrolyteSolution::getIonInletConcentration(unsigned i) const
{
  return ions[i].cInlet;
}
//---------------------------------------------------------------------------

#endif

