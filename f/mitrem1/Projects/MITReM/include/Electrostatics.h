//---------------------------------------------------------------------------

#ifndef ElectrostaticsH
#define ElectrostaticsH

//---------------------------------------------------------------------------

#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------

class Electrostatics
{
public :
  Electrostatics(ElectrolyteSolution* electrolyteSolution_) : electrolyteSolution(electrolyteSolution_) { };
  virtual ~Electrostatics() { };

  virtual double  calcConcentrationFactor(unsigned j) const = 0;
  virtual double  calcPotentialFactor() const = 0;

protected :
  ElectrolyteSolution*  electrolyteSolution;
};

//---------------------------------------------------------------------------

class Electrostatics_Electroneutrality : public Electrostatics
{
public :
  Electrostatics_Electroneutrality(ElectrolyteSolution* electrolyteSolution_) : Electrostatics(electrolyteSolution_) { };
  virtual ~Electrostatics_Electroneutrality() { };

  virtual double  calcConcentrationFactor(unsigned j) const;
  virtual double  calcPotentialFactor() const;
};

//---------------------------------------------------------------------------

class Electrostatics_Poisson : public Electrostatics
{
public :
  Electrostatics_Poisson(ElectrolyteSolution* electrolyteSolution_) : Electrostatics(electrolyteSolution_) { };
  virtual ~Electrostatics_Poisson() { };

  virtual double  calcConcentrationFactor(unsigned j) const;
  virtual double  calcPotentialFactor() const;
};

//---------------------------------------------------------------------------

#endif

