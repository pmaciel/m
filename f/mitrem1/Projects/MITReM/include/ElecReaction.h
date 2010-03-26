//---------------------------------------------------------------------------

#ifndef ElecReactionH
#define ElecReactionH

//---------------------------------------------------------------------------

#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------

class ElecReaction
{
public :
  ElecReaction(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_);
  virtual ~ElecReaction();

  // Set methods
  void    setLabel(const std::string &label);
  void    setKinParam(unsigned index, double value);
  void    setNElectrons(unsigned nElectrons);
  void    setAgentsRed(unsigned i, unsigned agent);
  void    setStoichRed(unsigned i, int stoich);
  void    setOrderRed(unsigned i, double order);
  void    setAgentsOxi(unsigned i, unsigned agent);
  void    setStoichOxi(unsigned i, int stoich);
  void    setOrderOxi(unsigned i, double order);

  // Get methods
  unsigned  getNAgentsRed() const;
  unsigned  getNAgentsOxi() const;
  std::string getLabel() const;
  double    getKinParam(unsigned index) const;
  unsigned  getNElectrons() const;
  unsigned  getAgentsRed(unsigned i) const;
  int      getStoichRed(unsigned i) const;
  double    getOrderRed(unsigned i) const;
  unsigned  getAgentsOxi(unsigned i) const;
  int      getStoichOxi(unsigned i) const;
  double    getOrderOxi(unsigned i) const;
  int      getStoichOf(unsigned i) const;

  // Methods
  virtual double  calcReactionRate(double V) const = 0;
  double      calcReactionCurrentDensity(double V) const;
  virtual double  calcReactionRateDerivativeU(double V) const = 0;
  virtual double  calcReactionRateDerivativeV(double V) const = 0;
  virtual double  calcReactionRateDerivativeCRed(double V, unsigned i) const = 0;
  virtual double  calcReactionRateDerivativeCOxi(double V, unsigned i) const = 0;
  virtual double  calcEquilibriumPotential() const = 0;

protected :
  // Members
  std::string  label;
  double*    kineticParameters;
  unsigned  nElectrons;             // number of electrons exchanged
  unsigned  nAgentsRed,nAgentsOxi;
  unsigned*  agentsRed;
  unsigned*  agentsOxi;
  int*    stoichRed;
  int*    stoichOxi;
  double*    orderRed;
  double*    orderOxi;
  ElectrolyteSolution*    electrolyteSolution;
};

//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
inline void ElecReaction::setLabel (const std::string &label)
{
  this->label = label;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setKinParam (unsigned index, double value)
{
  kineticParameters[index] = value;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setNElectrons (unsigned nElectrons)
{
  this->nElectrons = nElectrons;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setAgentsRed (unsigned i, unsigned agent)
{
  agentsRed[i] = agent;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setStoichRed (unsigned i, int stoich)
{
  stoichRed[i] = stoich;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setOrderRed (unsigned i, double order)
{
  orderRed[i] = order;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setAgentsOxi (unsigned i, unsigned agent)
{
  agentsOxi[i] = agent;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setStoichOxi (unsigned i, int stoich)
{
  stoichOxi[i] = stoich;
}
//---------------------------------------------------------------------------
inline void ElecReaction::setOrderOxi (unsigned i, double order)
{
  orderOxi[i] = order;
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
inline unsigned ElecReaction::getNAgentsRed() const
{
  return nAgentsRed;
}
//---------------------------------------------------------------------------
inline unsigned ElecReaction::getNAgentsOxi() const
{
  return nAgentsOxi;
}
//---------------------------------------------------------------------------
inline std::string ElecReaction::getLabel() const
{
  return label;
}
//---------------------------------------------------------------------------
inline double ElecReaction::getKinParam(unsigned index) const
{
  return kineticParameters[index];
}
//---------------------------------------------------------------------------
inline unsigned ElecReaction::getNElectrons() const
{
  return nElectrons;
}
//---------------------------------------------------------------------------
inline unsigned ElecReaction::getAgentsRed(unsigned i) const
{
  return agentsRed[i];
}
//---------------------------------------------------------------------------
inline int ElecReaction::getStoichRed(unsigned i) const
{
  return stoichRed[i];
}
//---------------------------------------------------------------------------
inline double ElecReaction::getOrderRed(unsigned i) const
{
  return orderRed[i];
}
//---------------------------------------------------------------------------
inline unsigned ElecReaction::getAgentsOxi(unsigned i) const
{
  return agentsOxi[i];
}
//---------------------------------------------------------------------------
inline int ElecReaction::getStoichOxi(unsigned i) const
{
  return stoichOxi[i];
}
//---------------------------------------------------------------------------
inline double ElecReaction::getOrderOxi(unsigned i) const
{
  return orderOxi[i];
}
//---------------------------------------------------------------------------

#endif

