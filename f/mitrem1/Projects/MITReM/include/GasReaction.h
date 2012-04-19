//---------------------------------------------------------------------------

#ifndef GasReactionH
#define GasReactionH

//---------------------------------------------------------------------------

#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------

class GasReaction
{
public :
  GasReaction(ElectrolyteSolution* electrolyteSolution_);
  virtual ~GasReaction();

  // Set methods
  void    setLabel(const std::string &label);
  void    setKinParam(unsigned index, double value);
  void    setDissolvedGas(unsigned dissolvedGas);

  // Get methods
  std::string getLabel() const;
  double    getKinParam(unsigned index) const;
  unsigned  getDissolvedGas() const;

  // Methods
  double  calcReactionRate() const;
  double  calcReactionRateDerivativeCDissGas() const;

protected :
  // Members
  std::string  label;
  double*    kineticParameters;
  unsigned  dissolvedGas;
  ElectrolyteSolution*    electrolyteSolution;
};

//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
inline void GasReaction::setLabel (const std::string &label_)
{
  this->label = label_;
}
//---------------------------------------------------------------------------
inline void GasReaction::setKinParam (unsigned index, double value)
{
  kineticParameters[index] = value;
}
//---------------------------------------------------------------------------
inline void GasReaction::setDissolvedGas(unsigned dissolvedGas_)
{
  this->dissolvedGas = dissolvedGas_;
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
inline std::string GasReaction::getLabel() const
{
  return label;
}
//---------------------------------------------------------------------------
inline double GasReaction::getKinParam(unsigned index) const
{
  return kineticParameters[index];
}
//---------------------------------------------------------------------------
inline unsigned GasReaction::getDissolvedGas() const
{
  return dissolvedGas;
}
//---------------------------------------------------------------------------

#endif

