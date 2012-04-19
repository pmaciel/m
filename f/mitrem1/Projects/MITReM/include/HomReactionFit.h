//---------------------------------------------------------------------------

#ifndef HomReactionFitH
#define HomReactionFitH

//---------------------------------------------------------------------------

#include <string>
#include <vector>
#include "Electrolyte.h"

//---------------------------------------------------------------------------

class HomReactionFit
{
public :
  HomReactionFit(const std::string &label_, Electrolyte** electrolytes, unsigned nElectrolytes);
  ~HomReactionFit();

  // Set methods
  void    setEquilibriumConstant(double K);
  void    setMinimumEquilibriumConstant(double KMin);
  void    setMaximumEquilibriumConstant(double KMax);

  // Get methods
  std::string  getLabel() const;
  double    getEquilibriumConstant() const;
  double    getMinimumEquilibriumConstant() const;
  double    getMaximumEquilibriumConstant() const;
  unsigned  getElectrolyte() const;

  // Methods
  double    calcOsmoticCoefficientDerivativeEquilibriumConstant(unsigned t);
  double    calcActivityCoefficientDerivativeEquilibriumConstant(unsigned t);
  double    calcConductivityDerivativeEquilibriumConstant(unsigned t);
  double    calcTransportNumberDerivativeEquilibriumConstant(unsigned t);
  double    calcDiffusionCoefficientDerivativeEquilibriumConstant(unsigned t);

private :
  class HomReactionReference
  {
  public :
    Electrolyte*  electrolyte;
    unsigned    electrolyteIndex;
    unsigned    reactionIndex;
  };

  std::string  label;
  double    K,KMin,KMax,epsilon;
  HomReactionReference homReactionReference;

  void errorHomReactionNotFound(const std::string label);
};

//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
inline void HomReactionFit::setMinimumEquilibriumConstant(double KMin_)
{
  this->KMin = KMin_;
}
//---------------------------------------------------------------------------
inline void HomReactionFit::setMaximumEquilibriumConstant(double KMax_)
{
  this->KMax = KMax_;
}
//---------------------------------------------------------------------------
inline void HomReactionFit::setEquilibriumConstant(double K_)
{
  this->K = K_;
  homReactionReference.electrolyte->setHomReactionEquilibriumConstant(homReactionReference.reactionIndex,K_);
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
inline std::string HomReactionFit::getLabel() const
{
  return label;
}
//---------------------------------------------------------------------------
inline double HomReactionFit::getEquilibriumConstant() const
{
  return K;
}
//---------------------------------------------------------------------------
inline double HomReactionFit::getMinimumEquilibriumConstant() const
{
  return KMin;
}
//---------------------------------------------------------------------------
inline double HomReactionFit::getMaximumEquilibriumConstant() const
{
  return KMax;
}
//---------------------------------------------------------------------------
inline unsigned HomReactionFit::getElectrolyte() const
{
  return homReactionReference.electrolyteIndex;
}
//---------------------------------------------------------------------------

#endif

