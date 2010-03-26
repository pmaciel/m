//---------------------------------------------------------------------------

#ifndef ElecReaction_LinearH
#define ElecReaction_LinearH

//---------------------------------------------------------------------------

#include "ElecReaction.h"

//---------------------------------------------------------------------------

class ElecReaction_Linear : public ElecReaction
{
public :
  ElecReaction_Linear(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_);
  virtual ~ElecReaction_Linear();

  // Methods
  virtual double  calcReactionRate(double V) const;
  virtual double  calcReactionRateDerivativeU(double V) const;
  virtual double  calcReactionRateDerivativeV(double V) const;
  virtual double  calcReactionRateDerivativeCRed(double V, unsigned i) const;
  virtual double  calcReactionRateDerivativeCOxi(double V, unsigned i) const;
  virtual double  calcEquilibriumPotential() const;
};

//---------------------------------------------------------------------------

#endif

