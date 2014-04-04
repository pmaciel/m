//---------------------------------------------------------------------------

#ifndef ElecReaction_BVH
#define ElecReaction_BVH

//---------------------------------------------------------------------------

#include "ElecReaction.h"

//---------------------------------------------------------------------------

class ElecReaction_BV : public ElecReaction
{
public :
  ElecReaction_BV(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_);
  virtual ~ElecReaction_BV();

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

