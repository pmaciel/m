//---------------------------------------------------------------------------

#ifndef ElecReaction_BVadsH
#define ElecReaction_BVadsH

//---------------------------------------------------------------------------

#include "ElecReaction.h"

//---------------------------------------------------------------------------

class ElecReaction_BVads : public ElecReaction
{
public :
  ElecReaction_BVads(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_);
  virtual ~ElecReaction_BVads();

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

