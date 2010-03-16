//---------------------------------------------------------------------------

#define U				electrolyteSolution->getSolutionPotential()
#define aaa				kineticParameters[0]
#define bbb				kineticParameters[1]

//---------------------------------------------------------------------------

#include "ElecReaction_Linear.h"
#include <math.h>
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElecReaction_Linear::ElecReaction_Linear(ElectrolyteSolution* electrolyteSolution_, unsigned nAgentsRed_, unsigned nAgentsOxi_) 
: ElecReaction(electrolyteSolution_, nAgentsRed_, nAgentsOxi_)
{
	kineticParameters = new double[2];
}
//---------------------------------------------------------------------------
ElecReaction_Linear::~ElecReaction_Linear()
{
	delete[] kineticParameters;
}
//---------------------------------------------------------------------------
double ElecReaction_Linear::calcReactionRate(double V) const
{
	return aaa*(V-U) + bbb;
}
//---------------------------------------------------------------------------
double ElecReaction_Linear::calcReactionRateDerivativeU(double V) const
{
	return -aaa;
}
//---------------------------------------------------------------------------
double ElecReaction_Linear::calcReactionRateDerivativeV(double V) const
{
	return aaa;
}
//---------------------------------------------------------------------------
double ElecReaction_Linear::calcReactionRateDerivativeCRed(double V, unsigned i) const
{
	return 0.;
}
//---------------------------------------------------------------------------
double ElecReaction_Linear::calcReactionRateDerivativeCOxi(double V, unsigned i) const
{
	return 0.;
}
//---------------------------------------------------------------------------
double ElecReaction_Linear::calcEquilibriumPotential() const
{
	return -bbb/aaa;
}
//---------------------------------------------------------------------------


