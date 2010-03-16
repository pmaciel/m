//---------------------------------------------------------------------------

#include "Electrode.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Electrode::Electrode (unsigned nBoundaries_, unsigned nElecReactions_, unsigned nGasReactions_)
	: nBoundaries(nBoundaries_), nElecReactions(nElecReactions_), nGasReactions(nGasReactions_)
{	
	boundaries = new unsigned[nBoundaries];
	elecReactions = new unsigned[nElecReactions];
	gasReactions = new unsigned[nGasReactions];
}
//---------------------------------------------------------------------------
Electrode::~Electrode ()
{
	delete[] boundaries;
	delete[] elecReactions;
	delete[] gasReactions;
}
//---------------------------------------------------------------------------
void Electrode::setBoundaries(unsigned b, unsigned boundary)
{
	if (b < nBoundaries)
	{
		boundaries[b] = boundary;
	}
	else
	{
		// error message
	}
}
//---------------------------------------------------------------------------
void Electrode::setElecReactions(unsigned r, unsigned elecReaction)
{
	if (r < nElecReactions)
	{
		elecReactions[r] = elecReaction;
	}
	else
	{
		// error message
	}
}
//---------------------------------------------------------------------------
void Electrode::setGasReactions(unsigned r, unsigned gasReaction)
{
	if (r < nGasReactions)
	{
		gasReactions[r] = gasReaction;
	}
	else
	{
		// error message
	}
}
//---------------------------------------------------------------------------