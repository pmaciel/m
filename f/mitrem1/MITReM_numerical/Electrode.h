//---------------------------------------------------------------------------

#ifndef ElectrodeH
#define ElectrodeH

//---------------------------------------------------------------------------

#include <string>
#include "TypeDefs.h"

//---------------------------------------------------------------------------

class Electrode
{
public :
	Electrode(unsigned nBoundaries_, unsigned nElecReactions_, unsigned nGasReactions_);
	~Electrode();

	void		setBoundaries(unsigned b, unsigned boundary);
	void		setElecReactions(unsigned r, unsigned elecReaction);
	void		setGasReactions(unsigned r, unsigned gasReaction);
	void		setLabel(const std::string &label);
	void		setPotential(double potential);
	
	unsigned	getNBoundaries() const;
	unsigned	getNElecReactions() const;
	unsigned	getNGasReactions() const;
	IndexList	getBoundaries() const;
	unsigned	getBoundaries(unsigned b) const;
	IndexList	getElecReactions() const;
	IndexList	getGasReactions() const;
	unsigned	getElecReactions(unsigned r) const;
	unsigned	getGasReactions(unsigned r) const;
	std::string	getLabel() const;
	double		getPotential() const;

private :
	unsigned		nBoundaries;
	EmptyIndexList	boundaries;				// array of boundary indices
	unsigned		nElecReactions;
	unsigned		nGasReactions;
	EmptyIndexList	elecReactions;			// array of electrode reactions indices
	EmptyIndexList	gasReactions;			// array of electrode reactions indices
	std::string		label;					// type of boundary
	double			potential;				// electrode potential
};

//---------------------------------------------------------------------------

inline void Electrode::setLabel(const std::string &label)
{
	this->label = label;
}

//---------------------------------------------------------------------------

inline void Electrode::setPotential(double potential)
{
	this->potential = potential;
}

//---------------------------------------------------------------------------

inline unsigned Electrode::getNBoundaries() const
{
	return nBoundaries;
}

//---------------------------------------------------------------------------

inline unsigned Electrode::getNElecReactions() const
{
	return nElecReactions;
}

//---------------------------------------------------------------------------

inline unsigned Electrode::getNGasReactions() const
{
	return nGasReactions;
}

//---------------------------------------------------------------------------

inline IndexList Electrode::getBoundaries() const
{
	return boundaries;
}

//---------------------------------------------------------------------------

inline unsigned Electrode::getBoundaries(unsigned b) const
{
	return boundaries[b];
}

//---------------------------------------------------------------------------

inline IndexList Electrode::getElecReactions() const
{
	return elecReactions;
}

//---------------------------------------------------------------------------

inline IndexList Electrode::getGasReactions() const
{
	return gasReactions;
}

//---------------------------------------------------------------------------

inline unsigned Electrode::getElecReactions(unsigned r) const
{
	return elecReactions[r];
}

//---------------------------------------------------------------------------

inline unsigned Electrode::getGasReactions(unsigned r) const
{
	return gasReactions[r];
}

//---------------------------------------------------------------------------

inline std::string Electrode::getLabel() const
{
	return label;
}

//---------------------------------------------------------------------------

inline double Electrode::getPotential() const
{
	return potential;
}

//---------------------------------------------------------------------------

#endif