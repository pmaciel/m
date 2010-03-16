//---------------------------------------------------------------------------

#ifndef ElementH
#define ElementH

//---------------------------------------------------------------------------

#include "TypeDefs.h"

//---------------------------------------------------------------------------

class Element
{
public :
	Element(unsigned nNodes_);
	~Element();

	void		setNodes(unsigned m, unsigned value);

	unsigned	getNNodes() const;
	IndexList	getNodes() const;
	unsigned	getNodes(unsigned m) const;

private :
	unsigned		nNodes;
	EmptyIndexList	nodes;				// array of node indices
};

//---------------------------------------------------------------------------

inline unsigned Element::getNNodes() const
{
	return nNodes;
}

//---------------------------------------------------------------------------

inline IndexList Element::getNodes() const
{
	return nodes;
}

//---------------------------------------------------------------------------

inline unsigned Element::getNodes(unsigned m) const
{
	return nodes[m];
}

//---------------------------------------------------------------------------

#endif