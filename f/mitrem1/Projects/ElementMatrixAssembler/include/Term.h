//---------------------------------------------------------------------------

#ifndef TermH
#define TermH

//---------------------------------------------------------------------------

#include "MITReM.h"
#include "MathematicsPhysics.h"
#include "TypeDefs.h"

//---------------------------------------------------------------------------

class Term
{
public :
	Term(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_);
	virtual ~Term();

	unsigned		var(unsigned n, unsigned j) const;
	unsigned		eq(unsigned m, unsigned i) const;

protected :
	unsigned		nDimensions;
	unsigned		nNodes;
	unsigned		nVariables;
	unsigned		nIons;
	double			elementSize;
	DoubleVector	normal;
	EmptyEmptyDoubleVectorList		normals;
	MITReM*			mitrem;
};

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
inline unsigned Term::var(unsigned n, unsigned j) const
{
	return n*nVariables+j;
}
//---------------------------------------------------------------------------
inline unsigned Term::eq(unsigned m, unsigned i) const
{
	//if (i == 0) return m*nVariables+nIons;
	//else if (i == nIons) return m*nVariables;
	//else return m*nVariables+i;
	return m*nVariables+i;
}
//---------------------------------------------------------------------------

#endif

