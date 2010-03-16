//---------------------------------------------------------------------------

#ifndef ElementPropsH
#define ElementPropsH

//---------------------------------------------------------------------------

#include "TypeDefs.h"

//---------------------------------------------------------------------------

class ElementProps
{
public :
	ElementProps(unsigned nDimensions_);
	virtual ~ElementProps();

	virtual double	calcSize(DoubleVectorList coordinates) const = 0;
	virtual DoubleVector	calcNormal(unsigned m, DoubleVectorList coordinates) const = 0; // internal normal

protected :
	unsigned		nDimensions;
	EmptyDoubleVector			normal;
};

//---------------------------------------------------------------------------

#endif

