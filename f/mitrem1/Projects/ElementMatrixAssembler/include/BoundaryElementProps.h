//---------------------------------------------------------------------------

#ifndef BoundaryElementPropsH
#define BoundaryElementPropsH

//---------------------------------------------------------------------------

#include "TypeDefs.h"

//---------------------------------------------------------------------------

class BoundaryElementProps
{
public :
	BoundaryElementProps(unsigned nDimensions_);
	virtual ~BoundaryElementProps();

	virtual double	calcSize(DoubleVectorList coordinates) const = 0;

protected :
	unsigned		nDimensions;
};

//---------------------------------------------------------------------------

#endif

