//---------------------------------------------------------------------------

#ifndef BoundaryElementProps_1DH
#define BoundaryElementProps_1DH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "BoundaryElementProps.h"

//---------------------------------------------------------------------------

class BoundaryElementProps_1D : public BoundaryElementProps
{
public :
	BoundaryElementProps_1D(unsigned nDimensions_);
	virtual ~BoundaryElementProps_1D();

	virtual double	calcSize(DoubleVectorList coordinates) const;
};

//---------------------------------------------------------------------------

#endif

