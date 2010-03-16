//---------------------------------------------------------------------------

#ifndef BoundaryElementProps_2DH
#define BoundaryElementProps_2DH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "BoundaryElementProps.h"

//---------------------------------------------------------------------------

class BoundaryElementProps_2D : public BoundaryElementProps
{
public :
	BoundaryElementProps_2D(unsigned nDimensions_);
	virtual ~BoundaryElementProps_2D();

	virtual double	calcSize(DoubleVectorList coordinates) const;
};

//---------------------------------------------------------------------------

#endif

