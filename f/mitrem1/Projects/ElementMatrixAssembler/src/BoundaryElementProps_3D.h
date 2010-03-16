//---------------------------------------------------------------------------

#ifndef BoundaryElementProps_3DH
#define BoundaryElementProps_3DH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "BoundaryElementProps.h"

//---------------------------------------------------------------------------

class BoundaryElementProps_3D : public BoundaryElementProps
{
public :
	BoundaryElementProps_3D(unsigned nDimensions_);
	virtual ~BoundaryElementProps_3D();

	virtual double	calcSize(DoubleVectorList coordinates) const;
};

//---------------------------------------------------------------------------

#endif

