//---------------------------------------------------------------------------

#include "BoundaryElementProps_1D.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
BoundaryElementProps_1D::BoundaryElementProps_1D(unsigned nDimensions_) 
	 : BoundaryElementProps(nDimensions_)
{
}
//---------------------------------------------------------------------------
BoundaryElementProps_1D::~BoundaryElementProps_1D()
{
}
//---------------------------------------------------------------------------
double BoundaryElementProps_1D::calcSize(DoubleVectorList coordinates) const
{
	return 1.;
}
//---------------------------------------------------------------------------

