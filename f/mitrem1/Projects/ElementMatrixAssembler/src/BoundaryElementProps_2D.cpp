//---------------------------------------------------------------------------

#include "BoundaryElementProps_2D.h"
#include <math.h>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
BoundaryElementProps_2D::BoundaryElementProps_2D(unsigned nDimensions_)
   : BoundaryElementProps(nDimensions_)
{
}
//---------------------------------------------------------------------------
BoundaryElementProps_2D::~BoundaryElementProps_2D()
{
}
//---------------------------------------------------------------------------
double BoundaryElementProps_2D::calcSize(DoubleVectorList coordinates) const
{
  double dx = coordinates[1][0] - coordinates[0][0];
  double dy = coordinates[1][1] - coordinates[0][1];
  return sqrt(dx*dx + dy*dy);
}
//---------------------------------------------------------------------------
