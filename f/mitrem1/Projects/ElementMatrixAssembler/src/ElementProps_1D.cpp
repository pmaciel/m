//---------------------------------------------------------------------------

#include "ElementProps_1D.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElementProps_1D::ElementProps_1D(unsigned nDimensions_)
   : ElementProps(nDimensions_)
{
}
//---------------------------------------------------------------------------
ElementProps_1D::~ElementProps_1D()
{
}
//---------------------------------------------------------------------------
double ElementProps_1D::calcSize(DoubleVectorList coordinates) const
{
  return coordinates[1][0] - coordinates[0][0];
}
//---------------------------------------------------------------------------
DoubleVector ElementProps_1D::calcNormal(unsigned m, DoubleVectorList coordinates) const
{
  // Confusing in 1D: the normal of node 1 is a vector on the side lying opposite of node 1, which is node 2!!!
  if (m == 0) normal[0] = -1;
  else normal[0] = 1;
  return normal;
}
//---------------------------------------------------------------------------

