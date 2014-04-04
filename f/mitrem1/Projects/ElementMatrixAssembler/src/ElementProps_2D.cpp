//---------------------------------------------------------------------------

#include "ElementProps_2D.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElementProps_2D::ElementProps_2D(unsigned nDimensions_)
   : ElementProps(nDimensions_)
{
}
//---------------------------------------------------------------------------
ElementProps_2D::~ElementProps_2D()
{
}
//---------------------------------------------------------------------------
double ElementProps_2D::calcSize(DoubleVectorList coordinates) const
{
  return 0.5*(  coordinates[0][0]*coordinates[1][1]
        + coordinates[1][0]*coordinates[2][1]
        + coordinates[2][0]*coordinates[0][1]
        - coordinates[0][0]*coordinates[2][1]
        - coordinates[1][0]*coordinates[0][1]
        - coordinates[2][0]*coordinates[1][1]);
}
//---------------------------------------------------------------------------
DoubleVector ElementProps_2D::calcNormal(unsigned m, DoubleVectorList coordinates) const
{
  unsigned p = (m+1)%3;
  unsigned q = (m+2)%3;
  normal[0] = coordinates[p][1] - coordinates[q][1];
  normal[1] = coordinates[q][0] - coordinates[p][0];
  return normal;
}
//---------------------------------------------------------------------------

