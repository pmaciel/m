//---------------------------------------------------------------------------

#ifndef ElementProps_2DH
#define ElementProps_2DH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "ElementProps.h"

//---------------------------------------------------------------------------

class ElementProps_2D : public ElementProps
{
public :
  ElementProps_2D(unsigned nDimensions_);
  virtual ~ElementProps_2D();

  virtual double  calcSize(DoubleVectorList coordinates) const;
  virtual DoubleVector  calcNormal(unsigned m, DoubleVectorList coordinates) const;

  // ElementProps nodes are numbered counterclockwise!!!
};

//---------------------------------------------------------------------------

#endif

