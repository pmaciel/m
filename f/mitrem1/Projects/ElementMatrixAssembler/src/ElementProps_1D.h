//---------------------------------------------------------------------------

#ifndef ElementProps_1DH
#define ElementProps_1DH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "ElementProps.h"

//---------------------------------------------------------------------------

class ElementProps_1D : public ElementProps
{
public :
  ElementProps_1D(unsigned nDimensions_);
  virtual ~ElementProps_1D();

  virtual double  calcSize(DoubleVectorList coordinates) const;
  virtual DoubleVector  calcNormal(unsigned m, DoubleVectorList coordinates) const;

  // ElementProps nodes are numbered such that x1 > x0!!!
};

//---------------------------------------------------------------------------

#endif

