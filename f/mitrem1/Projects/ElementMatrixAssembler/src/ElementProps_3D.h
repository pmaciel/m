//---------------------------------------------------------------------------

#ifndef ElementProps_3DH
#define ElementProps_3DH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "ElementProps.h"

//---------------------------------------------------------------------------

class ElementProps_3D : public ElementProps
{
public :
	ElementProps_3D(unsigned nDimensions_);
	virtual ~ElementProps_3D();

	virtual double	calcSize(DoubleVectorList coordinates) const;
	virtual DoubleVector	calcNormal(unsigned m, DoubleVectorList coordinates) const;

	// ElementProps nodes are numbered (how?) !!!
};

//---------------------------------------------------------------------------

#endif

