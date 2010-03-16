//---------------------------------------------------------------------------

#include "ElementProps_3D.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElementProps_3D::ElementProps_3D(unsigned nDimensions_) 
	 : ElementProps(nDimensions_)
{
}
//---------------------------------------------------------------------------
ElementProps_3D::~ElementProps_3D()
{
}
//---------------------------------------------------------------------------
double ElementProps_3D::calcSize(DoubleVectorList coordinates) const
{
  return -1./6.*(   coordinates[0][0]*coordinates[1][1]*coordinates[2][2] 
				  + coordinates[1][0]*coordinates[2][1]*coordinates[0][2]
				  + coordinates[2][0]*coordinates[0][1]*coordinates[1][2] 

				  + coordinates[3][0]*coordinates[1][1]*coordinates[0][2] 
				  + coordinates[1][0]*coordinates[0][1]*coordinates[3][2]
				  + coordinates[0][0]*coordinates[3][1]*coordinates[1][2] 

				  + coordinates[0][0]*coordinates[2][1]*coordinates[3][2]
				  + coordinates[2][0]*coordinates[3][1]*coordinates[0][2] 
				  + coordinates[3][0]*coordinates[0][1]*coordinates[2][2]

				  + coordinates[3][0]*coordinates[2][1]*coordinates[1][2]
				  + coordinates[2][0]*coordinates[1][1]*coordinates[3][2] 
				  + coordinates[1][0]*coordinates[3][1]*coordinates[2][2]
				  
				  - coordinates[2][0]*coordinates[1][1]*coordinates[0][2]
				  - coordinates[1][0]*coordinates[0][1]*coordinates[2][2]
				  - coordinates[0][0]*coordinates[2][1]*coordinates[1][2]

				  - coordinates[0][0]*coordinates[1][1]*coordinates[3][2]
				  - coordinates[1][0]*coordinates[3][1]*coordinates[0][2]
				  - coordinates[3][0]*coordinates[0][1]*coordinates[1][2]

				  - coordinates[3][0]*coordinates[2][1]*coordinates[0][2]
				  - coordinates[2][0]*coordinates[0][1]*coordinates[3][2]
				  - coordinates[0][0]*coordinates[3][1]*coordinates[2][2]

				  - coordinates[1][0]*coordinates[2][1]*coordinates[3][2]
				  - coordinates[2][0]*coordinates[3][1]*coordinates[1][2]
				  - coordinates[3][0]*coordinates[1][1]*coordinates[2][2]
				  );
}
//---------------------------------------------------------------------------
DoubleVector ElementProps_3D::calcNormal(unsigned m, DoubleVectorList coordinates) const
{
	double sign = -1. * ( 0.5-1.*(m%2) );
	unsigned p = (m+1)%4;
	unsigned q = (m+2)%4;
	unsigned r = (m+3)%4;
	normal[0] = sign*(coordinates[p][1]*coordinates[q][2] - coordinates[q][1]*coordinates[p][2] + coordinates[q][1]*coordinates[r][2] - coordinates[r][1]*coordinates[q][2] + coordinates[r][1]*coordinates[p][2] - coordinates[p][1]*coordinates[r][2]);
	normal[1] = sign*(coordinates[p][2]*coordinates[q][0] - coordinates[q][2]*coordinates[p][0] + coordinates[q][2]*coordinates[r][0] - coordinates[r][2]*coordinates[q][0] + coordinates[r][2]*coordinates[p][0] - coordinates[p][2]*coordinates[r][0]);
	normal[2] = sign*(coordinates[p][0]*coordinates[q][1] - coordinates[q][0]*coordinates[p][1] + coordinates[q][0]*coordinates[r][1] - coordinates[r][0]*coordinates[q][1] + coordinates[r][0]*coordinates[p][1] - coordinates[p][0]*coordinates[r][1]);
	return normal;
}
//---------------------------------------------------------------------------

