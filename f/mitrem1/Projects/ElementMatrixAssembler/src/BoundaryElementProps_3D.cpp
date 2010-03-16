//---------------------------------------------------------------------------

#include "BoundaryElementProps_3D.h"
#include <math.h>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
BoundaryElementProps_3D::BoundaryElementProps_3D(unsigned nDimensions_) 
	 : BoundaryElementProps(nDimensions_)
{
}
//---------------------------------------------------------------------------
BoundaryElementProps_3D::~BoundaryElementProps_3D()
{
}
//---------------------------------------------------------------------------
double BoundaryElementProps_3D::calcSize(DoubleVectorList coordinates) const
{
	double dxy =  coordinates[0][0]*coordinates[1][1]
				+ coordinates[1][0]*coordinates[2][1]
				+ coordinates[2][0]*coordinates[0][1]
				- coordinates[0][0]*coordinates[2][1] 
				- coordinates[1][0]*coordinates[0][1]
				- coordinates[2][0]*coordinates[1][1];
	double dyz =  coordinates[0][1]*coordinates[1][2]
				+ coordinates[1][1]*coordinates[2][2]
				+ coordinates[2][1]*coordinates[0][2]
				- coordinates[0][1]*coordinates[2][2] 
				- coordinates[1][1]*coordinates[0][2]
				- coordinates[2][1]*coordinates[1][2];
	double dzx =  coordinates[0][2]*coordinates[1][0]
				+ coordinates[1][2]*coordinates[2][0]
				+ coordinates[2][2]*coordinates[0][0]
				- coordinates[0][2]*coordinates[2][0] 
				- coordinates[1][2]*coordinates[0][0]
				- coordinates[2][2]*coordinates[1][0];
	return 0.5*sqrt(dxy*dxy + dyz*dyz + dzx*dzx);
}
//---------------------------------------------------------------------------
