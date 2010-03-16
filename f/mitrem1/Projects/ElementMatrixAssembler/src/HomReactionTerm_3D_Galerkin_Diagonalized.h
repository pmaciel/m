//---------------------------------------------------------------------------

#ifndef HomReactionTerm_3D_Galerkin_DiagonalizedH
#define HomReactionTerm_3D_Galerkin_DiagonalizedH

//---------------------------------------------------------------------------

#include "HomReactionTerm.h"

//---------------------------------------------------------------------------

class HomReactionTerm_3D_Galerkin_Diagonalized : public HomReactionTerm
{
public :
	HomReactionTerm_3D_Galerkin_Diagonalized(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
	virtual ~HomReactionTerm_3D_Galerkin_Diagonalized();

	virtual void	calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
	virtual void	calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

