//---------------------------------------------------------------------------

#ifndef MagneticTerm_3D_GalerkinH
#define MagneticTerm_3D_GalerkinH

//---------------------------------------------------------------------------

#include "MagneticTerm.h"

//---------------------------------------------------------------------------

class MagneticTerm_3D_Galerkin : public MagneticTerm
{
public :
	MagneticTerm_3D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
	virtual ~MagneticTerm_3D_Galerkin();

	virtual void	calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
	virtual void	calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

