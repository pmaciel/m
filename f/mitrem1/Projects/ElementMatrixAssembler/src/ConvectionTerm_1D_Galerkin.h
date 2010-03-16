//---------------------------------------------------------------------------

#ifndef ConvectionTerm_1D_GalerkinH
#define ConvectionTerm_1D_GalerkinH

//---------------------------------------------------------------------------

#include "ConvectionTerm.h"

//---------------------------------------------------------------------------

class ConvectionTerm_1D_Galerkin : public ConvectionTerm
{
public :
	ConvectionTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
	virtual ~ConvectionTerm_1D_Galerkin();

	virtual void	calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions);
	virtual void	calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions);

protected:
	double**		Nu;
};

//---------------------------------------------------------------------------

#endif

