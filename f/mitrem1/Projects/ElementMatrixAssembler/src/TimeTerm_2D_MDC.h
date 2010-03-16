//---------------------------------------------------------------------------

#ifndef TimeTerm_2D_MDCH
#define TimeTerm_2D_MDCH

//---------------------------------------------------------------------------

#include "TimeTerm.h"

//---------------------------------------------------------------------------

class TimeTerm_2D_MDC : public TimeTerm
{
public :
	TimeTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
	virtual ~TimeTerm_2D_MDC();

	virtual void	calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
	virtual void	calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

