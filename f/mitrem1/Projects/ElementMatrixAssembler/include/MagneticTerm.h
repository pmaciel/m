//---------------------------------------------------------------------------

#ifndef MagneticTermH
#define MagneticTermH

//---------------------------------------------------------------------------

#include "ElementTerm.h"

//---------------------------------------------------------------------------

class MagneticTerm : public ElementTerm
{
public :
	MagneticTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
	virtual ~MagneticTerm();

	virtual void	calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};
	virtual void	calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

	virtual void	calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

protected:
	double**		W;
	EmptyDoubleVector	vCrossB;
	double**		cCoefficients;
};

//---------------------------------------------------------------------------

#endif

