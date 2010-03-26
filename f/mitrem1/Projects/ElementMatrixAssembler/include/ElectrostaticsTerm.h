//---------------------------------------------------------------------------

#ifndef ElectrostaticsTermH
#define ElectrostaticsTermH

//---------------------------------------------------------------------------

#include "ElementTerm.h"

//---------------------------------------------------------------------------

class ElectrostaticsTerm : public ElementTerm
{
public :
  ElectrostaticsTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ElectrostaticsTerm();

  virtual void calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};
  virtual void calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

  virtual void calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

protected:
  double  K;
  double* Z;
};

//---------------------------------------------------------------------------

#endif

