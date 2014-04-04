//---------------------------------------------------------------------------

#ifndef ConvectionTermH
#define ConvectionTermH

//---------------------------------------------------------------------------

#include "ElementTerm.h"

//---------------------------------------------------------------------------

class ConvectionTerm : public ElementTerm
{
public :
  ConvectionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ConvectionTerm();

  virtual void calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};
  virtual void calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

  virtual void calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

protected:
  EmptyDoubleVector Vm;
  EmptyDoubleVector k;
  EmptyDoubleVector averageVelocity;
};

//---------------------------------------------------------------------------

#endif

