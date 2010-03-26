//---------------------------------------------------------------------------

#ifndef TimeTermH
#define TimeTermH

//---------------------------------------------------------------------------

#include "ElementTerm.h"

//---------------------------------------------------------------------------

class TimeTerm : public ElementTerm
{
public :
  TimeTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~TimeTerm();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

  virtual void  calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};
};

//---------------------------------------------------------------------------

#endif

