//---------------------------------------------------------------------------

#ifndef ElectrostaticsTerm_1D_MDCH
#define ElectrostaticsTerm_1D_MDCH

//---------------------------------------------------------------------------

#include "ElectrostaticsTerm.h"

//---------------------------------------------------------------------------

class ElectrostaticsTerm_1D_MDC : public ElectrostaticsTerm
{
public :
  ElectrostaticsTerm_1D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ElectrostaticsTerm_1D_MDC();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

