//---------------------------------------------------------------------------

#ifndef ElectrostaticsTerm_AX_MDCH
#define ElectrostaticsTerm_AX_MDCH

//---------------------------------------------------------------------------

#include "ElectrostaticsTerm.h"

//---------------------------------------------------------------------------

class ElectrostaticsTerm_AX_MDC : public ElectrostaticsTerm
{
public :
  ElectrostaticsTerm_AX_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ElectrostaticsTerm_AX_MDC();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

