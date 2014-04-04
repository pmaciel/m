//---------------------------------------------------------------------------

#ifndef ElectrostaticsTerm_1D_GalerkinH
#define ElectrostaticsTerm_1D_GalerkinH

//---------------------------------------------------------------------------

#include "ElectrostaticsTerm.h"

//---------------------------------------------------------------------------

class ElectrostaticsTerm_1D_Galerkin : public ElectrostaticsTerm
{
public :
  ElectrostaticsTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ElectrostaticsTerm_1D_Galerkin();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

