//---------------------------------------------------------------------------

#ifndef MagneticTerm_2D_GalerkinH
#define MagneticTerm_2D_GalerkinH

//---------------------------------------------------------------------------

#include "MagneticTerm.h"

//---------------------------------------------------------------------------

class MagneticTerm_2D_Galerkin : public MagneticTerm
{
public :
  MagneticTerm_2D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~MagneticTerm_2D_Galerkin();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

