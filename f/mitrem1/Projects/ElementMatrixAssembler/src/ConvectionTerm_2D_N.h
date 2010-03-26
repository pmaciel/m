//---------------------------------------------------------------------------

#ifndef ConvectionTerm_2D_NH
#define ConvectionTerm_2D_NH

//---------------------------------------------------------------------------

#include "ConvectionTerm.h"

//---------------------------------------------------------------------------

class ConvectionTerm_2D_N : public ConvectionTerm
{
public :
  ConvectionTerm_2D_N(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ConvectionTerm_2D_N();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

