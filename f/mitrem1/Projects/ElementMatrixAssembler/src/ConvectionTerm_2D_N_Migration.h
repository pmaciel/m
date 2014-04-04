//---------------------------------------------------------------------------

#ifndef ConvectionTerm_2D_N_MigrationH
#define ConvectionTerm_2D_N_MigrationH

//---------------------------------------------------------------------------

#include "ConvectionTerm.h"

//---------------------------------------------------------------------------

class ConvectionTerm_2D_N_Migration : public ConvectionTerm
{
public :
  ConvectionTerm_2D_N_Migration(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~ConvectionTerm_2D_N_Migration();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

protected:
  double** gradc;
  double* gradU;
};

//---------------------------------------------------------------------------

#endif

