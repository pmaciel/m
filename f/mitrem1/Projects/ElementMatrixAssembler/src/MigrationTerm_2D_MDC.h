//---------------------------------------------------------------------------

#ifndef MigrationTerm_2D_MDCH
#define MigrationTerm_2D_MDCH

//---------------------------------------------------------------------------

#include "MigrationTerm.h"

//---------------------------------------------------------------------------

class MigrationTerm_2D_MDC : public MigrationTerm
{
public :
  MigrationTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~MigrationTerm_2D_MDC();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

