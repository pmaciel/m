//---------------------------------------------------------------------------

#ifndef MigrationTerm_1D_ISGH
#define MigrationTerm_1D_ISGH

//---------------------------------------------------------------------------

#include "MigrationTerm.h"

//---------------------------------------------------------------------------

class MigrationTerm_1D_ISG : public MigrationTerm
{
public :
  MigrationTerm_1D_ISG(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~MigrationTerm_1D_ISG();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

