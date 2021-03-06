//---------------------------------------------------------------------------

#ifndef HomReactionTerm_2D_MDCH
#define HomReactionTerm_2D_MDCH

//---------------------------------------------------------------------------

#include "HomReactionTerm.h"

//---------------------------------------------------------------------------

class HomReactionTerm_2D_MDC : public HomReactionTerm
{
public :
  HomReactionTerm_2D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~HomReactionTerm_2D_MDC();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
};

//---------------------------------------------------------------------------

#endif

