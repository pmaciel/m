//---------------------------------------------------------------------------

#ifndef TimeTerm_1D_HCUH
#define TimeTerm_1D_HCUH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "TimeTerm.h"
#include "ElementProps.h"

//---------------------------------------------------------------------------

class TimeTerm_1D_HCU : public TimeTerm
{
public :
  TimeTerm_1D_HCU(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
  virtual ~TimeTerm_1D_HCU();

  virtual void  calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions);
  virtual void  calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions);

private:
  double*      alpha;
};

//---------------------------------------------------------------------------

#endif

