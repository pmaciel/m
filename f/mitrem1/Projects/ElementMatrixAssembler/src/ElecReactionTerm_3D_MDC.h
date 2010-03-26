//---------------------------------------------------------------------------

#ifndef ElecReactionTerm_3D_MDCH
#define ElecReactionTerm_3D_MDCH

//---------------------------------------------------------------------------

#include "ElecReactionTerm.h"

//---------------------------------------------------------------------------

class ElecReactionTerm_3D_MDC : public ElecReactionTerm
{
public :
  ElecReactionTerm_3D_MDC(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_);
  virtual ~ElecReactionTerm_3D_MDC();

  virtual void  calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
  virtual void  calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
};

//---------------------------------------------------------------------------

#endif

