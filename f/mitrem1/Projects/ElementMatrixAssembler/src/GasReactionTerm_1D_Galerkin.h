//---------------------------------------------------------------------------

#ifndef GasReactionTerm_1D_GalerkinH
#define GasReactionTerm_1D_GalerkinH

//---------------------------------------------------------------------------

#include "GasReactionTerm.h"

//---------------------------------------------------------------------------

class GasReactionTerm_1D_Galerkin : public GasReactionTerm
{
public :
  GasReactionTerm_1D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_, const bool isBubble_);
  virtual ~GasReactionTerm_1D_Galerkin();

  virtual void calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
  virtual void calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
};

//---------------------------------------------------------------------------

#endif

