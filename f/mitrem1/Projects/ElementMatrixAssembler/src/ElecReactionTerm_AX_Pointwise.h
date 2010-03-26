//---------------------------------------------------------------------------

#ifndef ElecReactionTerm_AX_PointwiseH
#define ElecReactionTerm_AX_PointwiseH

//---------------------------------------------------------------------------

#include "ElecReactionTerm.h"

//---------------------------------------------------------------------------

class ElecReactionTerm_AX_Pointwise : public ElecReactionTerm
{
public :
  ElecReactionTerm_AX_Pointwise(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_);
  virtual ~ElecReactionTerm_AX_Pointwise();

  virtual void  calcVec(EmptyDoubleVector boundaryElementVec, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions, double factor);
  virtual void  calcJac(EmptyDoubleMatrix boundaryElementJac, DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions, double factor);

  virtual double  calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
};

//---------------------------------------------------------------------------

#endif

