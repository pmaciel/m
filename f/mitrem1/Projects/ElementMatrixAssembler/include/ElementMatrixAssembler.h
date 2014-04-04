//---------------------------------------------------------------------------

#ifndef ElementMatrixAssemblerH
#define ElementMatrixAssemblerH

//---------------------------------------------------------------------------

#include "TypeDefs.h"

#include "MITReM.h"

#include "ElementProps.h"
#include "BoundaryElementProps.h"

#include "ConvectionTerm.h"
#include "DiffusionTerm.h"
#include "MigrationTerm.h"
#include "MagneticTerm.h"
#include "ElectrostaticsTerm.h"
#include "HomReactionTerm.h"
#include "TimeTerm.h"

#include "ElecReactionTerm.h"

#include "ElementContribution.h"
#include "BoundaryElementContribution.h"

//---------------------------------------------------------------------------

class ElementMatrixAssembler
{
public :
  ElementMatrixAssembler(
    const std::string &dimensions,
    MITReM* mitrem_,
    const std::string &convectionScheme,
    const std::string &diffusionScheme,
    const std::string &migrationScheme,
    const std::string &magneticScheme,
    const std::string &homReactionScheme,
    const std::string &electrostaticsScheme,
    const std::string &timeScheme,
    const std::string &elecReactionScheme,
    const std::string &gasReactionScheme,
    const bool _is_bubble=false,
    const bool _charge_flux=true,
    const bool _swap_first_and_last_equations=true);
  ~ElementMatrixAssembler();

  DoubleMatrix     calcElementMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  DoubleMatrix     calcElementJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

  DoubleVector     calcBoundaryElementVec(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);
  DoubleMatrix     calcBoundaryElementJac(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions);

  DoubleMatrix     calcElementTimeMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
  DoubleMatrix     calcElementTimeJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

  DoubleVectorList calcIonCurrentDensities(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

  DoubleListList   calcElecReactionCurrentDensities(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
  DoubleListList   calcGasReactionRates(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions);

  double           calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential);
  double           calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions);

  double calcESize(DoubleVectorList coordinates) const;
  double calcBESize(DoubleVectorList coordinates) const;

private:
  MITReM* mitrem;

  ElementProps*                elementProps;
  BoundaryElementProps*        boundaryElementProps;
  ElementContribution*         elementContribution;
  BoundaryElementContribution* boundaryElementContribution;

  ConvectionTerm*     convectionTerm;
  DiffusionTerm*      diffusionTerm;
  MigrationTerm*      migrationTerm;
  MagneticTerm*       magneticTerm;
  HomReactionTerm*    homReactionTerm;
  ElectrostaticsTerm* electrostaticsTerm;
  TimeTerm*           timeTerm;
  ElecReactionTerm*   elecReactionTerm;
  GasReactionTerm*    gasReactionTerm;

   unsigned nDimensions, nElementNodes, nBoundaryElementNodes, nVariables, nIons, nElecReactions, nGasReactions;

  // if charge conservation is to be assembled in place of one mass balance
  bool m_chargeconservation;

  // if first and last equations are to be swapped
  bool m_swap;

  void errorInvalidScheme(const std::string &scheme);
};

//---------------------------------------------------------------------------

inline DoubleVectorList ElementMatrixAssembler::calcIonCurrentDensities(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  return elementContribution->calcIonCurrentDensities(coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
}

//---------------------------------------------------------------------------

inline DoubleListList ElementMatrixAssembler::calcElecReactionCurrentDensities(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions_, IndexList elecReactions, double electrodePotential)
{
  return boundaryElementContribution->calcElecReactionCurrentDensities(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions_, elecReactions, electrodePotential);
}

//---------------------------------------------------------------------------

inline DoubleListList ElementMatrixAssembler::calcGasReactionRates(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions_, IndexList gasReactions)
{
  return boundaryElementContribution->calcGasReactionRates(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nGasReactions_, gasReactions);
}

//---------------------------------------------------------------------------

inline double ElementMatrixAssembler::calcCurrent(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions_, IndexList elecReactions, double electrodePotential)
{
  return boundaryElementContribution->calcCurrent(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions_, elecReactions, electrodePotential);
}

//---------------------------------------------------------------------------

inline double ElementMatrixAssembler::calcGasGeneration(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions_, IndexList gasReactions)
{
  return boundaryElementContribution->calcGasGeneration(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nGasReactions_, gasReactions);
}

//---------------------------------------------------------------------------

inline double ElementMatrixAssembler::calcESize(DoubleVectorList coordinates) const
{
  return elementProps->calcSize(coordinates);
}

//---------------------------------------------------------------------------

inline double ElementMatrixAssembler::calcBESize(DoubleVectorList coordinates) const
{
  return boundaryElementProps->calcSize(coordinates);
}

//---------------------------------------------------------------------------

#endif

