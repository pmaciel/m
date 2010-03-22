//---------------------------------------------------------------------------

#include <cstdlib>

#include "ElementMatrixAssembler.h"

#include "ElementProps_1D.h"
#include "BoundaryElementProps_1D.h"
#include "ConvectionTerm_1D_N.h"
#include "DiffusionTerm_1D_Galerkin.h"
#include "MigrationTerm_1D_Galerkin.h"
#include "HomReactionTerm_1D_Galerkin.h"
#include "HomReactionTerm_1D_Galerkin_Diagonalized.h"
#include "HomReactionTerm_1D_Galerkin_Diagonalized_M.h"
#include "ElectrostaticsTerm_1D_Galerkin.h"
#include "TimeTerm_1D_Galerkin.h"
#include "ElecReactionTerm_1D_Galerkin.h"
#include "GasReactionTerm_1D_Galerkin.h"

#include "DiffusionTerm_1D_MDC.h"
#include "MigrationTerm_1D_MDC.h"
#include "HomReactionTerm_1D_MDC.h"
#include "ElectrostaticsTerm_1D_MDC.h"
#include "TimeTerm_1D_MDC.h"
#include "ElecReactionTerm_1D_MDC.h"

#include "MigrationTerm_1D_ISG.h"

#include "ElementProps_2D.h"
#include "BoundaryElementProps_2D.h"
#include "ConvectionTerm_2D_N.h"
#include "ConvectionTerm_2D_N_Migration.h"
#include "DiffusionTerm_2D_Galerkin.h"
#include "MigrationTerm_2D_Galerkin.h"
#include "MigrationTerm_2D_N.h"
#include "MagneticTerm_2D_Galerkin.h"
#include "HomReactionTerm_2D_Galerkin.h"
#include "HomReactionTerm_2D_Galerkin_Diagonalized.h"
#include "HomReactionTerm_2D_Galerkin_Diagonalized_M.h"
#include "ElectrostaticsTerm_2D_Galerkin.h"
#include "TimeTerm_2D_Galerkin.h"
#include "ElecReactionTerm_2D_Galerkin.h"
#include "ElecReactionTerm_2D_Pointwise.h"
#include "GasReactionTerm_2D_Galerkin.h"

#include "DiffusionTerm_2D_MDC.h"
#include "MigrationTerm_2D_MDC.h"
#include "HomReactionTerm_2D_MDC.h"
#include "ElectrostaticsTerm_2D_MDC.h"
#include "TimeTerm_2D_MDC.h"
#include "ElecReactionTerm_2D_MDC.h"

#include "ConvectionTerm_AX_N.h"
#include "DiffusionTerm_AX_Galerkin.h"
#include "MigrationTerm_AX_Galerkin.h"
#include "HomReactionTerm_AX_Galerkin.h"
#include "HomReactionTerm_AX_Galerkin_Diagonalized.h"
#include "HomReactionTerm_AX_Galerkin_Diagonalized_M.h"
#include "ElectrostaticsTerm_AX_Galerkin.h"
#include "TimeTerm_AX_Galerkin.h"
#include "ElecReactionTerm_AX_Galerkin.h"
#include "GasReactionTerm_AX_Galerkin.h"

#include "DiffusionTerm_AX_MDC.h"
#include "MigrationTerm_AX_MDC.h"
#include "HomReactionTerm_AX_MDC.h"
#include "ElectrostaticsTerm_AX_MDC.h"
#include "TimeTerm_AX_MDC.h"
#include "ElecReactionTerm_AX_MDC.h"
#include "ElecReactionTerm_AX_Pointwise.h"

#include "ElementProps_3D.h"
#include "BoundaryElementProps_3D.h"
#include "ConvectionTerm_3D_N.h"
#include "DiffusionTerm_3D_Galerkin.h"
#include "MigrationTerm_3D_Galerkin.h"
#include "MagneticTerm_3D_Galerkin.h"
#include "HomReactionTerm_3D_Galerkin.h"
#include "HomReactionTerm_3D_Galerkin_Diagonalized.h"
#include "HomReactionTerm_3D_Galerkin_Diagonalized_M.h"
#include "ElectrostaticsTerm_3D_Galerkin.h"
#include "TimeTerm_3D_Galerkin.h"
#include "ElecReactionTerm_3D_Galerkin.h"
#include "ElecReactionTerm_3D_Pointwise.h"
#include "GasReactionTerm_3D_Galerkin.h"

#include "DiffusionTerm_3D_MDC.h"
#include "MigrationTerm_3D_MDC.h"
#include "HomReactionTerm_3D_MDC.h"
#include "ElectrostaticsTerm_3D_MDC.h"
#include "TimeTerm_3D_MDC.h"
#include "ElecReactionTerm_3D_MDC.h"

#include "ConvectionTerm_Template.h"

#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElementMatrixAssembler::ElementMatrixAssembler(
  const std::string &dimensions, MITReM* mitrem_,
  const std::string &convectionScheme,
  const std::string &diffusionScheme,
  const std::string &migrationScheme,
  const std::string &magneticScheme,
  const std::string &homReactionScheme,
  const std::string &electrostaticsScheme,
  const std::string &timeScheme,
  const std::string &elecReactionScheme,
  const std::string &gasReactionScheme,
  const bool _is_bubble,
  const bool _charge_conservation,
  const bool _swap_first_and_last_equations ) :
    mitrem(mitrem_)
{
  //std::string mitremFile_ = mitremFile;

  //mitrem = new MITReM(mitremFile_);

  nIons = mitrem->getNIons();
  nVariables = nIons+1;
  nElecReactions = mitrem->getNElecReactions();
  nGasReactions = mitrem->getNGasReactions();

  if (dimensions == "1D")
  {
    nDimensions = 1;
    nElementNodes = nDimensions+1;
    nBoundaryElementNodes = nDimensions;

    elementProps = new ElementProps_1D(nDimensions);
    boundaryElementProps = new BoundaryElementProps_1D(nDimensions);

    if (convectionScheme == "Empty")
    {
      convectionTerm = new ConvectionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    //else if (convectionScheme == "Galerkin")
    //{
    //	convectionTerm = new ConvectionTerm_1D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    //}
    //else if (convectionScheme == "MDC")
    //{
    //	convectionTerm = new ConvectionTerm_1D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    //}
    else if (convectionScheme == "N")
    {
      convectionTerm = new ConvectionTerm_1D_N(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("convectionScheme: "+convectionScheme);
    }

    if (diffusionScheme == "Empty")
    {
      diffusionTerm = new DiffusionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "Galerkin")
    {
      diffusionTerm = new DiffusionTerm_1D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "MDC")
    {
      diffusionTerm = new DiffusionTerm_1D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("diffusionScheme: "+diffusionScheme);
    }

    if (migrationScheme == "Empty")
    {
      migrationTerm = new MigrationTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "Galerkin")
    {
      migrationTerm = new MigrationTerm_1D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "MDC")
    {
      migrationTerm = new MigrationTerm_1D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "ISG")
    {
      migrationTerm = new MigrationTerm_1D_ISG(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("migrationScheme: "+migrationScheme);
    }

    if (magneticScheme == "Empty")
    {
      magneticTerm = new MagneticTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("magneticScheme: "+magneticScheme);
    }

    if (homReactionScheme == "Empty")
    {
      homReactionTerm = new HomReactionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin")
    {
      homReactionTerm = new HomReactionTerm_1D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal")
    {
      homReactionTerm = new HomReactionTerm_1D_Galerkin_Diagonalized(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal_M")
    {
      homReactionTerm = new HomReactionTerm_1D_Galerkin_Diagonalized_M(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "MDC")
    {
      homReactionTerm = new HomReactionTerm_1D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("homReactionScheme: "+homReactionScheme);
    }

    if (electrostaticsScheme == "Empty")
    {
      electrostaticsTerm = new ElectrostaticsTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "Galerkin")
    {
      electrostaticsTerm = new ElectrostaticsTerm_1D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "MDC")
    {
      electrostaticsTerm = new ElectrostaticsTerm_1D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("electrostaticsScheme: "+electrostaticsScheme);
    }

    if (timeScheme == "Empty")
    {
      timeTerm = new TimeTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "Galerkin")
    {
      timeTerm = new TimeTerm_1D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "MDC")
    {
      timeTerm = new TimeTerm_1D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("timeScheme: "+timeScheme);
    }

    if (elecReactionScheme == "Empty")
    {
      elecReactionTerm = new ElecReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Galerkin")
    {
      elecReactionTerm = new ElecReactionTerm_1D_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "MDC")
    {
      elecReactionTerm = new ElecReactionTerm_1D_MDC(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else
    {
      errorInvalidScheme("elecReactionScheme: "+elecReactionScheme);
    }

    if (gasReactionScheme == "Empty")
    {
      gasReactionTerm = new GasReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else if (gasReactionScheme == "Galerkin")
    {
      gasReactionTerm = new GasReactionTerm_1D_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else
    {
      errorInvalidScheme("gasReactionScheme: "+gasReactionScheme);
    }
  }
  else if (dimensions == "2D")
  {
    nDimensions = 2;
    nElementNodes = nDimensions+1;
    nBoundaryElementNodes = nDimensions;

    elementProps = new ElementProps_2D(nDimensions);
    boundaryElementProps = new BoundaryElementProps_2D(nDimensions);

    if (convectionScheme == "Empty")
    {
      convectionTerm = new ConvectionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (convectionScheme == "Galerkin")
    {
    //	convectionTerm = new ConvectionTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
      convectionTerm = new ConvectionTerm_Template<2,Scheme_Galerkin>(nVariables,mitrem,elementProps);
    }
    else if (convectionScheme == "LDA")
    {
      convectionTerm = new ConvectionTerm_Template<2,Scheme_LDA>(nVariables,mitrem,elementProps);
    }
    //else if (convectionScheme == "MDC")
    //{
    //	convectionTerm = new ConvectionTerm_2D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    //}
    else if (convectionScheme == "N" && migrationScheme == "N")
    {
      convectionTerm = new ConvectionTerm_2D_N_Migration(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (convectionScheme == "N")
    {
      convectionTerm = new ConvectionTerm_2D_N(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
//			convectionTerm = new ConvectionTerm_Template<2,Scheme_N>(nVariables,mitrem,elementProps);
    }
    else
    {
      errorInvalidScheme("convectionScheme: "+convectionScheme);
    }

    if (diffusionScheme == "Empty")
    {
      diffusionTerm = new DiffusionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "Galerkin")
    {
      diffusionTerm = new DiffusionTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "MDC")
    {
      diffusionTerm = new DiffusionTerm_2D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("diffusionScheme: "+diffusionScheme);
    }

    if (migrationScheme == "Empty")
    {
      migrationTerm = new MigrationTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "Galerkin")
    {
      migrationTerm = new MigrationTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "MDC")
    {
      migrationTerm = new MigrationTerm_2D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "N")
    {
      migrationTerm = new MigrationTerm_2D_N(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("migrationScheme: "+migrationScheme);
    }

    if (magneticScheme == "Empty")
    {
      magneticTerm = new MagneticTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
        else if (migrationScheme == "Galerkin")
    {
      magneticTerm = new MagneticTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("magneticScheme: "+magneticScheme);
    }

    if (homReactionScheme == "Empty")
    {
      homReactionTerm = new HomReactionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin")
    {
      homReactionTerm = new HomReactionTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal")
    {
      homReactionTerm = new HomReactionTerm_2D_Galerkin_Diagonalized(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal_M")
    {
      homReactionTerm = new HomReactionTerm_2D_Galerkin_Diagonalized_M(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "MDC")
    {
      homReactionTerm = new HomReactionTerm_2D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("homReactionScheme: "+homReactionScheme);
    }

    if (electrostaticsScheme == "Empty")
    {
      electrostaticsTerm = new ElectrostaticsTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "Galerkin")
    {
      electrostaticsTerm = new ElectrostaticsTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "MDC")
    {
      electrostaticsTerm = new ElectrostaticsTerm_2D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("electrostaticsScheme: "+electrostaticsScheme);
    }

    if (timeScheme == "Empty")
    {
      timeTerm = new TimeTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "Galerkin")
    {
      timeTerm = new TimeTerm_2D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "MDC")
    {
      timeTerm = new TimeTerm_2D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("timeScheme: "+timeScheme);
    }

    if (elecReactionScheme == "Empty")
    {
      elecReactionTerm = new ElecReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Galerkin")
    {
      elecReactionTerm = new ElecReactionTerm_2D_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "MDC")
    {
      elecReactionTerm = new ElecReactionTerm_2D_MDC(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Pointwise")
    {
      elecReactionTerm = new ElecReactionTerm_2D_Pointwise(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else
    {
      errorInvalidScheme("elecReactionScheme: "+elecReactionScheme);
    }

    if (gasReactionScheme == "Empty")
    {
      gasReactionTerm = new GasReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else if (gasReactionScheme == "Galerkin")
    {
      gasReactionTerm = new GasReactionTerm_2D_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else
    {
      errorInvalidScheme("gasReactionScheme: "+gasReactionScheme);
    }

  }
  else if (dimensions == "AX")
  {
    nDimensions = 2;
    nElementNodes = nDimensions+1;
    nBoundaryElementNodes = nDimensions;

    elementProps = new ElementProps_2D(nDimensions);
    boundaryElementProps = new BoundaryElementProps_2D(nDimensions);

    if (convectionScheme == "Empty")
    {
      convectionTerm = new ConvectionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    //else if (convectionScheme == "Galerkin")
    //{
    //	convectionTerm = new ConvectionTerm_AX_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    //}
    //else if (convectionScheme == "MDC")
    //{
    //	convectionTerm = new ConvectionTerm_AX_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    //}
    else if (convectionScheme == "N")
    {
      convectionTerm = new ConvectionTerm_AX_N(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("convectionScheme: "+convectionScheme);
    }

    if (diffusionScheme == "Empty")
    {
      diffusionTerm = new DiffusionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "Galerkin")
    {
      diffusionTerm = new DiffusionTerm_AX_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "MDC")
    {
      diffusionTerm = new DiffusionTerm_AX_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("diffusionScheme: "+diffusionScheme);
    }

    if (migrationScheme == "Empty")
    {
      migrationTerm = new MigrationTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "Galerkin")
    {
      migrationTerm = new MigrationTerm_AX_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "MDC")
    {
      migrationTerm = new MigrationTerm_AX_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("migrationScheme: "+migrationScheme);
    }

    if (magneticScheme == "Empty")
    {
      magneticTerm = new MagneticTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("magneticScheme: "+magneticScheme);
    }


    if (homReactionScheme == "Empty")
    {
      homReactionTerm = new HomReactionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin")
    {
      homReactionTerm = new HomReactionTerm_AX_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal")
    {
      homReactionTerm = new HomReactionTerm_AX_Galerkin_Diagonalized(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal_M")
    {
      homReactionTerm = new HomReactionTerm_AX_Galerkin_Diagonalized_M(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "MDC")
    {
      homReactionTerm = new HomReactionTerm_AX_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("homReactionScheme: "+homReactionScheme);
    }

    if (electrostaticsScheme == "Empty")
    {
      electrostaticsTerm = new ElectrostaticsTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "Galerkin")
    {
      electrostaticsTerm = new ElectrostaticsTerm_AX_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "MDC")
    {
      electrostaticsTerm = new ElectrostaticsTerm_AX_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("electrostaticsScheme: "+electrostaticsScheme);
    }

    if (timeScheme == "Empty")
    {
      timeTerm = new TimeTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "Galerkin")
    {
      timeTerm = new TimeTerm_AX_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "MDC")
    {
      timeTerm = new TimeTerm_AX_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("timeScheme: "+timeScheme);
    }

    if (elecReactionScheme == "Empty")
    {
      elecReactionTerm = new ElecReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Galerkin")
    {
      elecReactionTerm = new ElecReactionTerm_AX_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "MDC")
    {
      elecReactionTerm = new ElecReactionTerm_AX_MDC(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Pointwise")
    {
      elecReactionTerm = new ElecReactionTerm_AX_Pointwise(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else
    {
      errorInvalidScheme("elecReactionScheme: "+elecReactionScheme);
    }

    if (gasReactionScheme == "Empty")
    {
      gasReactionTerm = new GasReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else if (gasReactionScheme == "Galerkin")
    {
      gasReactionTerm = new GasReactionTerm_AX_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else
    {
      errorInvalidScheme("gasReactionScheme: "+gasReactionScheme);
    }
  }
  else if (dimensions == "3D")
  {
    nDimensions = 3;
    nElementNodes = nDimensions+1;
    nBoundaryElementNodes = nDimensions;

    elementProps = new ElementProps_3D(nDimensions);
    boundaryElementProps = new BoundaryElementProps_3D(nDimensions);

    if (convectionScheme == "Empty")
    {
      convectionTerm = new ConvectionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
//			convectionTerm = new ConvectionTerm_Template<3,Scheme_Null>(nVariables,mitrem,elementProps);
    }
    else if (convectionScheme == "Galerkin")
    {
      convectionTerm = new ConvectionTerm_Template<3,Scheme_Galerkin>(nVariables,mitrem,elementProps);
    }
    //else if (convectionScheme == "MDC")
    //{
    //	convectionTerm = new ConvectionTerm_3D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    //}
    else if (convectionScheme == "N")
    {
      convectionTerm = new ConvectionTerm_3D_N(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
//			convectionTerm = new ConvectionTerm_Template<3,Scheme_N>(nVariables,mitrem,elementProps);
    }
    else if (convectionScheme == "LDA")
    {
      convectionTerm = new ConvectionTerm_Template<3,Scheme_LDA>(nVariables,mitrem,elementProps);
    }
    else
    {
      errorInvalidScheme("convectionScheme: "+convectionScheme);
    }

    if (diffusionScheme == "Empty")
    {
      diffusionTerm = new DiffusionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "Galerkin")
    {
      diffusionTerm = new DiffusionTerm_3D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (diffusionScheme == "MDC")
    {
      diffusionTerm = new DiffusionTerm_3D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("diffusionScheme: "+diffusionScheme);
    }

    if (migrationScheme == "Empty")
    {
      migrationTerm = new MigrationTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "Galerkin")
    {
      migrationTerm = new MigrationTerm_3D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (migrationScheme == "MDC")
    {
      migrationTerm = new MigrationTerm_3D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("migrationScheme: "+migrationScheme);
    }

    if (magneticScheme == "Empty")
    {
      magneticTerm = new MagneticTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (magneticScheme == "Galerkin")
    {
      magneticTerm = new MagneticTerm_3D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("magneticScheme: "+magneticScheme);
    }

    if (homReactionScheme == "Empty")
    {
      homReactionTerm = new HomReactionTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin")
    {
      homReactionTerm = new HomReactionTerm_3D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal")
    {
      homReactionTerm = new HomReactionTerm_3D_Galerkin_Diagonalized(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "Galerkin_Diagonal_M")
    {
      homReactionTerm = new HomReactionTerm_3D_Galerkin_Diagonalized_M(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (homReactionScheme == "MDC")
    {
      homReactionTerm = new HomReactionTerm_3D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("homReactionScheme: "+homReactionScheme);
    }

    if (electrostaticsScheme == "Empty")
    {
      electrostaticsTerm = new ElectrostaticsTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "Galerkin")
    {
      electrostaticsTerm = new ElectrostaticsTerm_3D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (electrostaticsScheme == "MDC")
    {
      electrostaticsTerm = new ElectrostaticsTerm_3D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("electrostaticsScheme: "+electrostaticsScheme);
    }

    if (timeScheme == "Empty")
    {
      timeTerm = new TimeTerm(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "Galerkin")
    {
      timeTerm = new TimeTerm_3D_Galerkin(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else if (timeScheme == "MDC")
    {
      timeTerm = new TimeTerm_3D_MDC(nDimensions, nElementNodes, nVariables, mitrem, elementProps);
    }
    else
    {
      errorInvalidScheme("timeScheme: "+timeScheme);
    }

    if (elecReactionScheme == "Empty")
    {
      elecReactionTerm = new ElecReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Galerkin")
    {
      elecReactionTerm = new ElecReactionTerm_3D_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "MDC")
    {
      elecReactionTerm = new ElecReactionTerm_3D_MDC(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else if (elecReactionScheme == "Pointwise")
    {
      elecReactionTerm = new ElecReactionTerm_3D_Pointwise(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps);
    }
    else
    {
      errorInvalidScheme("elecReactionScheme: "+elecReactionScheme);
    }

    if (gasReactionScheme == "Empty")
    {
      gasReactionTerm = new GasReactionTerm(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else if (gasReactionScheme == "Galerkin")
    {
      gasReactionTerm = new GasReactionTerm_3D_Galerkin(nDimensions, nBoundaryElementNodes, nVariables, mitrem, boundaryElementProps, _is_bubble);
    }
    else
    {
      errorInvalidScheme("gasReactionScheme: "+gasReactionScheme);
    }
  }

  elementContribution = new ElementContribution(nDimensions, nElementNodes, nIons, convectionTerm, diffusionTerm, migrationTerm, magneticTerm, homReactionTerm, electrostaticsTerm, timeTerm);
  boundaryElementContribution = new BoundaryElementContribution(nBoundaryElementNodes, nIons, nElecReactions, nGasReactions, elecReactionTerm, gasReactionTerm);

  // if charge conservation is to be assembled in place of one mass balance
  m_chargeconservation = _charge_conservation;

  // if first and last equations are to be swapped
  m_swap = _swap_first_and_last_equations;
}

//---------------------------------------------------------------------------
ElementMatrixAssembler::~ElementMatrixAssembler()
{
  // delete in reverse order of construction!
  delete boundaryElementContribution;
  delete elementContribution;

  delete gasReactionTerm;
  delete elecReactionTerm;
  delete timeTerm;
  delete electrostaticsTerm;
  delete homReactionTerm;
  delete magneticTerm;
  delete migrationTerm;
  delete diffusionTerm;
  delete convectionTerm;

  delete boundaryElementProps;
  delete elementProps;
}
//---------------------------------------------------------------------------
DoubleMatrix ElementMatrixAssembler::calcElementMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  EmptyDoubleMatrix mat = elementContribution->calcMat(coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
  for (unsigned m=0; m<nElementNodes; m++)
  {
    // make linear combination
    if (m_chargeconservation)
    {
      double zi = (double) mitrem->getIonChargeNumber(0);
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        mat[m*nVariables][n] *= zi;
      }
      for (unsigned i=1; i<nIons; i++)
      {
        zi = (double) mitrem->getIonChargeNumber(i);
        for (unsigned n=0; n<nElementNodes*nVariables; n++)
        {
          mat[m*nVariables][n] += zi*mat[m*nVariables+i][n];
        }
      }
    }
    // swap first and last equation
    if (m_swap)
    {
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        double temp = mat[m*nVariables][n]; // copy first equation to temp
        mat[m*nVariables][n] = mat[m*nVariables+nIons][n]; // copy last equation to first
        mat[m*nVariables+nIons][n] = temp; // copy temp to last equation
      }
    }
  }
  return mat;
}
//---------------------------------------------------------------------------
DoubleMatrix ElementMatrixAssembler::calcElementJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  EmptyDoubleMatrix jac = elementContribution->calcJac(coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
  for (unsigned m=0; m<nElementNodes; m++)
  {
    // make linear combination
    if (m_chargeconservation)
    {
      double zi = (double) mitrem->getIonChargeNumber(0);
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        jac[m*nVariables][n] *= zi;
      }
      for (unsigned i=1; i<nIons; i++)
      {
        zi = (double) mitrem->getIonChargeNumber(i);
        for (unsigned n=0; n<nElementNodes*nVariables; n++)
        {
          jac[m*nVariables][n] += zi*jac[m*nVariables+i][n];
        }
      }
    }
    // swap first and last equation
    if (m_swap)
    {
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        double temp = jac[m*nVariables][n]; // copy first equation to temp
        jac[m*nVariables][n] = jac[m*nVariables+nIons][n]; // copy last equation to first
        jac[m*nVariables+nIons][n] = temp; // copy temp to last equation
      }
    }
  }
  return jac;
}
//---------------------------------------------------------------------------
DoubleVector ElementMatrixAssembler::calcBoundaryElementVec(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
  EmptyDoubleVector vec = boundaryElementContribution->calcVec(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);
  for (unsigned m=0; m<nBoundaryElementNodes; m++)
  {
    // make linear combination
    if (m_chargeconservation)
    {
      double zi = (double) mitrem->getIonChargeNumber(0);
      vec[m*nVariables] *= zi;
      for (unsigned i=1; i<nIons; i++)
      {
        zi = (double) mitrem->getIonChargeNumber(i);
        vec[m*nVariables] += zi*vec[m*nVariables+i];
      }
    }
    // swap first and last equation
    if (m_swap)
    {
      double temp = vec[m*nVariables]; // copy first equation to temp
      vec[m*nVariables] = vec[m*nVariables+nIons]; // copy last equation to first
      vec[m*nVariables+nIons] = temp; // copy temp to last equation
    }
  }
  return vec;
}
//---------------------------------------------------------------------------
DoubleMatrix ElementMatrixAssembler::calcBoundaryElementJac(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
  EmptyDoubleMatrix jac = boundaryElementContribution->calcJac(coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);
  for (unsigned m=0; m<nBoundaryElementNodes; m++)
  {
    // make linear combination
    if (m_chargeconservation)
    {
      double zi = (double) mitrem->getIonChargeNumber(0);
      for (unsigned n=0; n<nBoundaryElementNodes*nVariables; n++)
      {
        jac[m*nVariables][n] *= zi;
      }
      for (unsigned i=1; i<nIons; i++)
      {
        zi = (double) mitrem->getIonChargeNumber(i);
        for (unsigned n=0; n<nBoundaryElementNodes*nVariables; n++)
        {
          jac[m*nVariables][n] += zi*jac[m*nVariables+i][n];
        }
      }
    }
    // swap first and last equation
    if (m_swap)
    {
      for (unsigned n=0; n<nBoundaryElementNodes*nVariables; n++)
      {
        double temp = jac[m*nVariables][n]; // copy first equation to temp
        jac[m*nVariables][n] = jac[m*nVariables+nIons][n]; // copy last equation to first
        jac[m*nVariables+nIons][n] = temp; // copy temp to last equation
      }
    }
  }
  return jac;
}
//---------------------------------------------------------------------------
DoubleMatrix ElementMatrixAssembler::calcElementTimeMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  EmptyDoubleMatrix mat = elementContribution->calcTimeMat(coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
  for (unsigned m=0; m<nElementNodes; m++)
  {
    // make linear combination
    if (m_chargeconservation)
    {
      double zi = (double) mitrem->getIonChargeNumber(0);
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        mat[m*nVariables][n] *= zi;
      }
      for (unsigned i=1; i<nIons; i++)
      {
        zi = (double) mitrem->getIonChargeNumber(i);
        for (unsigned n=0; n<nElementNodes*nVariables; n++)
        {
          mat[m*nVariables][n] += zi*mat[m*nVariables+i][n];
        }
      }
    }
    // swap first and last equation
    if (m_swap)
    {
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        double temp = mat[m*nVariables][n]; // copy first equation to temp
        mat[m*nVariables][n] = mat[m*nVariables+nIons][n]; // copy last equation to first
        mat[m*nVariables+nIons][n] = temp; // copy temp to last equation
      }
    }
  }
  return mat;
}
//---------------------------------------------------------------------------
DoubleMatrix ElementMatrixAssembler::calcElementTimeJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  EmptyDoubleMatrix jac = elementContribution->calcTimeJac(coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
  for (unsigned m=0; m<nElementNodes; m++)
  {
    // make linear combination
    if (m_chargeconservation)
    {
      double zi = (double) mitrem->getIonChargeNumber(0);
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        jac[m*nVariables][n] *= zi;
      }
      for (unsigned i=1; i<nIons; i++)
      {
        zi = (double) mitrem->getIonChargeNumber(i);
        for (unsigned n=0; n<nElementNodes*nVariables; n++)
        {
          jac[m*nVariables][n] += zi*jac[m*nVariables+i][n];
        }
      }
    }
    // swap first and last equation
    if (m_swap)
    {
      for (unsigned n=0; n<nElementNodes*nVariables; n++)
      {
        double temp = jac[m*nVariables][n]; // copy first equation to temp
        jac[m*nVariables][n] = jac[m*nVariables+nIons][n]; // copy last equation to first
        jac[m*nVariables+nIons][n] = temp; // copy temp to last equation
      }
    }
  }
  return jac;
}
//---------------------------------------------------------------------------
void ElementMatrixAssembler::errorInvalidScheme(const std::string &scheme)
{
  std::cout << "ERROR : THE SCHEME " << scheme <<" DOES NOT EXIST." << std::endl;
  std::cin.get();
  exit(1);
}
//---------------------------------------------------------------------------

