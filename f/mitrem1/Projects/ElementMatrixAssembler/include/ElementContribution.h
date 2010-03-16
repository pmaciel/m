//---------------------------------------------------------------------------

#ifndef ElementContributionH
#define ElementContributionH

//---------------------------------------------------------------------------

#include "TypeDefs.h"
#include "ConvectionTerm.h"
#include "DiffusionTerm.h"
#include "MigrationTerm.h"
#include "MagneticTerm.h"
#include "HomReactionTerm.h"
#include "ElectrostaticsTerm.h"
#include "TimeTerm.h"

//---------------------------------------------------------------------------

class ElementContribution
{
public :
	ElementContribution(
		unsigned nDimensions_, 
		unsigned nElementNodes_, 
		unsigned nIons_, 
		ConvectionTerm* convectionTerm_, 
		DiffusionTerm* diffusionTerm_, 
		MigrationTerm* migrationTerm_, 
		MagneticTerm* magneticTerm_, 
		HomReactionTerm* homReactionTerm_, 
		ElectrostaticsTerm* electrostaticsTerm_, 
		TimeTerm* timeTerm_);
	~ElementContribution();

	EmptyDoubleMatrix	calcMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
	EmptyDoubleMatrix	calcJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

	EmptyDoubleMatrix	calcTimeMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);
	EmptyDoubleMatrix	calcTimeJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

	DoubleVectorList	calcIonCurrentDensities(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors);

protected:
	ConvectionTerm*		convectionTerm;
	DiffusionTerm*		diffusionTerm;
	MigrationTerm*		migrationTerm;
	MagneticTerm*			magneticTerm;
	HomReactionTerm*	homReactionTerm;
	ElectrostaticsTerm*	electrostaticsTerm;
	TimeTerm*			timeTerm;
	
	unsigned			nDimensions,nElementNodes,nIons,size;
	EmptyDoubleMatrix	elementMat;
	EmptyDoubleMatrix	elementTimeMat;
	EmptyDoubleMatrix	elementJac;
	EmptyDoubleMatrix	elementTimeJac;
	EmptyEmptyDoubleVectorList	elementCurr;
};

//---------------------------------------------------------------------------

#endif

