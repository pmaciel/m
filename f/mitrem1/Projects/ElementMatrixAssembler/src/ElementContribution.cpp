//---------------------------------------------------------------------------

#include "ElementContribution.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElementContribution::ElementContribution(
	unsigned nDimensions_, 
	unsigned nElementNodes_, 
	unsigned nIons_, 
	ConvectionTerm* convectionTerm_, 
	DiffusionTerm* diffusionTerm_, 
	MigrationTerm* migrationTerm_, 
	MagneticTerm* magneticTerm_, 
	HomReactionTerm* homReactionTerm_, 
	ElectrostaticsTerm* electrostaticsTerm_, 
	TimeTerm* timeTerm_) :
		convectionTerm(convectionTerm_),
		diffusionTerm(diffusionTerm_),
		migrationTerm(migrationTerm_),
		magneticTerm(magneticTerm_),
		homReactionTerm(homReactionTerm_),
		electrostaticsTerm(electrostaticsTerm_),
		timeTerm(timeTerm_),
		nDimensions(nDimensions_),
		nElementNodes(nElementNodes_),
		nIons(nIons_)
{
	size = nElementNodes*(nIons+1);
	elementMat = new double*[size];
	elementTimeMat = new double*[size];
	elementJac= new double*[size];
	elementTimeJac= new double*[size];
	elementCurr = new double*[nIons+1];
	for (unsigned m=0; m<size; m++)
	{
		elementMat[m] = new double[size];
		elementTimeMat[m] = new double[size];
		elementJac[m] = new double[size];
		elementTimeJac[m] = new double[size];
	}
	for (unsigned i=0; i<nIons+1; i++)
	{
		elementCurr[i] = new double[nDimensions];
	}
}	
//---------------------------------------------------------------------------
ElementContribution::~ElementContribution()
{
	for (unsigned m=0; m<size; m++)
	{
		delete[] elementMat[m];
		delete[] elementTimeMat[m];
		delete[] elementJac[m];
		delete[] elementTimeJac[m];
	}
	for (unsigned i=0; i<nIons+1; i++)
	{
		delete[] elementCurr[i];
	}
	delete[] elementMat;
	delete[] elementTimeMat;
	delete[] elementJac;
	delete[] elementTimeJac;
	delete[] elementCurr;
}	
//---------------------------------------------------------------------------
EmptyDoubleMatrix ElementContribution::calcMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	for (unsigned m=0; m<size; m++)
	{
		for (unsigned n=0; n<size; n++)
		{
			elementMat[m][n] = 0.;
		}
	}

	convectionTerm->calcMat(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	diffusionTerm->calcMat(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	migrationTerm->calcMat(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	magneticTerm->calcMat(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	homReactionTerm->calcMat(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	electrostaticsTerm->calcMat(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);

	return elementMat;
}
//---------------------------------------------------------------------------
EmptyDoubleMatrix ElementContribution::calcJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	for (unsigned m=0; m<size; m++)
	{
		for (unsigned n=0; n<size; n++)
		{
			elementJac[m][n] = 0.;
		}
	}
	
	convectionTerm->calcJac(elementJac, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	diffusionTerm->calcJac(elementJac, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	migrationTerm->calcJac(elementJac, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	magneticTerm->calcJac(elementMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	homReactionTerm->calcJac(elementJac, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);
	electrostaticsTerm->calcJac(elementJac, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);

	return elementJac;
}
//---------------------------------------------------------------------------
EmptyDoubleMatrix ElementContribution::calcTimeMat(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	for (unsigned m=0; m<size; m++)
	{
		for (unsigned n=0; n<size; n++)
		{
			elementTimeMat[m][n] = 0.;
		}
	}
	
	timeTerm->calcMat(elementTimeMat, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);

	return elementTimeMat;
}
//---------------------------------------------------------------------------
EmptyDoubleMatrix ElementContribution::calcTimeJac(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	for (unsigned m=0; m<size; m++)
	{
		for (unsigned n=0; n<size; n++)
		{
			elementTimeJac[m][n] = 0.;
		}
	}
	
	timeTerm->calcJac(elementTimeJac, coordinates, velocities, concentrations, potentials, temperatures, densities, volumeGasFractions, magneticFieldVectors);

	return elementTimeJac;
}
//---------------------------------------------------------------------------
DoubleVectorList ElementContribution::calcIonCurrentDensities(DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
	for (unsigned i=0; i<nIons+1; i++)
	{	
		for (unsigned c=0; c<nDimensions; c++)
		{			
			elementCurr[i][c] = 0.;
		}
	}

	convectionTerm->calcIonCurrentDensities(elementCurr,coordinates,velocities,concentrations,potentials,temperatures,densities,volumeGasFractions,magneticFieldVectors);
	diffusionTerm->calcIonCurrentDensities(elementCurr,coordinates,velocities,concentrations,potentials,temperatures,densities,volumeGasFractions,magneticFieldVectors);
	migrationTerm->calcIonCurrentDensities(elementCurr,coordinates,velocities,concentrations,potentials,temperatures,densities,volumeGasFractions,magneticFieldVectors);
	magneticTerm->calcIonCurrentDensities(elementCurr,coordinates,velocities,concentrations,potentials,temperatures,densities,volumeGasFractions,magneticFieldVectors);
	
	for (unsigned i=0; i<nIons; i++)
	{
		for (unsigned c=0; c<nDimensions; c++)
		{
			elementCurr[nIons][c] += elementCurr[i][c];
		}
	}

	return elementCurr;
}
//---------------------------------------------------------------------------

