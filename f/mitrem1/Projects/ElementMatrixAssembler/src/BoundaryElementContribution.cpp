//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

#include "BoundaryElementContribution.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
BoundaryElementContribution::BoundaryElementContribution(
	unsigned nBoundaryElementNodes_, 
	unsigned nIons_, 
	unsigned nElecReactionsMax_, 
	unsigned nGasReactionsMax_, 
	ElecReactionTerm* elecReactionTerm_,
	GasReactionTerm* gasReactionTerm_) : 
		elecReactionTerm(elecReactionTerm_), 
		gasReactionTerm(gasReactionTerm_),
		nBoundaryElementNodes(nBoundaryElementNodes_), 
		nIons(nIons_),
		nElecReactionsMax(nElecReactionsMax_), 
		nGasReactionsMax(nGasReactionsMax_)
{	
	size = nBoundaryElementNodes*(nIons+1);
	boundaryElementVec = new double[size];
	boundaryElementJac= new double*[size];
	boundaryElementElecReactionCurrentDensity = new double*[nBoundaryElementNodes];
	boundaryElementGasReactionRate = new double*[nBoundaryElementNodes];
	for (unsigned m=0; m<size; m++)
	{
		boundaryElementJac[m] = new double[size];
	}
	for (unsigned m=0; m<nBoundaryElementNodes; m++)
	{
		boundaryElementElecReactionCurrentDensity[m] = new double[nElecReactionsMax];
		boundaryElementGasReactionRate[m] = new double[nGasReactionsMax];	
	}
}	
//---------------------------------------------------------------------------
BoundaryElementContribution::~BoundaryElementContribution()
{
	for (unsigned m=0; m<nBoundaryElementNodes; m++)
	{
		delete[] boundaryElementGasReactionRate[m];	
		delete[] boundaryElementElecReactionCurrentDensity[m];
	}
	for (unsigned m=0; m<size; m++)
	{
		delete[] boundaryElementJac[m];
	}
	delete[] boundaryElementElecReactionCurrentDensity;	
	delete[] boundaryElementGasReactionRate;
	delete[] boundaryElementJac;
	delete[] boundaryElementVec;
}	
//---------------------------------------------------------------------------
EmptyDoubleVector BoundaryElementContribution::calcVec(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	for (unsigned m=0; m<size; m++)
	{
		boundaryElementVec[m] = 0.;
	}
	
	elecReactionTerm->calcVec(boundaryElementVec, coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);
	gasReactionTerm->calcVec(boundaryElementVec, coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);

	return boundaryElementVec;
}
//---------------------------------------------------------------------------
EmptyDoubleMatrix BoundaryElementContribution::calcJac(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential, unsigned nGasReactions, IndexList gasReactions)
{
	for (unsigned m=0; m<size; m++)
	{
		for (unsigned n=0; n<size; n++)
		{
			boundaryElementJac[m][n] = 0.;
		}
	}
	
	elecReactionTerm->calcJac(boundaryElementJac, coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);
	gasReactionTerm->calcJac(boundaryElementJac, coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential, nGasReactions, gasReactions);

	return boundaryElementJac;
}
//---------------------------------------------------------------------------
DoubleListList BoundaryElementContribution::calcElecReactionCurrentDensities(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nElecReactions, IndexList elecReactions, double electrodePotential)
{
	for (unsigned m=0; m<nBoundaryElementNodes; m++)
	{
		for (unsigned r=0; r<nElecReactionsMax; r++)
		{
			boundaryElementElecReactionCurrentDensity[m][r] = 0.;
		}
	}

	elecReactionTerm->calcElecReactionCurrentDensities(boundaryElementElecReactionCurrentDensity, coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nElecReactions, elecReactions, electrodePotential);

	return boundaryElementElecReactionCurrentDensity;
}
//---------------------------------------------------------------------------
DoubleListList BoundaryElementContribution::calcGasReactionRates(DoubleVectorList coordinates, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList surfaceGasFractions, unsigned nGasReactions, IndexList gasReactions)
{
	for (unsigned m=0; m<nBoundaryElementNodes; m++)
	{
		for (unsigned r=0; r<nGasReactionsMax; r++)
		{
			boundaryElementGasReactionRate[m][r] = 0.;
		}
	}

	gasReactionTerm->calcGasReactionRates(boundaryElementGasReactionRate, coordinates, concentrations, potentials, temperatures, densities, surfaceGasFractions, nGasReactions, gasReactions);

	return boundaryElementGasReactionRate;
}
//---------------------------------------------------------------------------

