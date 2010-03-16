//---------------------------------------------------------------------------

#include "TimeTerm_1D_HCU.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
TimeTerm_1D_HCU::TimeTerm_1D_HCU(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: TimeTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	alpha = new double[nNodes];
}
//---------------------------------------------------------------------------
TimeTerm_1D_HCU::~TimeTerm_1D_HCU()
{
	delete[] alpha;
}
//---------------------------------------------------------------------------
void TimeTerm_1D_HCU::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions)
{
	unsigned inlet = 0;
	const double averageVelocity = 0.5*(velocities[0][0] + velocities[1][0]);

	// Confusing in 1D: side 1 is the side lying opposite of node 1, which is node 2!!!
	for (unsigned m=0; m<nNodes; m++) 
	{
		DoubleVector normal = elementProps->calcNormal(m,coordinates);
		const double k = normal[0]*averageVelocity;
		alpha[m] = 0.;
		if (k > 0.) 
		{
			inlet = m;
			alpha[m] = 1.;
		}
	}
	
	elementSize = elementProps->calcSize(coordinates);
	double contribution1 = elementSize/6.;
	
	for (unsigned i=0; i<nIons; i++) 
	{
		double pecleti = averageVelocity*elementSize/mitrem->getIonDiffusionConstant(i);
		double zetai = pecleti/(pecleti+1.);
		for (unsigned m=0; m<nNodes; m++) 
		{
			double contribution2 = 0.5*elementSize*zetai*(alpha[m] - 0.5);
			unsigned eqmi = eq(m,i);
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,i)] += contribution1 + contribution2;
			}
			elementMat[eqmi][var(m,i)] += contribution1; // diagonal
		}
	}
}
//---------------------------------------------------------------------------
void TimeTerm_1D_HCU::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions)
{
}
//---------------------------------------------------------------------------

