//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------

#include "ElectrolyteModel_Exponential.h"
#include <iostream>
#include <math.h>
#include "MathematicsPhysics.h"

//--------------------------------------------------------------------------- 


//---------------------------------------------------------------------------
ElectrolyteModel_Exponential::ElectrolyteModel_Exponential(ElectrolyteSolution* electrolyteSolution_) 
	: ElectrolyteModel(electrolyteSolution_) 
{
	conductivityCorrectionFactor = 1.;
	nIons = electrolyteSolution->getNIons();
	z = new int[nIons];
	c = new double[nIons];
	D = new double[nIons];
	DLim = new double[nIons];
	DPow = new double[nIons];
}
//---------------------------------------------------------------------------
ElectrolyteModel_Exponential::~ElectrolyteModel_Exponential() 
{
	delete[] z;
	delete[] c;
	delete[] D;
	delete[] DLim;
	delete[] DPow;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_Exponential::init(bool verbose)
{
	cTotal = 0.;
	for (unsigned i=0; i<nIons-1; i++) // don't add the water concentration, which must have index nIons!
	{
		cTotal += c[i]; 
	}
	T = electrolyteSolution->getSolutionTemperature();
	for (unsigned i=0; i<nIons; i++)
	{
		z[i] = electrolyteSolution->getIonChargeNumber(i);
		c[i] = electrolyteSolution->getIonConcentration(i);
		D[i] = electrolyteSolution->getIonDiffusionConstant(i);
		DLim[i] = electrolyteSolution->getIonDiffusionConstantLimit(i);
		DPow[i] = electrolyteSolution->getIonDiffusionConstantPower(i);

		D[i] = (D[i]-DLim[i])*exp(-DPow[i]*cTotal)+DLim[i];
	}
}
//---------------------------------------------------------------------------
// Dij
double ElectrolyteModel_Exponential::calcDiffusionFactor(unsigned i, unsigned j) const
{
	return kroneck(i,j)*D[i]*conductivityCorrectionFactor;
}
//---------------------------------------------------------------------------
// Wi = zi*e*ui
double ElectrolyteModel_Exponential::calcMigrationFactor(unsigned i) const
{
	return F_CONST*z[i]*D[i]*conductivityCorrectionFactor/(R_CONST*T);
}
//---------------------------------------------------------------------------
// Conductivity
double ElectrolyteModel_Exponential::calcConductivity() const
{
	double kappa = 0.;
	for (unsigned i=0; i<nIons; i++) {
		kappa += c[i]*D[i]*z[i]*z[i];
	}
	return kappa*conductivityCorrectionFactor*F_CONST*F_CONST/(R_CONST*T);
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Exponential::calcOnsagerCoefficient_SF(unsigned i, unsigned j) const
{
	return c[i]*NA_CONST*kroneck(i,j)*D[i]*conductivityCorrectionFactor;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Exponential::calcActivityCoefficient_MM(unsigned i) const
{
	return 1.;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Exponential::calcActivityCoefficientDerivative_MM(unsigned i, unsigned j) const
{
	return 0.;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Exponential::calcOsmoticCoefficient_MM(double totalConcentration) const
{
	return 1.;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Exponential::calcMeanActivityCoefficient_MM(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const
{
	return 1.;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_Exponential::calcBinaryOnsagerCoefficients_SF(int* stoichCation, int* stoichAnion)
{
	Ls00 = D[0]*conductivityCorrectionFactor*c[0];
	Ls01 = 0.;
	Ls11 = D[1]*conductivityCorrectionFactor*c[1];
	unsigned nAssociatedIons = nIons-2;
	for (unsigned i=0; i<nAssociatedIons; i++)
	{
		unsigned i2 = i+2;
		Ls00 += D[i2]*conductivityCorrectionFactor*c[i2]*stoichCation[i]*stoichCation[i];
		Ls11 += D[i2]*conductivityCorrectionFactor*c[i2]*stoichAnion[i]*stoichAnion[i];
	}
	Ls00 /= R_CONST*T;
	Ls11 /= R_CONST*T;
}
//---------------------------------------------------------------------------
double ElectrolyteModel_Exponential::calcConductivityCorrectionFactor(double experimentalConductivity)
{
	conductivityCorrectionFactor = experimentalConductivity/calcConductivity();
	return conductivityCorrectionFactor;
}
//---------------------------------------------------------------------------

