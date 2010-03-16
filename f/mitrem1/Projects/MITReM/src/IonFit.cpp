//---------------------------------------------------------------------------

#include <cstdlib>

#include "IonFit.h"
#include <math.h>
#include <iostream>

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
IonFit::IonFit(const std::string &label_, Electrolyte** electrolytes, unsigned nElectrolytes)
: label(label_)
{
	bool found = false;
	for (unsigned e=0; e<nElectrolytes; e++)
	{
		for (unsigned i=0; i<electrolytes[e]->getNIons(); i++)
		{
			if (label == electrolytes[e]->getIonLabel(i))
			{
				found = true;
				IonReference ionReference;
				ionReference.electrolyte = electrolytes[e];
				ionReference.electrolyteIndex = e;
				ionReference.ionIndex = i;
				ionReferences.push_back(ionReference);
				break;
			}
		}
	}
	epsilon = 1e-3;
}
//---------------------------------------------------------------------------
IonFit::~IonFit()
{
}
//---------------------------------------------------------------------------


//--- METHODS ---------------------------------------------------------------
void IonFit::setDiameter(double d)
{
	this->d = d;
	std::vector<IonReference>::iterator i;
	for (i = ionReferences.begin(); i != ionReferences.end(); i++) 
	{
		i->electrolyte->setIonDiameter(i->ionIndex,d);
	}
}
//---------------------------------------------------------------------------
void IonFit::setDiffusionConstant(double D)
{ 
	this->D = D;
	std::vector<IonReference>::iterator i;
	for (i = ionReferences.begin(); i != ionReferences.end(); i++) 
	{
		i->electrolyte->setIonDiffusionConstant(i->ionIndex,D);
	}
}
//---------------------------------------------------------------------------
void IonFit::setMolarMass(double M)
{ 
	this->M = M;
	std::vector<IonReference>::iterator i;
	for (i = ionReferences.begin(); i != ionReferences.end(); i++) 
	{
		i->electrolyte->setIonMolarMass(i->ionIndex,M);
	}
}
//---------------------------------------------------------------------------
double IonFit::calcOsmoticCoefficientDerivativeDiameter(unsigned e,unsigned t)
{
	double d1 = (1.-epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d1);
	double y1 = ionReferences[e].electrolyte->calcOsmoticCoefficient(t);
	double d2 = (1.+epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d2);
	double y2 = ionReferences[e].electrolyte->calcOsmoticCoefficient(t);

	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d);
	return (y2-y1)/(2.*epsilon*d);
}
//---------------------------------------------------------------------------
double IonFit::calcActivityCoefficientDerivativeDiameter(unsigned e,unsigned t)
{
	double d1 = (1.-epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d1);
	double y1 = ionReferences[e].electrolyte->calcActivityCoefficient(t);
	double d2 = (1.+epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d2);
	double y2 = ionReferences[e].electrolyte->calcActivityCoefficient(t);

	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d);
	return (y2-y1)/(2.*epsilon*d);
}
//---------------------------------------------------------------------------
double IonFit::calcConductivityDerivativeDiameter(unsigned e,unsigned t)
{
	double d1 = (1.-epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d1);
	double y1 = ionReferences[e].electrolyte->calcConductivity(t);
	double d2 = (1.+epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d2);
	double y2 = ionReferences[e].electrolyte->calcConductivity(t);

	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d);
	return (y2-y1)/(2.*epsilon*d);
}
//---------------------------------------------------------------------------
double IonFit::calcTransportNumberDerivativeDiameter(unsigned e,unsigned t)
{
	double d1 = (1.-epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d1);
	double y1 = ionReferences[e].electrolyte->calcTransportNumber(t);
	double d2 = (1.+epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d2);
	double y2 = ionReferences[e].electrolyte->calcTransportNumber(t);

	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d);
	return (y2-y1)/(2.*epsilon*d);
}
//---------------------------------------------------------------------------
double IonFit::calcDiffusionCoefficientDerivativeDiameter(unsigned e,unsigned t)
{
	double d1 = (1.-epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d1);
	double y1 = ionReferences[e].electrolyte->calcDiffusionCoefficient(t);
	double d2 = (1.+epsilon)*d;
	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d2);
	double y2 = ionReferences[e].electrolyte->calcDiffusionCoefficient(t);

	ionReferences[e].electrolyte->setIonDiameter(ionReferences[e].ionIndex,d);
	return (y2-y1)/(2.*epsilon*d);
}
//---------------------------------------------------------------------------
double IonFit::calcOsmoticCoefficientDerivativeDiffusionConstant(unsigned e,unsigned t)
{
	return 0.;
}
//---------------------------------------------------------------------------
double IonFit::calcActivityCoefficientDerivativeDiffusionConstant(unsigned e,unsigned t)
{
	return 0.;
}
//---------------------------------------------------------------------------
double IonFit::calcConductivityDerivativeDiffusionConstant(unsigned e,unsigned t)
{
	double D1 = (1.-epsilon)*D;
	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D1);
	double y1 = ionReferences[e].electrolyte->calcConductivity(t);
	double D2 = (1.+epsilon)*D;
	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D2);
	double y2 = ionReferences[e].electrolyte->calcConductivity(t);

	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D);
	return (y2-y1)/(2.*epsilon*D);
}
//---------------------------------------------------------------------------
double IonFit::calcTransportNumberDerivativeDiffusionConstant(unsigned e,unsigned t)
{
	double D1 = (1.-epsilon)*D;
	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D1);
	double y1 = ionReferences[e].electrolyte->calcTransportNumber(t);
	double D2 = (1.+epsilon)*D;
	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D2);
	double y2 = ionReferences[e].electrolyte->calcTransportNumber(t);

	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D);
	return (y2-y1)/(2.*epsilon*D);
}
//---------------------------------------------------------------------------
double IonFit::calcDiffusionCoefficientDerivativeDiffusionConstant(unsigned e,unsigned t)
{
	double D1 = (1.-epsilon)*D;
	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D1);
	double y1 = ionReferences[e].electrolyte->calcDiffusionCoefficient(t);
	double D2 = (1.+epsilon)*D;
	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D2);
	double y2 = ionReferences[e].electrolyte->calcDiffusionCoefficient(t);

	ionReferences[e].electrolyte->setIonDiffusionConstant(ionReferences[e].ionIndex,D);
	return (y2-y1)/(2.*epsilon*D);
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void IonFit::errorIonNotFound(const std::string label)
{
	std::cout << "ERROR IN IonFit.cpp.\nTHE ION " << label 
		<< " WAS NOWHERE FOUND IN THE *.ions FILES." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------

