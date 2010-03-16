//---------------------------------------------------------------------------

#include "Electrolyte.h"
#include "DatFileReader.h"
#include "MathematicsPhysics.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Electrolyte::Electrolyte(const std::string &name_)
: MITReM(name_)
{
	name = name_;

	osmoticCoefficientsFile = name + ".osmoticcoeff";
	activityCoefficientsFile = name + ".activitycoeff";
	conductivitiesFile = name + ".conductivities";
	transportNumbersFile = name + ".transportnumbers";
	diffusionCoefficientsFile = name + ".diffusioncoeff";

	readOsmoticCoefficients();
	readActivityCoefficients();
	readConductivities();
	readTransportNumbers();
	readDiffusionCoefficients();
}
//---------------------------------------------------------------------------
Electrolyte::~Electrolyte()
{
}
//---------------------------------------------------------------------------


//--- READ METHODS ----------------------------------------------------------
void Electrolyte::readOsmoticCoefficients()
{
	DatFileReader datFile(osmoticCoefficientsFile);	

	nOsmoticCoefficients = datFile.readMultipleVector_nVectors("[nOsmoticCoefficients]");
	osmoticCoefficients = new Experiment[nOsmoticCoefficients];
	for (unsigned t=0; t<nOsmoticCoefficients; t++) {
		datFile.readMultipleVector("[nOsmoticCoefficients]", t, "<c>", osmoticCoefficients[t].c);
		datFile.readMultipleVector("[nOsmoticCoefficients]", t, "<y>", osmoticCoefficients[t].y);
	}
}
//---------------------------------------------------------------------------
void Electrolyte::readActivityCoefficients()
{
	DatFileReader datFile(activityCoefficientsFile);	

	nActivityCoefficients = datFile.readMultipleVector_nVectors("[nActivityCoefficients]");
	activityCoefficients = new Experiment[nActivityCoefficients];
	for (unsigned t=0; t<nActivityCoefficients; t++) {
		datFile.readMultipleVector("[nActivityCoefficients]", t, "<c>", activityCoefficients[t].c);
		datFile.readMultipleVector("[nActivityCoefficients]", t, "<y>", activityCoefficients[t].y);
	}
}
//---------------------------------------------------------------------------
void Electrolyte::readConductivities()
{
	DatFileReader datFile(conductivitiesFile);	

	nConductivities = datFile.readMultipleVector_nVectors("[nConductivities]");
	conductivities = new Experiment[nConductivities];
	for (unsigned t=0; t<nConductivities; t++) {
		datFile.readMultipleVector("[nConductivities]", t, "<c>", conductivities[t].c);
		datFile.readMultipleVector("[nConductivities]", t, "<y>", conductivities[t].y);
	}
}
//---------------------------------------------------------------------------
void Electrolyte::readTransportNumbers()
{
	DatFileReader datFile(transportNumbersFile);	

	nTransportNumbers = datFile.readMultipleVector_nVectors("[nTransportNumbers]");
	transportNumbers = new Experiment[nTransportNumbers];
	for (unsigned t=0; t<nTransportNumbers; t++) {
		datFile.readMultipleVector("[nTransportNumbers]", t, "<c>", transportNumbers[t].c);
		datFile.readMultipleVector("[nTransportNumbers]", t, "<y>", transportNumbers[t].y);
	}
}
//---------------------------------------------------------------------------
void Electrolyte::readDiffusionCoefficients()
{
	DatFileReader datFile(diffusionCoefficientsFile);	

	nDiffusionCoefficients = datFile.readMultipleVector_nVectors("[nDiffusionCoefficients]");
	diffusionCoefficients = new Experiment[nDiffusionCoefficients];
	for (unsigned t=0; t<nDiffusionCoefficients; t++) {
		datFile.readMultipleVector("[nDiffusionCoefficients]", t, "<c>", diffusionCoefficients[t].c);
		datFile.readMultipleVector("[nDiffusionCoefficients]", t, "<y>", diffusionCoefficients[t].y);
	}
}
//---------------------------------------------------------------------------


//--- METHODS ---------------------------------------------------------------
void Electrolyte::calcEquilibrium()
{
	MITReM::calcEquilibrium();
	int* stoichCation = new int[nHomReactions];
	int* stoichAnion = new int[nHomReactions];
	for (unsigned r=0; r<nHomReactions; r++)
	{
		stoichCation[r] = homReactions[r]->getStoichReag(0);
		stoichAnion[r] = homReactions[r]->getStoichReag(1);
	}
	electrolyteModel->calcBinaryOnsagerCoefficients_SF(stoichCation,stoichAnion);
	delete[] stoichCation;
	delete[] stoichAnion;
}
//---------------------------------------------------------------------------
void Electrolyte::initOsmoticCoefficient(unsigned t)
{
	if (t < nOsmoticCoefficients)
	{
		for (unsigned i=0; i<2; i++) 
		{
			electrolyteSolution->setIonConcentration(i,s[i]*osmoticCoefficients[t].c);
		}
		for (unsigned i=2; i<electrolyteSolution->getNIons(); i++) 
		{
			electrolyteSolution->setIonConcentration(i,0);
		}
		calcEquilibrium();
	}
	else errorExperimentDoesNotExist(t);
}
//---------------------------------------------------------------------------
void Electrolyte::initActivityCoefficient(unsigned t)
{
	if (t < nActivityCoefficients)
	{
		for (unsigned i=0; i<2; i++) 
		{
			electrolyteSolution->setIonConcentration(i,s[i]*activityCoefficients[t].c);
		}
		for (unsigned i=2; i<electrolyteSolution->getNIons(); i++) 
		{
			electrolyteSolution->setIonConcentration(i,0);
		}
		calcEquilibrium();
	}
	else errorExperimentDoesNotExist(t);
}
//---------------------------------------------------------------------------
void Electrolyte::initConductivity(unsigned t)
{
	if (t < nConductivities)
	{
		for (unsigned i=0; i<2; i++) 
		{
			electrolyteSolution->setIonConcentration(i,s[i]*conductivities[t].c);
		}
		for (unsigned i=2; i<electrolyteSolution->getNIons(); i++) 
		{
			electrolyteSolution->setIonConcentration(i,0);
		}
		calcEquilibrium();
	}
	else errorExperimentDoesNotExist(t);
}
//---------------------------------------------------------------------------
void Electrolyte::initTransportNumber(unsigned t)
{
	if (t < nTransportNumbers)
	{
		for (unsigned i=0; i<2; i++) 
		{
			electrolyteSolution->setIonConcentration(i,s[i]*transportNumbers[t].c);
		}
		for (unsigned i=2; i<electrolyteSolution->getNIons(); i++) 
		{
			electrolyteSolution->setIonConcentration(i,0);
		}
		calcEquilibrium();
	}
	else errorExperimentDoesNotExist(t);
}
//---------------------------------------------------------------------------
void Electrolyte::initDiffusionCoefficient(unsigned t)
{
	if (t < nDiffusionCoefficients)
	{
		for (unsigned i=0; i<2; i++) 
		{
			electrolyteSolution->setIonConcentration(i,s[i]*diffusionCoefficients[t].c);
		}
		for (unsigned i=2; i<electrolyteSolution->getNIons(); i++) 
		{
			electrolyteSolution->setIonConcentration(i,0);
		}
		calcEquilibrium();
	}
	else errorExperimentDoesNotExist(t);
}
//---------------------------------------------------------------------------
double Electrolyte::calcOsmoticCoefficient(unsigned t)
{
	initOsmoticCoefficient(t);
	double totalConcentration = (s[0]+s[1])*osmoticCoefficients[t].c;
	return electrolyteModel->calcOsmoticCoefficient_MM(totalConcentration);
}
//---------------------------------------------------------------------------
double Electrolyte::calcActivityCoefficient(unsigned t)
{
	initActivityCoefficient(t);
	return electrolyteModel->calcMeanActivityCoefficient_MM(s[0],s[1],activityCoefficients[t].c);
}
//---------------------------------------------------------------------------
double Electrolyte::calcConductivity(unsigned t)
{
	initConductivity(t);
	double normality = electrolyteSolution->getIonChargeNumber(0)*s[0]*conductivities[t].c;
	return electrolyteModel->calcEquivalentConductivity(normality);
}
//---------------------------------------------------------------------------
double Electrolyte::calcTransportNumber(unsigned t)
{
	initTransportNumber(t);
	return electrolyteModel->calcCationTransportNumber_SF();
}
//---------------------------------------------------------------------------
double Electrolyte::calcDiffusionCoefficient(unsigned t)
{
	initDiffusionCoefficient(t);
	return electrolyteModel->calcMolarThermodynamicDiffusionCoefficient_SF(s[0],s[1],diffusionCoefficients[t].c);
}
//---------------------------------------------------------------------------
double Electrolyte::calcL00N(unsigned t)
{
	initConductivity(t);
	double normality = electrolyteSolution->getIonChargeNumber(0)*s[0]*conductivities[t].c;
	return electrolyteModel->calcEquivalentBinaryOnsagerCoefficient_CC_SF(normality);
}
//---------------------------------------------------------------------------
double Electrolyte::calcL01N(unsigned t)
{
	initConductivity(t);
	double normality = electrolyteSolution->getIonChargeNumber(0)*s[0]*conductivities[t].c;
	return electrolyteModel->calcEquivalentBinaryOnsagerCoefficient_CA_SF(normality);
}
//---------------------------------------------------------------------------
double Electrolyte::calcL11N(unsigned t)
{
	initConductivity(t);
	double normality = electrolyteSolution->getIonChargeNumber(0)*s[0]*conductivities[t].c;
	return electrolyteModel->calcEquivalentBinaryOnsagerCoefficient_AA_SF(normality);
}
//---------------------------------------------------------------------------
double Electrolyte::calcLijN(unsigned t, unsigned i, unsigned j)
{
	initConductivity(t);
	double normality = electrolyteSolution->getIonChargeNumber(0)*s[0]*conductivities[t].c;
	return electrolyteModel->calcOnsagerCoefficient_SF(i,j)/normality;
}
//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
void Electrolyte::setHomReactionEquilibriumConstant(unsigned r, double K)
{ 
	homReactions[r]->setEquilibriumConstant(K);
	double kb = homReactions[r]->getBackwardRateConstant();
	double kf = K*kb;
	homReactions[r]->setForwardRateConstant(kf);
}
//---------------------------------------------------------------------------
void Electrolyte::setMaximumConcentration(double cMax)
{
	this->cMax = cMax;

	unsigned nOsmoticCoefficientsTemp = 0;
	for (unsigned t=0; t<nOsmoticCoefficients; t++) {
		if (osmoticCoefficients[t].c <= cMax) nOsmoticCoefficientsTemp++;
		else break;
	}
	nOsmoticCoefficients = nOsmoticCoefficientsTemp;

	unsigned nActivityCoefficientsTemp = 0;
	for (unsigned t=0; t<nActivityCoefficients; t++) {
		if (activityCoefficients[t].c <= cMax) nActivityCoefficientsTemp++;
		else break;
	}
	nActivityCoefficients = nActivityCoefficientsTemp;

	unsigned nConductivitiesTemp = 0;
	for (unsigned t=0; t<nConductivities; t++) {
		if (conductivities[t].c <= cMax) nConductivitiesTemp++;
		else break;
	}
	nConductivities = nConductivitiesTemp;

	unsigned nTransportNumbersTemp = 0;
	for (unsigned t=0; t<nTransportNumbers; t++) {
		if (transportNumbers[t].c <= cMax) nTransportNumbersTemp++;
		else break;
	}
	nTransportNumbers = nTransportNumbersTemp;

	unsigned nDiffusionCoefficientsTemp = 0;
	for (unsigned t=0; t<nDiffusionCoefficients; t++) {
		if (diffusionCoefficients[t].c <= cMax) nDiffusionCoefficientsTemp++;
		else break;
	}
	nDiffusionCoefficients = nDiffusionCoefficientsTemp;
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
double Electrolyte::getL00N(unsigned t) const
{
	double F2 = F_CONST*F_CONST;
	double RT = R_CONST*getSolutionTemperature();
	int z[2];
	z[0] = getIonChargeNumber(0);
	z[1] = getIonChargeNumber(1);
	double t0 = transportNumbers[t].y;
	return (t0*t0*conductivities[t].y/F2 - diffusionCoefficients[t].y/RT*z[0]*z[1]/(z[0]-z[1]))/(z[0]*z[0]);
}
//---------------------------------------------------------------------------
double Electrolyte::getL01N(unsigned t) const
{
	double F2 = F_CONST*F_CONST;
	double RT = R_CONST*getSolutionTemperature();
	int z[2];
	z[0] = getIonChargeNumber(0);
	z[1] = getIonChargeNumber(1);
	double t0 = transportNumbers[t].y;
	double t1 = 1.-t0;
	return (t0*t1*conductivities[t].y/F2 + diffusionCoefficients[t].y/RT*z[0]*z[1]/(z[0]-z[1]))/(z[0]*z[1]);
}
//---------------------------------------------------------------------------
double Electrolyte::getL11N(unsigned t) const
{
	double F2 = F_CONST*F_CONST;
	double RT = R_CONST*getSolutionTemperature();
	int z[2];
	z[0] = getIonChargeNumber(0);
	z[1] = getIonChargeNumber(1);
	double t1 = 1.-transportNumbers[t].y;
	return (t1*t1*conductivities[t].y/F2 - diffusionCoefficients[t].y/RT*z[0]*z[1]/(z[0]-z[1]))/(z[1]*z[1]);
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void Electrolyte::errorExperimentDoesNotExist(unsigned t)
{
	std::cout << "ERROR IN Electrolyte.cpp.\nTHE INDEX " << t 
		<< " EXCEEDS THE NUMBER OF EXPERIMENTS THAT WERE PERFORMED." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------

