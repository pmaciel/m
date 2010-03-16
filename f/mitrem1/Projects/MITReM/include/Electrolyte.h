//---------------------------------------------------------------------------

#ifndef ElectrolyteH
#define ElectrolyteH

//---------------------------------------------------------------------------

#include <string>
#include "MITReM.h"

//---------------------------------------------------------------------------

class Electrolyte : public MITReM 
{
public :
	Electrolyte(const std::string &name_);
	~Electrolyte();

	void		calcEquilibrium();
	void		initOsmoticCoefficient(unsigned t);
	void		initActivityCoefficient(unsigned t);
	void		initConductivity(unsigned t);
	void		initTransportNumber(unsigned t);
	void		initDiffusionCoefficient(unsigned t);
	double		calcOsmoticCoefficient(unsigned t);
	double		calcActivityCoefficient(unsigned t);
	double		calcConductivity(unsigned t);
	double		calcTransportNumber(unsigned t);
	double		calcDiffusionCoefficient(unsigned t);
	double		calcL00N(unsigned t);
	double		calcL01N(unsigned t);
	double		calcL11N(unsigned t);
	double		calcLijN(unsigned t, unsigned i, unsigned j);

	// Set methods
	void		setStoichCoefficients(unsigned sCation, unsigned sAnion);
	void		setMaximumConcentration(double cMax);
	void		setHomReactionEquilibriumConstant(unsigned r, double K);

	// Get methods
	std::string	getName() const;
	unsigned	getNOsmoticCoefficients() const;
	unsigned	getNActivityCoefficients() const;
	unsigned	getNConductivities() const;
	unsigned	getNTransportNumbers() const;
	unsigned	getNDiffusionCoefficients() const;
	double		getOsmoticCoefficient(unsigned t) const;
	double		getActivityCoefficient(unsigned t) const;
	double		getConductivity(unsigned t) const;
	double		getTransportNumber(unsigned t) const;
	double		getDiffusionCoefficient(unsigned t) const;
	double		getOsmoticCoefficientConcentration(unsigned t) const;
	double		getActivityCoefficientConcentration(unsigned t) const;
	double		getConductivityConcentration(unsigned t) const;
	double		getTransportNumberConcentration(unsigned t) const;
	double		getDiffusionCoefficientConcentration(unsigned t) const;
	double		getHomReactionEquilibriumConstant(unsigned r) const;
	double		getL00N(unsigned t) const;
	double		getL01N(unsigned t) const;
	double		getL11N(unsigned t) const;

protected :
	class Experiment
	{
	public :
		double			c;	// concentration
		double			y;	// experimental value
	};

	std::string			name;
	std::string			osmoticCoefficientsFile;
	std::string			activityCoefficientsFile;
	std::string			conductivitiesFile;
	std::string			transportNumbersFile;
	std::string			diffusionCoefficientsFile;
	unsigned			nOsmoticCoefficients;
	Experiment*			osmoticCoefficients;
	unsigned			nActivityCoefficients;
	Experiment*			activityCoefficients;
	unsigned			nConductivities;
	Experiment*			conductivities;
	unsigned			nTransportNumbers;
	Experiment*			transportNumbers;
	unsigned			nDiffusionCoefficients;						
	Experiment*			diffusionCoefficients;
	unsigned			s[2];
	double				cMax;

	void readOsmoticCoefficients();
	void readActivityCoefficients();
	void readConductivities();
	void readTransportNumbers();
	void readDiffusionCoefficients();
	void errorExperimentDoesNotExist(unsigned t);
};

//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
inline void Electrolyte::setStoichCoefficients(unsigned sCation, unsigned sAnion)
{
	s[0] = sCation;
	s[1] = sAnion;
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
inline std::string Electrolyte::getName() const
{
	return name;
}
//---------------------------------------------------------------------------
inline unsigned Electrolyte::getNOsmoticCoefficients() const
{
	return nOsmoticCoefficients;
}
//---------------------------------------------------------------------------
inline unsigned Electrolyte::getNActivityCoefficients() const
{
	return nActivityCoefficients;
}
//---------------------------------------------------------------------------
inline unsigned Electrolyte::getNConductivities() const
{
	return nConductivities;
}
//---------------------------------------------------------------------------
inline unsigned Electrolyte::getNTransportNumbers() const
{
	return nTransportNumbers;
}
//---------------------------------------------------------------------------
inline unsigned Electrolyte::getNDiffusionCoefficients() const
{
	return nDiffusionCoefficients;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getOsmoticCoefficient(unsigned t) const
{
	return osmoticCoefficients[t].y;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getActivityCoefficient(unsigned t) const
{
	return activityCoefficients[t].y;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getConductivity(unsigned t) const
{
	return conductivities[t].y;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getTransportNumber(unsigned t) const
{
	return transportNumbers[t].y;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getDiffusionCoefficient(unsigned t) const
{
	return diffusionCoefficients[t].y;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getOsmoticCoefficientConcentration(unsigned t) const
{
	return osmoticCoefficients[t].c;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getActivityCoefficientConcentration(unsigned t) const
{
	return activityCoefficients[t].c;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getConductivityConcentration(unsigned t) const
{
	return conductivities[t].c;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getTransportNumberConcentration(unsigned t) const
{
	return transportNumbers[t].c;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getDiffusionCoefficientConcentration(unsigned t) const
{
	return diffusionCoefficients[t].c;
}
//---------------------------------------------------------------------------
inline double Electrolyte::getHomReactionEquilibriumConstant(unsigned r) const
{
	return homReactions[r]->getEquilibriumConstant();
}
//---------------------------------------------------------------------------

#endif

