//---------------------------------------------------------------------------

#ifndef ElectrolyteModelH
#define ElectrolyteModelH

//---------------------------------------------------------------------------

#include <vector>
#include <string>
#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------

class ElectrolyteModel
{
public :
	ElectrolyteModel(ElectrolyteSolution* electrolyteSolution_);
	virtual ~ElectrolyteModel();

	// Methods
	virtual void	init(bool verbose) = 0;
	virtual double	calcConductivityCorrectionFactor(double experimentalConductivity) = 0;	// Correction factor for electrolyteModel properties to satisfy experimental conductivity
	virtual double	calcDiffusionFactor(unsigned i, unsigned j) const = 0;
	virtual double	calcMigrationFactor(unsigned i) const = 0;
	virtual double	calcConductivity() const = 0;
	virtual double	calcOnsagerCoefficient_SF(unsigned i, unsigned j) const = 0;	// Solvent-fixed Onsager coefficient
	virtual double	calcActivityCoefficient_MM(unsigned i) const = 0;	// McMillan-Mayer activity coefficient
	virtual double	calcActivityCoefficientDerivative_MM(unsigned i, unsigned j) const = 0;
	
	// Binary electrolytes
	virtual double	calcOsmoticCoefficient_MM(double totalConcentration) const = 0; //McMillan-Mayer osmotic coefficient
	virtual double	calcMeanActivityCoefficient_MM(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const = 0; // McMillan-Mayer mean activity coefficient
	virtual void	calcBinaryOnsagerCoefficients_SF(int* stoichCation, int* stoichAnion) = 0;
	double			calcEquivalentBinaryOnsagerCoefficient_CC_SF(double normality) const;
	double			calcEquivalentBinaryOnsagerCoefficient_CA_SF(double normality) const;
	double			calcEquivalentBinaryOnsagerCoefficient_AA_SF(double normality) const;
	double			calcEquivalentConductivity(double normality) const;
	double			calcMolarThermodynamicDiffusionCoefficient_SF(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const; // Solvent-fixed thermodynamic diffusion coefficient divided by concentration
	double			calcCationTransportNumber_SF() const; // solvent fixed cation electrolyteModel number

protected :
	// Members
	ElectrolyteSolution*	electrolyteSolution;
	double					conductivityCorrectionFactor;
	double					Ls00,Ls11,Ls01;			// binary lij in reference frame of solvent
	double					T;

	bool verbose;
};

//---------------------------------------------------------------------------

#endif

