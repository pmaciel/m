//---------------------------------------------------------------------------

#ifndef ElectrolyteModel_MSAH
#define ElectrolyteModel_MSAH

//---------------------------------------------------------------------------

#include <vector>
#include <string>
#include "ElectrolyteModel_DH.h"
#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------

class ElectrolyteModel_MSA : public ElectrolyteModel_DH
{
public :
	ElectrolyteModel_MSA(ElectrolyteSolution* electrolyteSolution_);
	virtual ~ElectrolyteModel_MSA();

	// Methods
	void			init(bool verbose);
	//double			calcConductivityBest() const;

	// Binary electrolytes
	double		calcOsmoticCoefficient_MM(double totalConcentration) const;

protected :
	// Members
	double		Omega, Psi, Gamma;				// Coulomb parameters
	double		X0, X1, X2, X3, Delta, exclVol;	// hard sphere parameters
	double*		FGamma;
	double*		FPsi;
	double*		alpha;
	double*		beta;
	double		theta,phi;
	double		Delta_2, Delta_3, Delta_4, 
				X2DivX3, X2DivX3_2, X2DivX3_3,
				X2LnDeltaDivX3X3, FLong1, FLong2;
	double		dav;					// "electrostatic" average diameter
	double*		d;

	// Methods
	void			calcHardSphereParameters();
	double			calcXp(unsigned exp) const;
	void			calcCoulombParameters();
	void			calcOmega(double g);
	void			calcPsi(double g);
	double			calcGammaFunction(double g) const;
	double			calcGammaFunctionDerivative(double g) const;
	double			calcElectrophoreticCorrection(unsigned i, unsigned j) const;
	//double			calcElectrophoreticCorrection(unsigned i) const;
	double			calcRelaxationCorrection(unsigned i, unsigned j) const;	
	//double			calcRelaxationCorrection(unsigned i) const;
};

//---------------------------------------------------------------------------

#endif

