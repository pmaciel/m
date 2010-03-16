//---------------------------------------------------------------------------

#define _USE_MATH_DEFINES

//---------------------------------------------------------------------------

#include "ElectrolyteModel_MSA.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "MathematicsPhysics.h"

//--------------------------------------------------------------------------- 


//---------------------------------------------------------------------------
ElectrolyteModel_MSA::ElectrolyteModel_MSA(ElectrolyteSolution* electrolyteSolution_)
: ElectrolyteModel_DH(electrolyteSolution_)
{	
	conductivityCorrectionFactor = 1.;
	FPsi = new double[nIons];
	FGamma = new double[nIons];
	alpha = new double[nIons];
	beta = new double[nIons];
	d = new double[nIons];
}
//---------------------------------------------------------------------------
ElectrolyteModel_MSA::~ElectrolyteModel_MSA()
{
	delete[] FPsi;
	delete[] FGamma;
	delete[] alpha;
	delete[] beta;
	delete[] d;
}
//---------------------------------------------------------------------------
void ElectrolyteModel_MSA::init(bool verbose)
{
	this->verbose = verbose;
	T = electrolyteSolution->getSolutionTemperature();
	viscosity = electrolyteSolution->getSolventDynamicViscosity();
	dielectricConstant = electrolyteSolution->getSolventDielectricConstant();
	density = electrolyteSolution->getSolutionDensity();
	for (unsigned i=0; i<nIons; i++)
	{
		z[i] = electrolyteSolution->getIonChargeNumber(i);
		c[i] = electrolyteSolution->getIonConcentration(i);
		D[i] = electrolyteSolution->getIonDiffusionConstant(i);
		M[i] = electrolyteSolution->getIonMolarMass(i);
		d[i] = electrolyteSolution->getIonDiameter(i);
	}
	
	LB = e_CONST*e_CONST/(4.*M_PI*kB_CONST*T*dielectricConstant);

	// calculate parameters
	calcConcentrationParameters();
	calcHardSphereParameters();
	calcCoulombParameters();
	calcEigenproblem();

	std::ofstream output;
	/*if (verbose) 
	{
		output.setf(std::ios::scientific);
		output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
		double totalConcentration = c[0]+c[1]   +c[2];
		//if (nIons == 3) totalConcentration += 2*c[2];
		output << totalConcentration << '\t';
		for (unsigned i=0; i<nIons; i++)
			output << c[i] << '\t';
		output << Gamma << '\t';
		double osmoticPressureC = (-Gamma*Gamma*Gamma/3. - 2.*LB*Psi*Psi)/M_PI;
		double osmoticPressureIdHS = (X0/Delta + 3.*X1*X2/Delta_2 + (3.-X3)*X2*X2*X2/Delta_3)*6./M_PI;
		output << (osmoticPressureC + osmoticPressureIdHS)/(totalConcentration*NA_CONST) << '\t';
		output.close();
	}*/

	double F0 = -log(Delta);
	double F1 = 3.*X2/Delta;
	double F2 = 3.*X1/Delta
				+ 3.*X2*X2/(X3*Delta*Delta)
				+ 3.*X2*X2*log(Delta)/(X3*X3);
	double F3 = (X0 - X2*X2*X2/(X3*X3))/Delta
				+ 3.*X1*X2/(Delta*Delta)
				- 2.*X2*X2*X2*log(Delta)/(X3*X3*X3)
				+ X2*X2*X2/(X3*Delta*Delta)*(2./Delta - 1./X3);

	for (unsigned i=0; i<nIons; i++) 
	{	
		// Activity coefficients in the McMillan-Mayer reference frame
		// THEY SHOULD BE CONVERTED TO THE LEWIS-RANDALL REFERENCE FRAME!!!
		double lnfE = -LB*(Gamma*z[i]*z[i]/FGamma[i] + Psi*d[i]*
					  ((2.*z[i] - Psi*d[i]*d[i])/FGamma[i] + Psi*d[i]*d[i]/3.));
		double lnfHS = F0 + F1*d[i] + F2*d[i]*d[i] + F3*d[i]*d[i]*d[i];
		f[i] = exp(lnfE + lnfHS);

		/*if (verbose) 
		{
			output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
			output << f[i] << '\t';
			output.close();
		}*/

		double dX0 = M_PI/6.;
		double dX1 = dX0*d[i];
		double dX2 = dX1*d[i];
		double dX3 = dX2*d[i];

		double G0 = dX3/Delta;
		double G1 = dX2*3./Delta
					+ dX3*3.*X2/Delta_2;
		double G2 = dX1*3./Delta
					+ dX2*(6.*X2DivX3/Delta_2 + 6.*X2LnDeltaDivX3X3)
					+ dX3*FLong1;
		double G3 = dX0/Delta
					+ dX1*3.*X2/Delta_2
					+ dX2*FLong1
					+ dX3*FLong2;
		
		for (unsigned j=i; j<nIons; j++) 
		{
			// Partial derivatives of the activity coefficients --> d(lnfi)/d(cj) = d(lnfj)/d(ci)
			// THEY SHOULD BE CONVERTED TO THE LEWIS-RANDALL REFERENCE FRAME!!!
			double dlnfC = -M_PI*LB*(alpha[i]*alpha[j]/(Delta*Omega) + LB/2.*beta[i]*beta[j]/
							(Gamma + M_PI*LB*(phi - M_PI*theta*theta/(2.*Delta*Omega))));
			double dlnfHS = G0 + G1*d[j] + G2*d[j]*d[j] + G3*d[j]*d[j]*d[j];
			dlnf[i][j] = (dlnfC + dlnfHS)*NA_CONST;
			dlnf[j][i] = dlnf[i][j];
			
			// Lij in the solvent-fixed reference frame --> Lij=Lji
			double EC = calcElectrophoreticCorrection(i,j);
			double RC = calcRelaxationCorrection(i,j);
			ls[i][j] = (kroneck(i,j)*c[i]*D[i]/(R_CONST*T) + EC + RC)*conductivityCorrectionFactor;
			ls[j][i] = ls[i][j];

			/*if (verbose) 
			{
				output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
				output << ls[i][j] << '\t';
				output.close();
			}*/
		}
	}	
	
	double contrib3 = 0;
	for (unsigned k=0; k<nIons; k++) 
	{	
		for (unsigned l=0; l<nIons; l++) 
		{	
			contrib3 += M[k]*M[l]*ls[k][l];
		}
	}
	for (unsigned i=0; i<nIons; i++) 
	{	
		for (unsigned j=i; j<nIons; j++) 
		{			
			// Lij in the mass-fixed reference frame
			lm[i][j] = ls[i][j];
			double contrib1 = 0;
			double contrib2 = 0;
			for (unsigned k=0; k<nIons; k++) 
			{
				contrib1 += M[k]*ls[k][j];
				contrib2 += M[k]*ls[i][k];
			}
			lm[i][j] = ls[i][j] 
								 - c[i]/density*contrib1 
								 - c[j]/density*contrib2
								 + c[i]*c[j]/(density*density)*contrib3;
			lm[j][i] = lm[i][j];
		}
	}

	/*if (verbose) 
	{
		output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
		output << std::endl;
		output.close();
	}
	exit(1);*/
}
//---------------------------------------------------------------------------
// Hard Sphere parameters (X0, X1, X2 and X3 and also Delta)
void ElectrolyteModel_MSA::calcHardSphereParameters()
{
	X0 = calcXp(0);
	X1 = calcXp(1);
	X2 = calcXp(2);
	X3 = calcXp(3);

	Delta = 1. - X3;

	exclVol = (3.*X1*X2/X0+X3)/4.;

	X2DivX3 = X2/X3;
	X2DivX3_2 = X2DivX3*X2DivX3;
	X2DivX3_3 = X2DivX3_2*X2DivX3;
	Delta_2 = Delta*Delta;
	Delta_3 = Delta_2*Delta;
	Delta_4 = Delta_3*Delta;
	X2LnDeltaDivX3X3 = X2DivX3*log(Delta)/X3;
	FLong1 = -3.*X2DivX3_2/Delta
			 + 3.*(X1-X2DivX3_2)/Delta_2
			 + 6.*X2*X2DivX3/Delta_3
			 - 6.*X2DivX3*X2LnDeltaDivX3X3;
	FLong2 = 5.*X2DivX3_3/Delta
			 + (X0+X2DivX3_3)/Delta_2
			 + 2.*X2*(3.*X1-2.*X2DivX3_2)/Delta_3
			 + 6.*X2*X2*X2DivX3/Delta_4
			 + 6.*X2DivX3_2*X2LnDeltaDivX3X3;
}
//---------------------------------------------------------------------------
// Xp
double ElectrolyteModel_MSA::calcXp(unsigned exp) const
{
	double Xp = 0.;
	for (unsigned i=0; i<nIons; i++) 
	{
		double factor = 1.;
		for (unsigned p=0; p<exp; p++)
		{
			factor *= d[i];
		}
		Xp += c[i]*factor;
	}
	Xp *= NA_CONST*M_PI/6.;
	return Xp;
}
//---------------------------------------------------------------------------
// Coulomb parameters (Omega, Psi and Gamma; Delta is calculated in calcHardSphereParameters)
void ElectrolyteModel_MSA::calcCoulombParameters()
{
	double gammaOld;

	// Starting value is the half inverse Debye length
	Gamma = HIDL;
	
	// Newton iterations
	unsigned iter = 0;
	const unsigned iterMax = 10;
	do {
		gammaOld = Gamma;
		calcOmega(gammaOld);
		calcPsi(gammaOld);
		Gamma = gammaOld-calcGammaFunction(gammaOld)/calcGammaFunctionDerivative(gammaOld);
		iter++;
	} while ((fabs(1.-Gamma/gammaOld) > 1e-9) && (iter < iterMax)); 
	// Stop when relative change is small enough. 

	theta = 0.;
	phi = 0.;
	for (unsigned i=0; i<nIons; i++) 
	{
		FGamma[i] = 1. + Gamma*d[i];
		FPsi[i] = (z[i] - Psi*d[i]*d[i])/FGamma[i];
		theta += c[i]*d[i]*d[i]*FPsi[i]/FGamma[i];
		phi += c[i]*d[i]*FPsi[i]*FPsi[i]/FGamma[i];
	}
	theta *= NA_CONST;
	phi *= NA_CONST;

	for (unsigned i=0; i<nIons; i++) 
	{
		alpha[i] = d[i]*(FPsi[i] + Psi*d[i]*d[i]/3.);
		beta[i] = FPsi[i]*FPsi[i] - M_PI*alpha[i]*theta/(Delta*Omega);
	}

	//dav = 0;
	//for (unsigned i=0; i<nIons; i++) 
	//{
	//	dav += z[i]*z[i]*c[i]*d[i];
	//}
	//dav /= I;
	dav = (HIDL/Gamma - 1.)/Gamma;
}
//---------------------------------------------------------------------------
// Omega
void ElectrolyteModel_MSA::calcOmega(double g)
{
	Omega = 0.;
	for (unsigned i=0; i<nIons; i++) 
	{
		Omega += c[i]*d[i]*d[i]*d[i]/(1.+g*d[i]);
	}
	Omega *= NA_CONST*M_PI/(2.*Delta);
	Omega += 1.;
}
//---------------------------------------------------------------------------
// Psi
void ElectrolyteModel_MSA::calcPsi(double g)
{
	Psi = 0.;
	for (unsigned i=0; i<nIons; i++) 
	{
		Psi += c[i]*d[i]*z[i]/(1.+g*d[i]);
	}
	Psi *= NA_CONST*M_PI/(2.*Delta*Omega);
}
//---------------------------------------------------------------------------
// Gamma function
double ElectrolyteModel_MSA::calcGammaFunction(double g) const
{
	double sum = 0.;
	for (unsigned i=0; i<nIons; i++) 
	{
		double factor = (z[i]-Psi*d[i]*d[i])/(1.+g*d[i]);
		sum += c[i]*factor*factor;
	}
	return g*g - M_PI*LB*NA_CONST*sum;
}
//---------------------------------------------------------------------------
// Derivative Gamma function (actually an approximation with d(Psi)/d(Gamma)=0 )
double ElectrolyteModel_MSA::calcGammaFunctionDerivative(double g) const
{
	double sum = 0.;
	for (unsigned i=0; i<nIons; i++) 
	{
		double factor = (z[i]-Psi*d[i]*d[i])/(1.+g*d[i]);
		sum += c[i]*d[i]*factor*factor/(1.+g*d[i]);
	}
	double dgFunct = 2.*(g+M_PI*LB*NA_CONST*sum);
	//if (dgFunct <= 0) errorZero("ElectrolyteModel.cpp","calcGammaFunctionDerivative");
	return dgFunct;
}
//---------------------------------------------------------------------------
// Electrophoretic correction for the Onsager coefficients
double ElectrolyteModel_MSA::calcElectrophoreticCorrection(unsigned i, unsigned j) const
{
	double ECorrC1, ECorrC2, ECorrHS;
	double dij = 0.5*(d[i]+d[j]);
	
	// Coulomb (1)
	double sumA = 0.;
	for (unsigned k=0; k<nIons; k++) 
	{
		sumA += c[k]*z[k]*z[k]*d[k]/(FGamma[k]*FGamma[k]);
	}
	sumA *= M_PI*LB*NA_CONST;
	if (sumA <= 0) errorZero("ElectrolyteModel.cpp","calcElectrophoreticCorrection");
	ECorrC1 = -LB*z[i]*z[j]/(FGamma[i]*FGamma[j]*(Gamma+sumA));

	// Coulomb (2) IS ALWAYS > 0
	// ----------------------------------------------------------------------------	//
	// This is a second order term and SHOULD be smaller than the first order term.	//
	// Unfortunately, this may not always be so for |z[i]|>1.						//
	// ----------------------------------------------------------------------------	//
	double x = 4.*HIDL*dij; 
	ECorrC2 = LB*LB*z[i]*z[i]*z[j]*z[j]*exp(x)*expInt(x)/(FGamma[i]*FGamma[i]*FGamma[j]*FGamma[j]);

	// Hard Sphere IS ALWAYS < 0
	if (X0 <= 0) errorZero("ElectrolyteModel.cpp","calcElectrophoreticCorrection");
	ECorrHS = -dij*dij*(1.-exclVol/5.+exclVol*exclVol/10.)/(1.+2.*exclVol);

	/*if (verbose) 
	{
		std::ofstream output;
		output.setf(std::ios::scientific);
		output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
		output << ECorrC1*c[j]/(3.*viscosity) << '\t' << ECorrC2*c[j]/(3.*viscosity) << '\t' << ECorrHS*c[j]/(3.*viscosity) << '\t';
		output.close();
	}*/

	// Total
	return (ECorrC1 + ECorrC2 + ECorrHS)*c[i]*c[j]/(3.*viscosity);
}
//---------------------------------------------------------------------------
// Electrophoretic correction for the conductivity (Psi!=0)
//double ElectrolyteModel_MSA::calcElectrophoreticCorrection(unsigned i) const
//{
//	// THIS FORMULA IS NEVER USED
//	
//	double ECorrC1,/*ECorrC2,*//*ECorrHS;
//	
//	// Coulomb (1)
//	ECorrC1 = -(z[i]*Gamma + Psi*d[i])/(M_PI*NA_CONST*FGamma[i]);
//
//	// Coulomb (2) and Hard Sphere
//	//ECorrC2 = 0.;
//	ECorrHS = 0.;
//	for (unsigned j=0; j<nIons; j++) {
//		double dij = (d[i]+d[j])/2.;
//		//double x = 4.*HIDL*dij;
//		//ECorrC2 += z[j]*z[j]*z[j]*ions[j].c/(FGamma[j]*FGamma[j])*exp(x)*expInt(x);
//		ECorrHS -= z[j]*c[j]*dij*dij;
//	}
//	//ECorrC2 *= LB*LB*z[i]*z[i]/(FGamma[i]*FGamma[i]);
//	ECorrHS *= (1.-exclVol/5.+exclVol*exclVol/10.)/(1.+2.*exclVol);
//
//	return (ECorrC1 + /*ECorrC2 + *//*ECorrHS)/(3.*viscosity);
//}
//---------------------------------------------------------------------------
// Relaxation correction for the Onsager coefficients
double ElectrolyteModel_MSA::calcRelaxationCorrection(unsigned i, unsigned j) const
{
	double RCorr = 0.;
	double q,sqrtq,expfac; 
	// no need to include p=0, because it produces a zero term
	for (unsigned p=1; p<nIons; p++) 
	{
		q = Eigenvalues[p];
		sqrtq = sqrt(Eigenvalues[p]);
		expfac = exp(-2.*HIDL*sqrtq*dav)/((1.+Gamma*dav)*(1.+Gamma*dav));
		RCorr += sinh(2.*HIDL*sqrtq*dav)/dav
				 *(1.-q)*expfac*Eigenvectors[p][i]*Eigenvectors[p][j]
				 /(q + Gamma/HIDL*sqrtq + Gamma*Gamma/(2.*HIDL*HIDL) - expfac/2.);
	}

	/*if (verbose) 
	{
		std::ofstream output;
		output.setf(std::ios::scientific);
		output.open("Sensitivity.xls", std::ofstream::out | std::ofstream::app);
		output << -LB/(3*cond0*R_CONST*T)*RCorr/c[i] << '\t';
		output.close();
	}*/

	return -LB/(3*cond0*R_CONST*T)*RCorr;
}
//---------------------------------------------------------------------------
// Relaxation correction for the conductivity
/*double ElectrolyteModel_MSA::calcRelaxationCorrection(unsigned i) const
{
	// THIS FORMULA IS NEVER USED
	
	double RCorr = 0.;
	double q,sqrtq,expfac,SumEigenvectors; 
	// no need to include p=0, because it produces a zero term
	for (unsigned p=1; p<nIons; p++) 
	{
		q = Eigenvalues[p];
		sqrtq = sqrt(Eigenvalues[p]);
		expfac = exp(-2.*HIDL*sqrtq*dav)/((1.+Gamma*dav)*(1.+Gamma*dav));
		SumEigenvectors = 0.;
		for (unsigned j=0; j<nIons; j++)
		{
			SumEigenvectors += z[j]*Eigenvectors[p][j];
		}
		RCorr += sinh(2.*HIDL*sqrtq*dav)/dav
				 *(1.-q)*expfac*Eigenvectors[p][i]*SumEigenvectors
				 /(q + Gamma/HIDL*sqrtq + Gamma*Gamma/(2.*HIDL*HIDL) - expfac/2.);
	}
	return -LB/(3*cond0*R_CONST*T)*RCorr;
}*/
//---------------------------------------------------------------------------
double ElectrolyteModel_MSA::calcOsmoticCoefficient_MM(double totalConcentration) const 
{
	// Coulomb
	double osmoticPressureC = (-Gamma*Gamma*Gamma/3. - 2.*LB*Psi*Psi)/M_PI;
	
	// Ideal + Hard Sphere
	double osmoticPressureIdHS = (X0/Delta + 3.*X1*X2/Delta_2 + (3.-X3)*X2*X2*X2/Delta_3)*6./M_PI;
	
	return (osmoticPressureC + osmoticPressureIdHS)/(totalConcentration*NA_CONST);
}
//---------------------------------------------------------------------------
// Conductivity
/*double ElectrolyteModel_MSA::calcConductivityBest() const
{
	// THIS FORMULA IS NEVER USED
	
	double kappa = cond0/(R_CONST*T);
	for (unsigned i=0; i<nIons; i++) {
		if (z[i] != 0) {
			kappa += z[i]*(c[i]*calcElectrophoreticCorrection(i) + calcRelaxationCorrection(i));
		}
	}
	return kappa*F_CONST*F_CONST;
}*/
//---------------------------------------------------------------------------

