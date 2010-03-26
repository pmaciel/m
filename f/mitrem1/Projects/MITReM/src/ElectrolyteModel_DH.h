//---------------------------------------------------------------------------

#ifndef ElectrolyteModel_DHH
#define ElectrolyteModel_DHH

//---------------------------------------------------------------------------

#include <vector>
#include <string>
#include "ElectrolyteModel_Ideal.h"
#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

class ElectrolyteModel_DH : public ElectrolyteModel_Ideal
{
public :
  ElectrolyteModel_DH(ElectrolyteSolution* electrolyteSolution_);
  virtual ~ElectrolyteModel_DH();

  // Methods
  virtual void  init(bool verbose);
  virtual double  calcDiffusionFactor(unsigned i, unsigned j) const;
  virtual double  calcMigrationFactor(unsigned i) const;
  virtual double  calcConductivity() const;
  virtual double  calcOnsagerCoefficient_SF(unsigned i, unsigned j) const;
  virtual double  calcActivityCoefficient_MM(unsigned i) const;
  virtual double  calcActivityCoefficientDerivative_MM(unsigned i, unsigned j) const;

  // Binary electrolytes
  virtual double  calcOsmoticCoefficient_MM(double totalConcentration) const;
  virtual double  calcMeanActivityCoefficient_MM(unsigned sCation, unsigned sAnion, double electrolyteConcentration) const;
  virtual void  calcBinaryOnsagerCoefficients_SF(int* stoichCation, int* stoichAnion);

protected :
  // Members
  double    LB;            // Bjerrum length
  double    I;            // 2x ionic strength
  double    HIDL;          // half inverse Debye length
  double    cond0;          // sum zi^2*ci*Di, is proportional to the conductivity at infinite dilution
  double*    Alphas;
  double*    NormFactors;
  double*    Eigenvalues;
  double**  Eigenvectors;
  double**  ls;            // array of lij in reference frame of solvent
  double**  lm;            // array of lij in reference frame of mass average
  double*    f;            // array of fi
  double**  dlnf;          // array of d(ln(fi))/d(rhoj)
  double    viscosity;
  double    dielectricConstant;
  double    density;
  double*    M;
  double    Mscs;          // M(solvent)*c(solvent)

  // Methods
  void      calcConcentrationParameters();
  virtual double  calcElectrophoreticCorrection(unsigned i, unsigned j) const;
  //virtual double  calcElectrophoreticCorrection(unsigned i) const;
  virtual double  calcRelaxationCorrection(unsigned i, unsigned j) const;
  //virtual double  calcRelaxationCorrection(unsigned i) const;
  void      calcEigenproblem();
  double      calcAlphaFunction(double a) const;
  double      calcAlphaFunctionDerivative(double a) const;

  // Error messages
  void      errorZero(const std::string &cppfile, const std::string &function) const;
};

//---------------------------------------------------------------------------

#endif

