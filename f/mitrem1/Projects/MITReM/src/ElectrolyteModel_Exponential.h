//---------------------------------------------------------------------------

#ifndef ElectrolyteModel_ExponentialH
#define ElectrolyteModel_ExponentialH

//---------------------------------------------------------------------------

#include <vector>
#include <string>
#include "ElectrolyteModel.h"
#include "ElectrolyteSolution.h"

//---------------------------------------------------------------------------

class ElectrolyteModel_Exponential : public ElectrolyteModel
{
public :
  ElectrolyteModel_Exponential(ElectrolyteSolution* electrolyteSolution_);
  virtual ~ElectrolyteModel_Exponential();

  // Methods
  virtual void  init(bool verbose);
  virtual double  calcConductivityCorrectionFactor(double experimentalConductivity);
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
  unsigned  nIons;
  int*    z;
  double*    c;
  double*    D;
  double*    DLim;
  double*    DPow;
  double    cTotal;
};

//---------------------------------------------------------------------------

#endif

