//---------------------------------------------------------------------------

#ifndef IonFitH
#define IonFitH

//---------------------------------------------------------------------------

#include <string>
#include <vector>
#include "Electrolyte.h"

//---------------------------------------------------------------------------

class IonFit
{
public :
  IonFit(const std::string &label_, Electrolyte** electrolytes, unsigned nElectrolytes);
  ~IonFit();

  // Set methods
  void    setDiameter(double d);
  void    setMinimumDiameter(double dMin);
  void    setMaximumDiameter(double dMax);
  void    setDiffusionConstant(double D);
  void    setMinimumDiffusionConstant(double DMin);
  void    setMaximumDiffusionConstant(double DMax);
  void    setMolarMass(double M);

  // Get methods
  std::string  getLabel() const;
  double    getDiameter() const;
  double    getMinimumDiameter() const;
  double    getMaximumDiameter() const;
  double    getDiffusionConstant() const;
  double    getMinimumDiffusionConstant() const;
  double    getMaximumDiffusionConstant() const;
  double    getMolarMass() const;
  unsigned  getNElectrolytes() const;
  unsigned  getElectrolytes(unsigned e) const;

  // Methods
  double    calcOsmoticCoefficientDerivativeDiameter(unsigned e,unsigned t);
  double    calcActivityCoefficientDerivativeDiameter(unsigned e,unsigned t);
  double    calcConductivityDerivativeDiameter(unsigned e,unsigned t);
  double    calcTransportNumberDerivativeDiameter(unsigned e,unsigned t);
  double    calcDiffusionCoefficientDerivativeDiameter(unsigned e,unsigned t);
  double    calcOsmoticCoefficientDerivativeDiffusionConstant(unsigned e,unsigned t);
  double    calcActivityCoefficientDerivativeDiffusionConstant(unsigned e,unsigned t);
  double    calcConductivityDerivativeDiffusionConstant(unsigned e,unsigned t);
  double    calcTransportNumberDerivativeDiffusionConstant(unsigned e,unsigned t);
  double    calcDiffusionCoefficientDerivativeDiffusionConstant(unsigned e,unsigned t);

private :
  class IonReference
  {
  public :
    Electrolyte*  electrolyte;
    unsigned    electrolyteIndex;
    unsigned    ionIndex;
  };

  std::string  label;
  double    d,dMin,dMax,D,DMin,DMax,M,epsilon;
  std::vector<IonReference> ionReferences;

  void errorIonNotFound(const std::string label);
};

//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
inline void IonFit::setMinimumDiameter(double dMin_)
{
  this->dMin = dMin_;
}
//---------------------------------------------------------------------------
inline void IonFit::setMaximumDiameter(double dMax_)
{
  this->dMax = dMax_;
}
//---------------------------------------------------------------------------
inline void IonFit::setMinimumDiffusionConstant(double DMin_)
{
  this->DMin = DMin_;
}
//---------------------------------------------------------------------------
inline void IonFit::setMaximumDiffusionConstant(double DMax_)
{
  this->DMax = DMax_;
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
inline std::string IonFit::getLabel() const
{
  return label;
}
//---------------------------------------------------------------------------
inline double IonFit::getDiameter() const
{
  return d;
}
//---------------------------------------------------------------------------
inline double IonFit::getMinimumDiameter() const
{
  return dMin;
}
//---------------------------------------------------------------------------
inline double IonFit::getMaximumDiameter() const
{
  return dMax;
}
//---------------------------------------------------------------------------
inline double IonFit::getDiffusionConstant() const
{
  return D;
}
//---------------------------------------------------------------------------
inline double IonFit::getMinimumDiffusionConstant() const
{
  return DMin;
}
//---------------------------------------------------------------------------
inline double IonFit::getMaximumDiffusionConstant() const
{
  return DMax;
}
//---------------------------------------------------------------------------
inline double IonFit::getMolarMass() const
{
  return M;
}
//---------------------------------------------------------------------------
inline unsigned IonFit::getNElectrolytes() const
{
  return unsigned(ionReferences.size());
}
//---------------------------------------------------------------------------
inline unsigned IonFit::getElectrolytes(unsigned e) const
{
  return ionReferences[e].electrolyteIndex;
}
//---------------------------------------------------------------------------

#endif

