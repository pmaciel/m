//---------------------------------------------------------------------------

#ifndef HomReactionH
#define HomReactionH

//---------------------------------------------------------------------------

#include "ElectrolyteSolution.h"
#include "ElectrolyteModel.h"

//---------------------------------------------------------------------------

class HomReaction
{
public :
  HomReaction(ElectrolyteSolution* electrolyteSolution_, ElectrolyteModel* electrolyteModel_, unsigned nReagents_, unsigned nProducts_);
  ~HomReaction();

  // Set methods
  void    setLabel(const std::string &label);
  void    setForwardRateConstant(double kf);
  void    setBackwardRateConstant(double kb);
  void    setReagents(unsigned i, unsigned agent);
  void    setStoichReag(unsigned i, int stoich);
  void    setProducts(unsigned i, unsigned agent);
  void    setStoichProd(unsigned i, int stoich);
  void    setEquilibriumConstant(double K);

  // Get methods
  unsigned    getNReagents() const;
  unsigned    getNProducts() const;
  std::string  getLabel() const;
  double      getForwardRateConstant() const;
  double      getBackwardRateConstant() const;
  unsigned    getReagents(unsigned i) const;
  int          getStoichReag(unsigned i) const;
  unsigned    getProducts(unsigned i) const;
  int          getStoichProd(unsigned i) const;
  double      getEquilibriumConstant() const;
  int          getStoichOf(unsigned i) const;

  // Methods
  double    calcForwardRateConstant() const;
  double    calcBackwardRateConstant() const;
  double    calcReactionRate() const;
  double    progressToEquilibrium() const;
  double    calcDeviationFromEquilibrium() const;
  double    calcRelativeDeviationFromEquilibrium() const;
  double    calcDeviationFromEquilibriumDerivative(int* stoichList) const;

protected :
  // Members
  std::string  label;
  double      kf,kb,K;
  unsigned    nReagents,nProducts;
  unsigned*    reagents;
  unsigned*    products;
  int*        stoichReag;
  int*        stoichProd;
  ElectrolyteSolution*  electrolyteSolution;
  ElectrolyteModel*      electrolyteModel;
  double*      cReagentsSave;
  double*      cProductsSave;

  // Methods
  double    calcMaximumBackwardProgress() const;
  double    calcMaximumForwardProgress() const;
  double    calcEquilibriumProgress(double xMin, double xMax) const;
  void      progress(double x) const;
  double    calcActivityProduct() const;
};

//---------------------------------------------------------------------------


//--- SET METHODS -----------------------------------------------------------
inline void HomReaction::setLabel(const std::string &label_)
{
  this->label = label_;
}
//---------------------------------------------------------------------------
inline void HomReaction::setForwardRateConstant(double kf_)
{
  this->kf = kf_;
}
//---------------------------------------------------------------------------
inline void HomReaction::setBackwardRateConstant(double kb_)
{
  this->kb = kb_;
}
//---------------------------------------------------------------------------
inline void HomReaction::setReagents(unsigned i, unsigned agent)
{
  reagents[i] = agent;
}
//---------------------------------------------------------------------------
inline void HomReaction::setStoichReag(unsigned i, int stoich)
{
  stoichReag[i] = stoich;
}
//---------------------------------------------------------------------------
inline void HomReaction::setProducts(unsigned i, unsigned agent)
{
  products[i] = agent;
}
//---------------------------------------------------------------------------
inline void HomReaction::setStoichProd(unsigned i, int stoich)
{
  stoichProd[i] = stoich;
}
//---------------------------------------------------------------------------
inline void HomReaction::setEquilibriumConstant(double K_)
{
  this->K = K_;
}
//---------------------------------------------------------------------------


//--- GET METHODS -----------------------------------------------------------
inline unsigned HomReaction::getNReagents() const
{
  return nReagents;
}
//---------------------------------------------------------------------------
inline unsigned HomReaction::getNProducts() const
{
  return nProducts;
}
//---------------------------------------------------------------------------
inline std::string HomReaction::getLabel() const
{
  return label;
}
//---------------------------------------------------------------------------
inline double HomReaction::getForwardRateConstant() const
{
  return kf;
}
//---------------------------------------------------------------------------
inline double HomReaction::getBackwardRateConstant() const
{
  return kb;
}
//---------------------------------------------------------------------------
inline unsigned HomReaction::getReagents(unsigned i) const
{
  return reagents[i];
}
//---------------------------------------------------------------------------
inline int HomReaction::getStoichReag(unsigned i) const
{
  return stoichReag[i];
}
//---------------------------------------------------------------------------
inline unsigned HomReaction::getProducts(unsigned i) const
{
  return products[i];
}
//---------------------------------------------------------------------------
inline int HomReaction::getStoichProd(unsigned i) const
{
  return stoichProd[i];
}
//---------------------------------------------------------------------------
inline double HomReaction::getEquilibriumConstant() const
{
  return K;
}
//---------------------------------------------------------------------------

#endif

