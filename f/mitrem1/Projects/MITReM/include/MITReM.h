//---------------------------------------------------------------------------

#ifndef MITReMH
#define MITReMH

//---------------------------------------------------------------------------

#include <string>
#include <vector>

#include "xmlParser.h"
#include "ElectrolyteSolution.h"
#include "ElectrolyteModel.h"
#include "Electrostatics.h"
#include "HomReaction.h"
#include "ElecReaction.h"
#include "GasReaction.h"

//---------------------------------------------------------------------------

class MITReM
{
public :
  MITReM(const std::string &name);
  MITReM(const std::string& s_database, const std::string& s_ec_label);
  ~MITReM();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Methods
  void			init(const double* c, double U, double T, double density);
  void			calcEquilibrium(std::vector<double>& c);
  void			calcEquilibrium();
  void			checkElectroneutrality() const;
  //void		swapIons(unsigned i, unsigned j);
  void			correctVForPotentialDifference(double &Vwe, double &Vce, const double Awe, const double Ace);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ElectrolyteModel equations
  double		calcTransportDiffusionFactor(unsigned i, unsigned j) const;	// returns D(i,j)
  double		calcTransportMigrationFactor(unsigned i) const;				// returns z(i)*F*u(i)
  double		calcTransportConductivity() const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Electrostatics equation
  double		calcElectrostaticsConcentrationFactor(unsigned j) const;
  double		calcElectrostaticsPotentialFactor() const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solution
  void			setSolutionDensity(double solutionDensity);
  void			setSolutionKinematicViscosity(double solutionKinematicViscosity);
  void			setSolutionTemperature(double solutionTemperature);
  void			setSolutionPotential(double solutionPotential);
  double		getSolutionDensity() const;
  double		getSolutionKinematicViscosity() const;
  double		getSolutionTemperature() const;
  double		getSolutionPotential() const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solvent
  void			setSolventDielectricConstant(double solventDielectricConst);
  void			setSolventDynamicViscosity(double solventDynamicViscosity);
  double		getSolventDielectricConstant() const;
  double		getSolventDynamicViscosity() const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Ions
  void			setIonLabel(unsigned i, const std::string &label);
  void			setIonChargeNumber(unsigned i, int chargeNumber);
  void			setIonDiffusionConstant(unsigned i, double diffusionConstant);
  void			setIonDiameter(unsigned i, double diameter);
  void			setIonMolarMass(unsigned i, double molarMass);
  void			setIonConcentration(unsigned i, double concentration);
  void			setIonInletConcentration(unsigned i, double inletConcentration);
  unsigned	getNIons() const;
  std::string	getIonLabel(unsigned i) const;
  int				getIonChargeNumber(unsigned i) const;
  double		getIonDiffusionConstant(unsigned i) const;
  double		getIonDiameter(unsigned i) const;
  double		getIonMolarMass(unsigned i) const;
  double		getIonConcentration(unsigned i) const;
  double		getIonInletConcentration(unsigned i) const;
  inline double getIonTVExpansionCoefficient(unsigned i) const     {  return electrolyteSolution->getIonTVExpansionCoefficient(i);      }
  inline double getIonCDensificationCoefficient(unsigned i) const  {  return electrolyteSolution->getIonCDensificationCoefficient(i);   }
  inline double getIonMMagneticSusceptibility(unsigned i) const    {  return electrolyteSolution->getIonMMagneticSusceptibility(i);     }
  double		calcIonActivity_MM(unsigned i) const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Homogeneous reactions
  unsigned		getNHomReactions() const;
  std::string	getHomReactionLabel(unsigned r) const;
  unsigned		getHomReactionNReagents(unsigned r) const;
  unsigned		getHomReactionNProducts(unsigned r) const;
  unsigned		getHomReactionReagents(unsigned r, unsigned i) const;
  unsigned		getHomReactionProducts(unsigned r, unsigned i) const;
  int					getHomReactionStoichReag(unsigned r, unsigned i) const;
  int					getHomReactionStoichProd(unsigned r, unsigned i) const;
  double			getHomReactionForwardRateConstant(unsigned r) const;
  double			getHomReactionBackwardRateConstant(unsigned r) const;
  double			calcHomReactionRate(unsigned r) const;
  //double			calcHomReactionForwardRateDividedByCReag(unsigned r, unsigned divideCReag) const;
  //double			calcHomReactionBackwardRateDividedByCProd(unsigned r, unsigned divideCProd) const;
  double			calcHomReactionForwardRateConstant(unsigned r) const;
  double			calcHomReactionBackwardRateConstant(unsigned r) const;
  double			calcHomReactionRelativeDeviationFromEquilibrium(unsigned r) const;
  //double calcHomReactionForwardRateDividedByCReagDerivativeCReag(unsigned r, unsigned divideCReag, unsigned derivativeCReag) const;
  //double calcHomReactionBackwardRateDividedByCProdDerivativeCProd(unsigned r, unsigned divideCProd, unsigned derivativeCProd) const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Electrode reactions
  unsigned		getNElecReactions() const;
  std::string	getElecReactionLabel(unsigned r) const;
  unsigned		getElecReactionNElectrons(unsigned r) const;
  unsigned		getElecReactionNAgentsRed(unsigned r) const;
  unsigned		getElecReactionNAgentsOxi(unsigned r) const;
  unsigned		getElecReactionAgentsRed(unsigned r, unsigned i) const;
  unsigned		getElecReactionAgentsOxi(unsigned r, unsigned i) const;
  int					getElecReactionStoichRed(unsigned r, unsigned i) const;
  int					getElecReactionStoichOxi(unsigned r, unsigned i) const;
  double			calcElecReactionRate(unsigned r, double V) const;
  double			calcElecReactionCurrentDensity(unsigned r, double V) const;
  double			calcElecReactionRateDerivativeV(unsigned r, double V) const;
  double			calcElecReactionRateDerivativeU(unsigned r, double V) const;
  double			calcElecReactionRateDerivativeCRed(unsigned r, double V, unsigned i) const;
  double			calcElecReactionRateDerivativeCOxi(unsigned r, double V, unsigned i) const;
  double			calcElecReactionEquilibriumPotential(unsigned r) const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Gas reactions
  unsigned		getNGasReactions() const;
  std::string	getGasReactionLabel(unsigned r) const;
  unsigned		getGasReactionDissolvedGas(unsigned r) const;
  double			calcGasReactionRate(unsigned r) const;
  double			calcGasReactionRateDerivativeCDissGas(unsigned r) const;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void		setConductivity(double conductivity);
  //double		getConductivity() const;
  void		setElectrolyteModel(std::string EM);

protected :
  // Input filenames
  std::string		electrolyteSolutionFile;
  std::string		modelsFile;
  std::string		homReactionsFile;
  std::string		elecReactionsFile;

  ElectrolyteSolution*	electrolyteSolution;
  ElectrolyteModel*			electrolyteModel;
  Electrostatics*				electrostatics;
  unsigned				nHomReactions;
  HomReaction**		homReactions;				// array of homogenous reactions
  unsigned				nElecReactions;
  ElecReaction**	elecReactions;			// array of electrode reactions
  unsigned				nGasReactions;
  GasReaction**		gasReactions;				// array of electrode reactions

  double*			f;
  double*			x;
  double**		dfdx;
  int**				homReactionStoichMat;
  double*			cSave;
  double			conductivity;

  // Read methods
  void		setConductivityCorrectionFactor();

  // Methods
  void		assembleSystem() const;
  double	calcStepLength() const;

  // Error methods
  void		errorInvalidModel(const std::string &file, const std::string &model) const;
  void		errorConflictingData(const std::string &file1, const std::string &file2, const std::string &conflict) const;
  void		errorZero(const std::string &cppfile, const std::string &parameter) const;
};

//---------------------------------------------------------------------------


//--- SOLUTION --------------------------------------------------------------
inline void MITReM::setSolutionDensity(double solutionDensity)
{
  electrolyteSolution->setSolutionDensity(solutionDensity);
}
//---------------------------------------------------------------------------
inline void MITReM::setSolutionKinematicViscosity(double solutionKinematicViscosity)
{
  electrolyteSolution->setSolutionKinematicViscosity(solutionKinematicViscosity);
}
//---------------------------------------------------------------------------
inline void MITReM::setSolutionTemperature(double solutionTemperature)
{
  electrolyteSolution->setSolutionTemperature(solutionTemperature);
}
//---------------------------------------------------------------------------
inline void MITReM::setSolutionPotential(double solutionPotential)
{
  electrolyteSolution->setSolutionPotential(solutionPotential);
}
//---------------------------------------------------------------------------
inline double MITReM::getSolutionDensity() const
{
  return electrolyteSolution->getSolutionDensity();
}
//---------------------------------------------------------------------------
inline double MITReM::getSolutionKinematicViscosity() const
{
  return electrolyteSolution->getSolutionKinematicViscosity();
}
//---------------------------------------------------------------------------
inline double MITReM::getSolutionTemperature()const
{
  return electrolyteSolution->getSolutionTemperature();
}
//---------------------------------------------------------------------------
inline double MITReM::getSolutionPotential() const
{
  return electrolyteSolution->getSolutionPotential();
}
//---------------------------------------------------------------------------
inline void MITReM::setConductivity(double conductivity)
{
  this->conductivity = conductivity;
  setConductivityCorrectionFactor();
}
//---------------------------------------------------------------------------


//--- SOLVENT ---------------------------------------------------------------
inline void MITReM::setSolventDielectricConstant(double solventDielectricConst)
{
  electrolyteSolution->setSolventDielectricConstant(solventDielectricConst);
}
//---------------------------------------------------------------------------
inline void MITReM::setSolventDynamicViscosity(double solventDynamicViscosity)
{
  electrolyteSolution->setSolventDynamicViscosity(solventDynamicViscosity);
}
//---------------------------------------------------------------------------
inline double MITReM::getSolventDielectricConstant() const
{
  return electrolyteSolution->getSolventDielectricConstant();
}
//---------------------------------------------------------------------------
inline double MITReM::getSolventDynamicViscosity() const
{
  return electrolyteSolution->getSolventDynamicViscosity();
}
//---------------------------------------------------------------------------


//--- IONS ------------------------------------------------------------------
inline void MITReM::setIonLabel(unsigned i, const std::string &label)
{
  electrolyteSolution->setIonLabel(i,label);
}
//---------------------------------------------------------------------------
inline void MITReM::setIonChargeNumber(unsigned i, int chargeNumber)
{
  electrolyteSolution->setIonChargeNumber(i,chargeNumber);
}
//---------------------------------------------------------------------------
inline void MITReM::setIonDiffusionConstant(unsigned i, double diffusionConstant)
{
  electrolyteSolution->setIonDiffusionConstant(i,diffusionConstant);
}
//---------------------------------------------------------------------------
inline void MITReM::setIonDiameter(unsigned i, double diameter)
{
  electrolyteSolution->setIonDiameter(i,diameter);
}
//---------------------------------------------------------------------------
inline void MITReM::setIonMolarMass(unsigned i, double molarMass)
{
  electrolyteSolution->setIonMolarMass(i,molarMass);
}
//---------------------------------------------------------------------------
inline void MITReM::setIonConcentration(unsigned i, double concentration)
{
  electrolyteSolution->setIonConcentration(i,concentration);
}
//---------------------------------------------------------------------------
inline void MITReM::setIonInletConcentration(unsigned i, double inletConcentration)
{
  electrolyteSolution->setIonInletConcentration(i,inletConcentration);
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getNIons() const
{
  return electrolyteSolution->getNIons();
}
//---------------------------------------------------------------------------
inline std::string MITReM::getIonLabel(unsigned i) const
{
  return electrolyteSolution->getIonLabel(i);
}
//---------------------------------------------------------------------------
inline int MITReM::getIonChargeNumber(unsigned i) const
{
  return electrolyteSolution->getIonChargeNumber(i);
}
//---------------------------------------------------------------------------
inline double MITReM::getIonDiffusionConstant(unsigned i) const
{
  return electrolyteSolution->getIonDiffusionConstant(i);
}
//---------------------------------------------------------------------------
inline double MITReM::getIonDiameter(unsigned i) const
{
  return electrolyteSolution->getIonDiameter(i);
}
//---------------------------------------------------------------------------
inline double MITReM::getIonMolarMass(unsigned i) const
{
  return electrolyteSolution->getIonMolarMass(i);
}
//---------------------------------------------------------------------------
inline double MITReM::getIonConcentration(unsigned i) const
{
  return electrolyteSolution->getIonConcentration(i);
}
//---------------------------------------------------------------------------
inline double MITReM::getIonInletConcentration(unsigned i) const
{
  return electrolyteSolution->getIonInletConcentration(i);
}
//---------------------------------------------------------------------------
inline double MITReM::calcIonActivity_MM(unsigned i) const
{
  return electrolyteModel->calcActivityCoefficient_MM(i)*electrolyteSolution->getIonConcentration(i);
}
//---------------------------------------------------------------------------


//--- HOMOGENEOUS REACTIONS -------------------------------------------------
inline unsigned MITReM::getNHomReactions() const
{
  return nHomReactions;
}
//---------------------------------------------------------------------------
inline std::string MITReM::getHomReactionLabel(unsigned r) const
{
  return homReactions[r]->getLabel();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getHomReactionNReagents(unsigned r) const
{
  return homReactions[r]->getNReagents();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getHomReactionNProducts(unsigned r) const
{
  return homReactions[r]->getNProducts();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getHomReactionReagents(unsigned r, unsigned i) const
{
  return homReactions[r]->getReagents(i);
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getHomReactionProducts(unsigned r, unsigned i) const
{
  return homReactions[r]->getProducts(i);
}
//---------------------------------------------------------------------------
inline int MITReM::getHomReactionStoichReag(unsigned r, unsigned i) const
{
  return homReactions[r]->getStoichReag(i);
}
//---------------------------------------------------------------------------
inline int MITReM::getHomReactionStoichProd(unsigned r, unsigned i) const
{
  return homReactions[r]->getStoichProd(i);
}
//---------------------------------------------------------------------------
inline double MITReM::getHomReactionForwardRateConstant(unsigned r) const
{
  return homReactions[r]->getForwardRateConstant();
}
//---------------------------------------------------------------------------
inline double MITReM::getHomReactionBackwardRateConstant(unsigned r) const
{
  return homReactions[r]->getBackwardRateConstant();
}
//---------------------------------------------------------------------------
inline double MITReM::calcHomReactionRate(unsigned r) const
{
  return homReactions[r]->calcReactionRate();
}
//---------------------------------------------------------------------------
inline double MITReM::calcHomReactionForwardRateConstant(unsigned r) const
{
  return homReactions[r]->calcForwardRateConstant();
}
//---------------------------------------------------------------------------
inline double MITReM::calcHomReactionBackwardRateConstant(unsigned r) const
{
  return homReactions[r]->calcBackwardRateConstant();
}
//---------------------------------------------------------------------------
inline double MITReM::calcHomReactionRelativeDeviationFromEquilibrium(unsigned r) const
{
  return homReactions[r]->calcRelativeDeviationFromEquilibrium();
}
//---------------------------------------------------------------------------


//--- ELECTRODE REACTIONS ---------------------------------------------------
inline unsigned MITReM::getNElecReactions() const
{
  return nElecReactions;
}
//---------------------------------------------------------------------------
inline std::string MITReM::getElecReactionLabel(unsigned r) const
{
  return elecReactions[r]->getLabel();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getElecReactionNElectrons(unsigned r) const
{
  return elecReactions[r]->getNElectrons();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getElecReactionNAgentsRed(unsigned r) const
{
  return elecReactions[r]->getNAgentsRed();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getElecReactionNAgentsOxi(unsigned r) const
{
  return elecReactions[r]->getNAgentsOxi();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getElecReactionAgentsRed(unsigned r, unsigned i) const
{
  return elecReactions[r]->getAgentsRed(i);
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getElecReactionAgentsOxi(unsigned r, unsigned i) const
{
  return elecReactions[r]->getAgentsOxi(i);
}
//---------------------------------------------------------------------------
inline int MITReM::getElecReactionStoichRed(unsigned r, unsigned i) const
{
  return elecReactions[r]->getStoichRed(i);
}
//---------------------------------------------------------------------------
inline int MITReM::getElecReactionStoichOxi(unsigned r, unsigned i) const
{
  return elecReactions[r]->getStoichOxi(i);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionRate(unsigned r, double V) const
{
  return elecReactions[r]->calcReactionRate(V);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionCurrentDensity(unsigned r, double V) const
{
  return elecReactions[r]->calcReactionCurrentDensity(V);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionRateDerivativeV(unsigned r, double V) const
{
  return elecReactions[r]->calcReactionRateDerivativeV(V);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionRateDerivativeU(unsigned r, double V) const
{
  return elecReactions[r]->calcReactionRateDerivativeU(V);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionRateDerivativeCRed(unsigned r, double V, unsigned i) const
{
  return elecReactions[r]->calcReactionRateDerivativeCRed(V,i);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionRateDerivativeCOxi(unsigned r, double V, unsigned i) const
{
  return elecReactions[r]->calcReactionRateDerivativeCOxi(V,i);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElecReactionEquilibriumPotential(unsigned r) const
{
  return elecReactions[r]->calcEquilibriumPotential();
}
//---------------------------------------------------------------------------


//--- GAS REACTIONS ---------------------------------------------------
inline unsigned MITReM::getNGasReactions() const
{
  return nGasReactions;
}
//---------------------------------------------------------------------------
inline std::string MITReM::getGasReactionLabel(unsigned r) const
{
  return gasReactions[r]->getLabel();
}
//---------------------------------------------------------------------------
inline unsigned MITReM::getGasReactionDissolvedGas(unsigned r) const
{
  return gasReactions[r]->getDissolvedGas();
}
//---------------------------------------------------------------------------
inline double MITReM::calcGasReactionRate(unsigned r) const
{
  return gasReactions[r]->calcReactionRate();
}
//---------------------------------------------------------------------------
inline double MITReM::calcGasReactionRateDerivativeCDissGas(unsigned r) const
{
  return gasReactions[r]->calcReactionRateDerivativeCDissGas();
}
//---------------------------------------------------------------------------


//--- TRANSPORT EQUATIONS ---------------------------------------------------
inline double MITReM::calcTransportDiffusionFactor(unsigned i, unsigned j) const
{
  return electrolyteModel->calcDiffusionFactor(i,j);
}
//---------------------------------------------------------------------------
inline double MITReM::calcTransportMigrationFactor(unsigned i) const
{
  return electrolyteModel->calcMigrationFactor(i);
}
//---------------------------------------------------------------------------
inline double MITReM::calcTransportConductivity() const
{
  return electrolyteModel->calcConductivity();
}
//---------------------------------------------------------------------------


//--- ELECTROSTATICS EQUATION -----------------------------------------------
inline double MITReM::calcElectrostaticsConcentrationFactor(unsigned j) const
{
  return electrostatics->calcConcentrationFactor(j);
}
//---------------------------------------------------------------------------
inline double MITReM::calcElectrostaticsPotentialFactor() const
{
  return electrostatics->calcPotentialFactor();
}
//---------------------------------------------------------------------------

#endif



