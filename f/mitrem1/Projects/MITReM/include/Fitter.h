//---------------------------------------------------------------------------

#ifndef FitterH
#define FitterH

//---------------------------------------------------------------------------

#include <string>
#include <vector>
#include "Electrolyte.h"
#include "IonFit.h"
#include "HomReactionFit.h"

//---------------------------------------------------------------------------

class Fitter
{
public :
	Fitter(const std::string &name_);
	~Fitter();

	// Methods
	void		solve();

protected :
	std::string		name;

	std::string		electrolytesFile;
	std::string		ionsFile;
	std::string		homReactionsFile;

	unsigned			nElectrolytes;
	Electrolyte**		electrolytes;
	unsigned			nIons;
	IonFit**			ions;
	unsigned			nHomReactions;
	HomReactionFit**	homReactions;

	unsigned				nDiametersToFit, nDiffusionConstantsToFit, nEquilibriumConstantsToFit;
	std::vector<unsigned>	diametersToFit;
	std::vector<unsigned>	diffusionConstantsToFit;
	std::vector<unsigned>	equilibriumConstantsToFit;

	unsigned			nExperimentalValues,nParametersToFit;
	double*				R;
	double*				JtR;
	double**			J;
	double**			JtJ;
	double				maxStep,residu;
	unsigned			maxIter;

	// Read methods
	void		readElectrolytes();
	void		readIons();
	void		readHomReactions();

	// Methods
	unsigned	calcRowOsmoticCoefficient(unsigned e, unsigned t) const;	
	unsigned	calcRowActivityCoefficient(unsigned e, unsigned t) const;
	unsigned	calcRowConductivity(unsigned e, unsigned t) const;
	unsigned	calcRowTransportNumber(unsigned e, unsigned t) const;
	unsigned	calcRowDiffusionCoefficient(unsigned e, unsigned t) const;
	double		calcLimitingOsmoticCoefficient(unsigned e) const;	
	double		calcLimitingActivityCoefficient(unsigned e) const;
	double		calcLimitingConductivity(unsigned e) const;
	double		calcLimitingTransportNumber(unsigned e) const;
	double		calcLimitingDiffusionCoefficient(unsigned e) const;
	double		calcStepLength() const;
			
	// Error methods
	void		errorInvalidModel(const std::string &file, const std::string &model) const;
	void		errorConflictingData(const std::string &file1, const std::string &file2, const std::string &conflict) const;
	void		errorZero(const std::string &cppfile, const std::string &parameter) const;
	void		errorSingularMatrix() const;
};

//---------------------------------------------------------------------------

#endif

