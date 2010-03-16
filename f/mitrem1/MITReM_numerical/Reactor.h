//---------------------------------------------------------------------------

#ifndef ReactorH
#define ReactorH

//---------------------------------------------------------------------------

#include <string>
#include <fstream>
#include "MITReM.h"
#include "Vector.h"
#include "Element.h"
#include "Boundary.h"
#include "Electrode.h"
#include "ElementMatrixAssembler.h"
#include "SolverLinear.h"

//---------------------------------------------------------------------------

class Reactor {
public :

	Reactor(const std::string &name, const std::string &database, const std::string &ecLabel, const std::string &outputFile, const std::string &outputFormat, bool DoE, const std::string &initialConditionsFile);
	virtual ~Reactor();

	void			solve();


	void		setIonDiffusionConstant(unsigned i, double diffusionConstant) 
	{mitrem->setIonDiffusionConstant(i, diffusionConstant);};

	void		setIonDiameter(unsigned i, double diameter)
	{mitrem->setIonDiameter(i, diameter);};

	void		setIonMolarMass(unsigned i, double molarMass)
	{mitrem->setIonMolarMass(i, molarMass);};

	void		setIonConcentration(unsigned i, double concentration)
	{mitrem->setIonConcentration(i, concentration); mitrem->setIonInletConcentration(i, concentration);};

	void		setConductivity(double conductivity)
	{mitrem->calcEquilibrium(); mitrem->setConductivity(conductivity);};

	int			getIonChargeNumber(unsigned i) const
	{return mitrem->getIonChargeNumber(i);};

	double		getConductivity() const
	{mitrem->calcEquilibrium(); return mitrem->calcTransportConductivity();}

	void		setElectrolyteModel(std::string EM)
	{mitrem->setElectrolyteModel(EM);};

	void		init();


	double*		xVec;
	double***				ionCurrentDensities;
	double*****				gasReactionRates;
	double*****				elecReactionCurrentDensities;
	double		residu;

	std::string		outputFile;

protected :
	bool DoE;

	// Input filenames
	std::string reactorFile;
	std::string nodesFile;
	std::string elementsFile;
	std::string boundaryElementsFile;
	std::string boundariesFile;
	std::string electrodesFile;
	std::string flowFieldFile;
	std::string magneticFieldFile;
	std::string solverFile;
	std::string electrodePotentialsFile;
	std::string miotrasDatFile;
	std::string miotrasFlowFile;
	std::string initialConditionsFile;

	// Input file
	std::ifstream	input;

	// Output file
	//std::string		outputFile;
	std::ofstream	output;
	std::string		outputFormat;

	MITReM*					mitrem;
	ElementMatrixAssembler*	elementMatrixAssembler;
	SolverLinear*			AMat;
	
	unsigned		nDimensions, nElements, nBoundaryElements, nBoundaries, nElectrodes, nNodes, nIons, nVariables, nElementNodes, nBoundaryElementNodes, size;
	std::string		dimensions;
	Element**		elements;				// array of elements
	Element**		boundaryElements;		// array of boundary elements
	Boundary**		boundaries;				// array of boundaries
	Electrode**		electrodes;				// array of electrodes
	Vector**		nodes;					// array of nodes
	Vector**		flowField;				// array of velocity vectors
	Vector**		magneticField;				// array of magnetic field vectors
	
	bool			timeAccurate;
	unsigned		nElectrodePotentials;
	class ElectrodePotential
	{
	public:
		double	t;
		double*	V;
	};
	ElectrodePotential*	electrodePotentials;

	double		/*residu, */residuMax, relaxFacConc, relaxFacConcInc, relaxFacPot, relaxFacPotInc;
	unsigned	iterMax;
	bool		converged;
	double*		inletVec;
	//double*		xVec;
	double*		xVecOld;
	double*		xVecPtr;
	double*		bVec;
	double*		bVecOld;
	double		dt, dtOld, timeFactor1, timeFactor2, timeFactor3;

	double*		dualMesh;

	EmptyDoubleVectorList	coordinates;
	EmptyDoubleVectorList	velocities;
	EmptyDoubleVectorList	magneticFieldVectors;
	EmptyDoubleVectorList	concentrations;
	EmptyDoubleList			potentials;
	EmptyDoubleList			temperatures; 
	EmptyDoubleList			densities; 
	EmptyDoubleList			voidFractions;
	double*					currents;

	// Read methods
	void readReactor();
	void readNodes();
	void readElements();
	void readBoundaryElements();
	void readBoundaries();
	void readElectrodes();
	void readFlowField();
	void readMagneticField();
	void readSolver();
	void readElectrodePotentials();

	// Error methods
	void errorConflictingData(const std::string &file1, const std::string &file2, const std::string &conflict);
	void errorFileDoesNotExist(const std::string &file);

	// Solver methods
	unsigned		var(unsigned n, unsigned j) const;
	unsigned		eq(unsigned m, unsigned i) const;
	void			solve_Stationary(unsigned p);
	void			solve_TimeAccurate(unsigned p);
	void			write_Stationary_Excel(unsigned p);
	void			write_TimeAccurate_Excel(unsigned p);
	void			write_TecPlot(unsigned p);
	void			write_TecPlot_Residual(unsigned p);
	void			initLinearSystem_Stationary();
	void			initLinearSystem_TimeAccurate();
	void			solveLinearSystem_Stationary(unsigned iter);
	void			solveLinearSystem_TimeAccurate();
	void			addElementMatrices_Stationary();
	void			addElementMatrices_TimeAccurate();
	void			addElementJacobians_Stationary();
	void			addElementJacobians_TimeAccurate();
	void			addBoundaryElementVectors();
	void			addBoundaryElementJacobians();
	void			addElementTimeMatrices();
	void			imposeEssentialBoundaryConditions();
	void			calcCurrentDensities();

	void			skip(std::ifstream &istrm, unsigned nSkips, char delimiter);

};

//---------------------------------------------------------------------------

inline unsigned Reactor::var(unsigned n, unsigned j) const
{
	return n*nVariables+j;
}

//---------------------------------------------------------------------------

inline unsigned Reactor::eq(unsigned m, unsigned i) const
{
	//if (i == 0) return m*nVariables+nIons;
	//else if (i == nIons) return m*nVariables;
	//else return m*nVariables+i;
	return m*nVariables+i;
}
//---------------------------------------------------------------------------

#endif