//---------------------------------------------------------------------------

#ifndef HomReactionTermH
#define HomReactionTermH

//---------------------------------------------------------------------------

#include "ElementTerm.h"

//---------------------------------------------------------------------------

class HomReactionTerm : public ElementTerm
{
public :
	HomReactionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_);
	virtual ~HomReactionTerm();

	virtual void	calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};
	virtual void	calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

	virtual void	calcIonCurrentDensities(EmptyEmptyDoubleVectorList elementCurr, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors) {};

protected:
	unsigned		nHomReactions;
	double**		kf;
	double**		kb;
	double*			Hfmj;
	double*			Hbmj;
	double**		Hfmjk;
	double**		Hbmjk;
};

//---------------------------------------------------------------------------

#endif

