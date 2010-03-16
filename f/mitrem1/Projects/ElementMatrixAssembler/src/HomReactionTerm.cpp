//---------------------------------------------------------------------------

#include "HomReactionTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
HomReactionTerm::HomReactionTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: ElementTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
	nHomReactions = mitrem->getNHomReactions();
	kf = new double*[nNodes];
	kb = new double*[nNodes];
	Hfmj = new double[nNodes];
	Hbmj = new double[nNodes];
	Hfmjk = new double*[nNodes];
	Hbmjk = new double*[nNodes];
	for (unsigned n=0; n<nNodes; n++)
	{
		kf[n] = new double[nHomReactions];
		kb[n] = new double[nHomReactions];
		Hfmjk[n] = new double[nNodes];
		Hbmjk[n] = new double[nNodes];
	}
}
//---------------------------------------------------------------------------
HomReactionTerm::~HomReactionTerm()
{
	for (unsigned n=0; n<nNodes; n++)
	{
		delete[] kf[n];
		delete[] kb[n];
		delete[] Hfmjk[n];
		delete[] Hbmjk[n];
	}
	delete[] kf;
	delete[] kb;
	delete[] Hfmj;
	delete[] Hbmj;
	delete[] Hfmjk;
	delete[] Hbmjk;
}
//---------------------------------------------------------------------------

