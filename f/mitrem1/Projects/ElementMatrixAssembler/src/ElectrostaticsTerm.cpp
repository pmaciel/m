//---------------------------------------------------------------------------

#include "ElectrostaticsTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElectrostaticsTerm::ElectrostaticsTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : ElementTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
  Z = new double[nIons];
}
//---------------------------------------------------------------------------
ElectrostaticsTerm::~ElectrostaticsTerm()
{
  delete[] Z;
}
//---------------------------------------------------------------------------

