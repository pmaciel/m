//---------------------------------------------------------------------------

#include "ElementTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
ElementTerm::ElementTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
	: Term(nDimensions_, nNodes_, nVariables_, mitrem_), elementProps(elementProps_)
{
}
//---------------------------------------------------------------------------
ElementTerm::~ElementTerm() 
{
}
//---------------------------------------------------------------------------

