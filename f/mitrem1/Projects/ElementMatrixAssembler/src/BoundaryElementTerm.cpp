//---------------------------------------------------------------------------

#include "BoundaryElementTerm.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
BoundaryElementTerm::BoundaryElementTerm(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, BoundaryElementProps* boundaryElementProps_)
  : Term(nDimensions_, nNodes_, nVariables_, mitrem_), boundaryElementProps(boundaryElementProps_)
{
}
//---------------------------------------------------------------------------
BoundaryElementTerm::~BoundaryElementTerm()
{
}
//---------------------------------------------------------------------------

