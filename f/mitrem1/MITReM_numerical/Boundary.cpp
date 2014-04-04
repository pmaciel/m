//---------------------------------------------------------------------------

#include "Boundary.h"

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Boundary::Boundary (unsigned nBoundaryElements_)
  : nBoundaryElements(nBoundaryElements_)
{
  boundaryElements = new unsigned[nBoundaryElements];
}
//---------------------------------------------------------------------------
Boundary::~Boundary ()
{
  delete[] boundaryElements;
}
//---------------------------------------------------------------------------
void Boundary::setBoundaryElements(unsigned b, unsigned boundaryElement)
{
  if (b < nBoundaryElements)
  {
    boundaryElements[b] = boundaryElement;
  }
  else
  {
    // error message
  }
}
//---------------------------------------------------------------------------
