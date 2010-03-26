//---------------------------------------------------------------------------

#include "Term.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Term::Term(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_)
  : nDimensions(nDimensions_), nNodes(nNodes_), nVariables(nVariables_), mitrem(mitrem_)
{
  nIons = mitrem->getNIons();
  normals = new EmptyDoubleVector[nNodes];
  for (unsigned m=0; m<nNodes; m++)
  {
    normals[m] = new double[nDimensions];
  }
}
//---------------------------------------------------------------------------
Term::~Term()
{
  for (unsigned m=0; m<nNodes; m++)
  {
    delete[] normals[m];
  }
  delete[] normals;
}
//---------------------------------------------------------------------------

