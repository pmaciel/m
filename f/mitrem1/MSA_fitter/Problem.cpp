//---------------------------------------------------------------------------

#include "Problem.h"
#include "DatFileReader.h"

//---------------------------------------------------------------------------


//--- READ METHODS ----------------------------------------------------------
Problem::Problem ()
{
  DatFileReader datFile("MSA_fitter.dat");

  datFile.readScalar("[fitterName]",fitterName);

  fitter = new Fitter(fitterName);
  fitter->solve();
}
//---------------------------------------------------------------------------
Problem::~Problem ()
{
  delete fitter;
}
//---------------------------------------------------------------------------
