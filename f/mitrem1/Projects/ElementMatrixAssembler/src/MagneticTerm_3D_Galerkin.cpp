//---------------------------------------------------------------------------

#include "MagneticTerm_3D_Galerkin.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MagneticTerm_3D_Galerkin::MagneticTerm_3D_Galerkin(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
  : MagneticTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MagneticTerm_3D_Galerkin::~MagneticTerm_3D_Galerkin()
{
}
//---------------------------------------------------------------------------
void MagneticTerm_3D_Galerkin::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  // !!! STILL HAS TO BE IMPLEMENTED !!!
}
//---------------------------------------------------------------------------
void MagneticTerm_3D_Galerkin::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
}
//---------------------------------------------------------------------------

