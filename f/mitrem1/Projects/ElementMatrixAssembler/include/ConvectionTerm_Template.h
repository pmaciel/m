//---------------------------------------------------------------------------

#ifndef ConvectionTerm_TemplateH
#define ConvectionTerm_TemplateH

//---------------------------------------------------------------------------

#include <iostream>
#include "ConvectionTerm.h"

//---------------------------------------------------------------------------

enum ListConvectionSchemes {
  Scheme_Null,
  Scheme_Galerkin,
  Scheme_LDA,
  Scheme_N
};

//---------------------------------------------------------------------------

template< unsigned DIM, ListConvectionSchemes CSCHEME >
class ConvectionTerm_Template : public ConvectionTerm
{
 public:

  ConvectionTerm_Template(unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_)
    : ConvectionTerm(DIM,DIM+1,nVariables_,mitrem_,elementProps_)
  {
  }

  ~ConvectionTerm_Template()
  {
  }

  void calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
  {
    std::cerr << "convection scheme not implemented!" << std::endl;
    throw 42;
  }

  void calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
  {
  }

 private:

  void calculate_k(DoubleVectorList velocities, DoubleVectorList coordinates)
  {
    m_kp.assign(nNodes,0.);
    m_km.assign(nNodes,0.);
    for (unsigned d=0; d<nDimensions; ++d) {
      averageVelocity[d] = 0.;
      for (unsigned m=0; m<nNodes; ++m)
        averageVelocity[d] += velocities[m][d];
      averageVelocity[d] /= (double) nNodes;
    }
    m_sumkp = 0.;
    for (unsigned m=0; m<nNodes; ++m) {
      normal = elementProps->calcNormal(m,coordinates);
      k[m] = 0.;
      for (unsigned d=0; d<nDimensions; ++d)
        k[m] += normal[d]*averageVelocity[d];
      k[m] /= (double) nDimensions;
      m_kp[m] = std::max(0.,k[m]);
      m_km[m] = std::min(0.,k[m]);
      m_sumkp += m_kp[m];
    }
  }

  double calculate_factor(DoubleList volumeGasFractions)
  {
    double f = 0.;
    for (unsigned m=0; m<nNodes; ++m)
      f += volumeGasFractions[m];
    return (1. - f/(double) nNodes);
  }

  std::vector< double > m_k;      // k+
  std::vector< double > m_kp;     // k+
  std::vector< double > m_km;     // k-
  double                m_sumkp;  // sum(k+)
};

//---------------------------------------------------------------------------

template<>
void ConvectionTerm_Template< 2,Scheme_Null >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
};

template<>
void ConvectionTerm_Template< 3,Scheme_Null >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
};

//---------------------------------------------------------------------------

template<>
void ConvectionTerm_Template< 2,Scheme_Galerkin >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  calculate_k(velocities,coordinates);
  const double F = calculate_factor(volumeGasFractions);
  const double W = 1./((double) nNodes);
  for (unsigned i=0; i<nIons; ++i)
    for (unsigned m=0; m<nNodes; ++m)
      for (unsigned n=0; n<nNodes; ++n)
        elementMat[eq(m,i)][var(n,i)] -= W*k[n]*F;
};

template<>
void ConvectionTerm_Template< 3,Scheme_Galerkin >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  calculate_k(velocities,coordinates);
  const double F = calculate_factor(volumeGasFractions);
  const double W = 1./((double) nNodes);
  for (unsigned i=0; i<nIons; ++i)
    for (unsigned m=0; m<nNodes; ++m)
      for (unsigned n=0; n<nNodes; ++n)
        elementMat[eq(m,i)][var(n,i)] -= W*k[n]*F;
};

//---------------------------------------------------------------------------

template<>
void ConvectionTerm_Template< 2,Scheme_LDA >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  calculate_k(velocities,coordinates);
  const double F = calculate_factor(volumeGasFractions);
  for (unsigned i=0; i<nIons; ++i)
    for (unsigned m=0; m<nNodes; ++m)
      for (unsigned n=0; n<nNodes; ++n)
        elementMat[eq(m,i)][var(n,i)] -= m_kp[m]*k[n]/m_sumkp*F;
};

template<>
void ConvectionTerm_Template< 3,Scheme_LDA >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  calculate_k(velocities,coordinates);
  const double F = calculate_factor(volumeGasFractions);
  for (unsigned i=0; i<nIons; ++i)
    for (unsigned m=0; m<nNodes; ++m)
      for (unsigned n=0; n<nNodes; ++n)
        elementMat[eq(m,i)][var(n,i)] -= m_kp[m]*k[n]/m_sumkp*F;
};

//---------------------------------------------------------------------------

template<>
void ConvectionTerm_Template< 2,Scheme_N >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  calculate_k(velocities,coordinates);
  const double F = calculate_factor(volumeGasFractions);
  for (unsigned i=0; i<nIons; ++i)
    for (unsigned m=0; m<nNodes; ++m)
      for (unsigned n=0; n<nNodes; ++n)
        elementMat[eq(m,i)][var(n,i)] -= m_kp[m]*(m==n? 1:m_km[n]/m_sumkp)*F;
};

template<>
void ConvectionTerm_Template< 3,Scheme_N >::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList volumeGasFractions, DoubleVectorList magneticFieldVectors)
{
  calculate_k(velocities,coordinates);
  const double F = calculate_factor(volumeGasFractions);
  for (unsigned i=0; i<nIons; ++i)
    for (unsigned m=0; m<nNodes; ++m)
      for (unsigned n=0; n<nNodes; ++n)
        elementMat[eq(m,i)][var(n,i)] -= m_kp[m]*(m==n? 1:m_km[n]/m_sumkp)*F;
};

//---------------------------------------------------------------------------

#endif

