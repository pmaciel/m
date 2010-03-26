//---------------------------------------------------------------------------

#ifndef SolverBoundaryElementH
#define SolverBoundaryElementH

//---------------------------------------------------------------------------

#include "MatrixStructure.h"

//---------------------------------------------------------------------------

class SolverBoundaryElement : public MatrixStructure
{
public :
  SolverBoundaryElement(MITReM* mitrem_, Reactor* reactor_);
  virtual ~SolverBoundaryElement();

  double* calcBoundaryElementVec(unsigned be, unsigned e, double** boundaryElementVariables_);
  double** calcBoundaryElementJac(unsigned be, unsigned e, double** boundaryElementVariables_);

protected:
  virtual void calcElecReactionVec() { };
  virtual void calcElecReactionJac() { };

  unsigned size;
  double* boundaryElementVec;
  double** boundaryElementJac;
  double** boundaryElementVariables;
  Electrode* electrode;
};

//---------------------------------------------------------------------------

class SolverBoundaryElement_1D : public SolverBoundaryElement
{
public :
  SolverBoundaryElement_1D(MITReM* mitrem_, Reactor* reactor_) : SolverBoundaryElement(mitrem_, reactor_) { };
  virtual ~SolverBoundaryElement_1D() { };

protected:
  virtual void calcElecReactionVec();
  virtual void calcElecReactionJac();
};

//---------------------------------------------------------------------------

#endif
