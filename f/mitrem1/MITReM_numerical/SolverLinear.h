//---------------------------------------------------------------------------

#ifndef SolverLinearH
#define SolverLinearH

//---------------------------------------------------------------------------

#include <vector>
#include "Element.h"

//---------------------------------------------------------------------------

class SolverLinear
{
public :
  SolverLinear (const Element* const* elements, unsigned nElements, unsigned nNodes, unsigned nVariables);
  virtual ~SolverLinear () {};

  virtual void init () = 0;
  virtual void add (unsigned row, unsigned col, double value) = 0;
  virtual void multiplyWithVector (double* X, double* Result) = 0;
  virtual void residu (double* X, double* B) = 0;
  virtual double getResidual (double* X, double* B, unsigned int row) = 0;
  virtual void imposeVar (unsigned row) = 0;
  virtual void solve (double* B) = 0;
  virtual void print (const std::string &filename) = 0;
  virtual void clearRow (unsigned row) = 0;
  virtual void linearCombination (unsigned toAddTo, unsigned toAdd, double factor) = 0;

protected :
  unsigned  nRows;
  std::vector< std::vector<unsigned> >  matStruct;
};

//---------------------------------------------------------------------------

class SolverLinear_BandLU : public SolverLinear
{
public :
  SolverLinear_BandLU (const Element* const* elements, unsigned nElements, unsigned nNodes, unsigned nVariables);
  virtual ~SolverLinear_BandLU();

  virtual void init ();
  virtual void add (unsigned row, unsigned col, double value);
  virtual void multiplyWithVector (double* X, double* Result);
  virtual void residu (double* X, double* B);
  virtual double getResidual (double* X, double* B, unsigned int row);
  virtual void imposeVar (unsigned row);
  virtual void solve (double* B);
  virtual void print (const std::string &filename);
  virtual void clearRow (unsigned row);
  virtual void linearCombination (unsigned toAddTo, unsigned toAdd, double factor);

private :
  void LUDecomp (double* B);
  void errorSingularMatrix();

  unsigned  ld, ud, LD, UD, nColumns;
  double**  A;
};

//---------------------------------------------------------------------------

class SolverLinear_GMRES : public SolverLinear
{
public :
  SolverLinear_GMRES (const Element* const* elements, unsigned nElements, unsigned nNodes, unsigned nVariables);
  virtual ~SolverLinear_GMRES();

  virtual void init ();
  virtual void add (unsigned row, unsigned col, double value);
  virtual void multiplyWithVector (double* X, double* Result);
  virtual void residu (double* X, double* B);
  virtual double getResidual (double* X, double* B, unsigned int row);
  virtual void imposeVar (unsigned row);
  virtual void solve (double* B);
  virtual void print (const std::string &filename);
  virtual void clearRow (unsigned row);
  virtual void linearCombination (unsigned toAddTo, unsigned toAdd, double factor);

private :
  int iluk(int *n, double *a, int *ja, int *ia, int *lfil,
            double*& aluold, int*& jluold, int *ju, int*& levsold, int *iwk,
            double *w, int *jw, int *ierr);
  double ddot(int *n, double *dx, int *incx, double *dy, int *incy);
  int daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
  double dnrm2(int *n, double *dx, int *incx);
  void amux(int *n, double *x, double *y, double *a, int *ja, int *ia);
  void lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju);
  void pgmres(int *n, int *im, double *rhs, double *sol, double *vv,
      double *eps, int *maxits, int*iout,
      double *aa, int *ja, int *ia,
      double *alu, int *jlu, int *ju,
      int *ierr);
  unsigned getIndex(unsigned row, unsigned col);

  int                               c__1;
  unsigned                          N,nVariables;
  int                          *Rows, *Columns;
  double*                           A;
  void*  pt_[64];
  void solveWithPardiso(double B[]);
  void solveWithGMRES(double B[]);
};

//---------------------------------------------------------------------------

#endif
