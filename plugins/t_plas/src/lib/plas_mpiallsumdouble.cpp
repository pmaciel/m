
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * MPI Allreduce of the data type MPI_DOUBLE with MPI_SUM as
 * an operator, applicable for scalar data.
 */

double plas_MpiAllSumDouble(double val)
{
#ifdef MPI
  double recval;
  MPI_Allreduce(&val,&recval,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return recval;
#else
  return val;
#endif
}

