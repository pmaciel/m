
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * MPI Allreduce of the data type MPI_INT with MPI_SUM as an
 * operator, applicable for scalar data.
 */

int plas_MpiAllSumInt(int val)
{
#ifdef MPI
  int recval;
  MPI_Allreduce(&val,&recval,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return recval;
#else
  return val;
#endif
}

