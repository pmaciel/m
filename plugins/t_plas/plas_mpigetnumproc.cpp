
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * This function returns the number of processes.
 */

int plas_MpiGetNumProc()
{
#ifdef MPI
  int numProc;
  MPI_Comm_size(MPI_COMM_WORLD,&numProc);
  return numProc;
#else
  return 1;
#endif
}

