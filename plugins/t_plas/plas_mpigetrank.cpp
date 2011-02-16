
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * This function returns the rank of the calling process.
 */

int plas_MpiGetRank()
{
#ifdef MPI
  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  return irank;
#else
  return 0;
#endif
}

