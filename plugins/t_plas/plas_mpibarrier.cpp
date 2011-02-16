
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * MPI Barrier
 */

void plas_MpiBarrier()
{
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#else
  return;
#endif
}

