
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * MPI Broadcast of the data type MPI_DOUBLE
 */

void plas_MpiBroadcastDouble(double *variable, int size, int root)
{
#ifdef MPI
  MPI_Bcast(variable,size,MPI_DOUBLE,root,MPI_COMM_WORLD);
#else
  return;
#endif
}

