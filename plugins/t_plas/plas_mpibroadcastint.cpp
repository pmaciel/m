
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * MPI Broadcast of the data type MPI_INT
 */

void plas_MpiBroadcastInt(int *variable, int size, int root)
{
#ifdef MPI
  MPI_Bcast(variable,size,MPI_INT,root,MPI_COMM_WORLD);
#else
  return;
#endif
}

