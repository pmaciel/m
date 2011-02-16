
#include "common.h"


/*
 * This file contains MPI functionality for parallel
 * computation of dispersed two-phase flow.
 *
 * MPI Allreduce of the data type MPI_INT with MPI_SUM as an
 * operator, extended to arrays.
 */

void plas_MpiAllSumIntArray(int *val, int size)
{
#ifdef MPI
  int i,sendval[size],recval[size];
  for(i=0; i<size; i++){sendval[i] = val[i];}
  MPI_Allreduce(sendval,recval,size,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  for(i=0; i<size; i++){val[i] = recval[i];}
#else
  return;
#endif
}

