
#include "common.h"


/*
 * This file includes the functionality to perform  trajectory
 * integrations of dispersed entities.
 *
 * This function solves the trajectory equation.
 */

void plas_SolveGaussSeidel(int numDim, double **mat, double *s, double *rhs)
{
  int idim,jdim;
  double errLoc,errMax,sOld[numDim];

  for(idim=0; idim<numDim; idim++){
    sOld[idim] = s[idim];
  }

  do{
    errMax = 0.0;
    for(idim=0; idim<numDim; idim++){
      s[idim] = 0.0;
      for(jdim=0; jdim<numDim; jdim++){
        if(jdim<idim){
          s[idim] += mat[idim][jdim]*s[jdim];
        }
        if(jdim>idim){
          s[idim] += mat[idim][jdim]*sOld[jdim];
        }
      }
      s[idim] = (rhs[idim]-s[idim])/mat[idim][idim];
      errLoc = fabs(s[idim]-sOld[idim]);
      if(errLoc>errMax){errMax = errLoc;}
      sOld[idim] = s[idim];
    }
  }while(errMax>1e-6);
}

