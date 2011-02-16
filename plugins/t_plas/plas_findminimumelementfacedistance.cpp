
#include "common.h"


/*
 * This file includes all functinality to perform an element
 * search for a dispersed entity.
 *
 * This function finds the smallest element face distance to
 * an entity.
 */

void plas_FindMinimumElementFaceDistance(int numDim, LOCAL_ENTITY_VARIABLES *ent, int *idx, double *dmin)
{
  int idim,jdim;
  double idist,posvec[numDim],normvec[numDim];

  for(idim=0; idim<ent->edata.numElmFaces; idim++){
    for(jdim=0; jdim<numDim; jdim++){
      posvec[jdim] = ent->pos[jdim]-ent->edata.elmFaceVectors[idim][jdim];
      normvec[jdim] = ent->edata.elmNorms[idim][jdim];
    }
    plas_NormalizeVector(numDim,normvec);
    idist = plas_CalcVectScalarProduct(numDim,posvec,normvec);
    if(idim==0){
      *idx = idim;
      *dmin = idist;
    } else{
      if(idist<*dmin){
        *idx = idim;
        *dmin = idist;
      }
    }
  }
}

