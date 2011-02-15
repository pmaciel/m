
#include "common.h"

/*
 * This function computes coalescense of bubbles according
 * to the thin film theory.
 *
 * The routine is not used for the time being. Before using
 * it, please check carefully the implementation.
 */

int plas_Coalescence(PLAS_DATA *data, double dj, double di, double *uijRelPrPr, double *x, double *y, double *z,
  double Mi, double Mj, double *uiPrPr, double *uiPrPrNew, double *uiNew, double *ui)
{
  double h0 = 1e-4;
  double hf = 1e-8;
  double Cc = 0.5;
  double sigma = 0.06;
  double Rij = 2.0*1/(2/dj+2/di);
  double T = sqrt(Rij*Rij*Rij*data->fp.rhoCont/(16*sigma))*log(h0/hf);
  double Tau_ij = Cc*Rij/uijRelPrPr[0];
  int idim,isCoal;
  double diStar;

  //***Coalescense model***//

  if(Tau_ij>T){

    isCoal = 1;
    data->sd.coalesc++;

    //***Compute new entity diameter***//

    diStar = pow((dj*dj*dj+di*di*di),1/3);

    if(data->fp.numDim==2){

      uiPrPrNew[0] = uiPrPr[0]+Mi/(Mi+Mj);
      uiPrPrNew[1] = uiPrPr[1];

      uiNew[0] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,x);
      uiNew[1] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,y);

    } else if(data->fp.numDim==3){

      uiPrPrNew[0] = uiPrPr[0]+Mi/(Mi+Mj);
      uiPrPrNew[1] = uiPrPr[1];
      uiPrPrNew[2] = uiPrPr[2];

      uiNew[0] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,x);
      uiNew[1] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,y);
      uiNew[2] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,z);
    }

    di = diStar;

    //***Update entity velocity***//

    for(idim=0; idim<data->fp.numDim; idim++){
      ui[idim] = uiNew[idim];
    }

  } else{
    isCoal = 0;
  }

  return isCoal;
}

