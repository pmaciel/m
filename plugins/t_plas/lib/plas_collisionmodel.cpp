
#include "common.h"

/*
 * This function computes inter-entity collisions based on the
 * stochastic model of Sommerfeld.
 */

void plas_CollisionModel(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, int numDens, double dtLagr)
{
  int idim;
  double
    c_T = 0.3,
    e = 1.0,
    mu_static = 0.0,
    mu_dynamic = 0.0,
    ui[data->fp.numDim],
    uiFluct[data->fp.numDim],
    uj[data->fp.numDim],
    ujFluct[data->fp.numDim],
    uiNew[data->fp.numDim],
    uijRel[data->fp.numDim],
    Pcoll,
    FreqColl,
    conc,
    di,
    dj,
    tau_f,
    R,
    stokes,
    epsilon,
    rms,
    pos1[data->fp.numDim],
    pos1Pr[data->fp.numDim],
    pos2[data->fp.numDim],
    pos2Pr[data->fp.numDim],
    posRel[data->fp.numDim],
    x[data->fp.numDim],
    y[data->fp.numDim],
    z[data->fp.numDim],
    xPr[data->fp.numDim],
    yPr[data->fp.numDim],
    zPr[data->fp.numDim],
    xPrPr[data->fp.numDim],
    yPrPr[data->fp.numDim],
    zPrPr[data->fp.numDim],
    uiPrPr[data->fp.numDim],
    uiPrPrNew[data->fp.numDim],
    ujPrPr[data->fp.numDim],
    uijRelPrPr[data->fp.numDim],
    Phi,
    Psi,
    Mi,
    Mj,
    Jx,
    Jy,
    Jz,
    Y = 0.,
    Z = 0.,
    L = 0.,
    length;


  /*
  double **M = new double*[data->fp.numDim];
  for (int r=0; r<data->fp.numDim; ++r)
    M[r] = new double[data->fp.numDim];
  */
  double **M = (double **) malloc( sizeof(double*)*(data->fp.numDim) );
  for (idim=0; idim<data->fp.numDim; ++idim)
    M[idim] = (double *) malloc( sizeof(double)*(data->fp.numDim) );


  //***Set local entity properties (subscript i)***//

  for(idim=0; idim<data->fp.numDim; idim++){
    ui[idim] = ent->vel[idim];
  }
  di = ent->diam;

  //***Set fictitious collision partner properties (subscript j)***//

  for(idim=0; idim<data->fp.numDim; idim++){
    uj[idim] = plas_RandomGaussian(data->pd[ent->node].avgVel[idim],data->pd[ent->node].stdVel[idim]);
  }
  dj = plas_RandomGaussian(data->pd[ent->node].avgDiam,data->pd[ent->node].stdDiam);

  //***Correlation between the fluctuations of the real and fictitious particle***//

  if(data->ip.collisionModel==2){

    for(idim=0; idim<data->fp.numDim; idim++){
      uiFluct[idim] = ui[idim]-data->pd[ent->node].avgVel[idim];
    }

    tau_f = c_T*plasinterface_getEulerianTimeScale(ent->node);
    stokes = ent->kinRespTime/tau_f;
    epsilon = plas_RandomGaussian(0.0,1.0);
    R = exp(-0.55*pow(stokes,0.4));

    for(idim=0; idim<data->fp.numDim; idim++){
      rms = sqrt(pow(data->pd[ent->node].avgVel[idim],2.0)+pow(data->pd[ent->node].stdVel[idim],2.0));
      ujFluct[idim] = R*uiFluct[idim]+rms*sqrt(1-R*R)*epsilon;
      uj[idim] += ujFluct[idim];
    }
  }

  //***Calculate relative velocity***//

  for(idim=0; idim<data->fp.numDim; idim++){
    uijRel[idim] = ui[idim]-uj[idim];
  }

  //***Collision probability criterion***//

  conc = numDens/plasinterface_getNodVol(ent->node);
  FreqColl = PI/4.0*(di+dj)*(di+dj)*plas_CalcVectorLength(data->fp.numDim,uijRel)*conc;
  Pcoll = FreqColl*dtLagr;

  //***Perform collision***//

  if(plas_RandomDouble()<Pcoll){

    //***Transform into 2nd coordinate system***//

    x[0] = y[1] = z[2] = 1.0;
    x[1] = x[2] = y[0] = y[2] = z[0] = z[1] = 0.0;

    for(idim=0; idim<data->fp.numDim; idim++){
      xPr[idim] = uijRel[idim]/plas_CalcVectorLength(data->fp.numDim,uijRel);
    }

    plas_CalcCrossProduct_3D(yPr,xPr,x);
    plas_CalcCrossProduct_3D(zPr,xPr,yPr);

    pos1Pr[0] = plas_CalcVectScalarProduct(data->fp.numDim,ent->pos,xPr);
    pos1Pr[1] = plas_CalcVectScalarProduct(data->fp.numDim,ent->pos,yPr);
    pos1Pr[2] = plas_CalcVectScalarProduct(data->fp.numDim,ent->pos,zPr);

    while(L>1){
      Y = plas_RandomDouble();
      Z = plas_RandomDouble();
      L = sqrt(Y*Y+Z*Z);
    }

    Phi = asin(L);
    Psi = atan(Y/Z);

    plas_CalcRotationMatrix_3D(Phi,M,2);
    plas_CalcMatVectScalarProduct_3D(pos2Pr,M,xPr);
    plas_CalcVectorRotation_3D(Psi,M,xPr);
    plas_CalcMatVectScalarProduct_3D(pos2Pr,M,pos2Pr);

    for(idim=0; idim<data->fp.numDim; idim++){
      pos2Pr[idim] *= (dj+di)/2.0;
    }

    pos2[0] = plas_CalcVectScalarProduct(data->fp.numDim,pos2Pr,x);
    pos2[1] = plas_CalcVectScalarProduct(data->fp.numDim,pos2Pr,y);
    pos2[2] = plas_CalcVectScalarProduct(data->fp.numDim,pos2Pr,z);

    //***Transform into 3rd coordinate system***//

    for(idim=0; idim<data->fp.numDim; idim++){
      posRel[idim] = pos1[idim]-pos2[idim];
    }

    for(idim=0; idim<data->fp.numDim; idim++){
      xPrPr[idim] = (posRel[idim])/plas_CalcVectorLength(data->fp.numDim,posRel);
    }

    plas_CalcCrossProduct_3D(zPrPr,xPrPr,uijRel);

    for(idim=0; idim<data->fp.numDim; idim++){
      zPrPr[idim] /= plas_CalcVectorLength(data->fp.numDim,zPrPr);
    }

    plas_CalcCrossProduct_3D(yPrPr,zPrPr,xPrPr);

    uiPrPr[0] = plas_CalcVectScalarProduct(data->fp.numDim,ui,xPrPr);
    uiPrPr[1] = plas_CalcVectScalarProduct(data->fp.numDim,ui,yPrPr);
    uiPrPr[2] = plas_CalcVectScalarProduct(data->fp.numDim,ui,zPrPr);
    ujPrPr[0] = plas_CalcVectScalarProduct(data->fp.numDim,uj,xPrPr);
    ujPrPr[1] = plas_CalcVectScalarProduct(data->fp.numDim,uj,yPrPr);
    ujPrPr[2] = plas_CalcVectScalarProduct(data->fp.numDim,uj,zPrPr);
    uijRelPrPr[0] = plas_CalcVectScalarProduct(data->fp.numDim,uijRel,xPrPr);
    uijRelPrPr[1] = plas_CalcVectScalarProduct(data->fp.numDim,uijRel,yPrPr);
    uijRelPrPr[2] = plas_CalcVectScalarProduct(data->fp.numDim,uijRel,zPrPr);

    //***Sliding or non-sliding collision***//

    Mi = data->md.rhoDisp*PI*(dj/2.0)*(dj/2.0);
    Mj = data->md.rhoDisp*PI*(dj/2.0)*(dj/2.0);
    Jx = -(1.0-e)*uiPrPr[1]*(Mi*Mj)/(Mi+Mj);
    length = plas_CalcVectorLength(data->fp.numDim,uijRel);

    if(fabs(length)<7.0/2.0*mu_static*(1.0+e)*fabs(uiPrPr[0])){
      Jy = -2.0/7.0*uijRelPrPr[1]*Mi*Mj/(Mi+Mj);
      Jz = -2.0/7.0*uijRelPrPr[2]*Mi*Mj/(Mi+Mj);
    }else{
      Jy = -mu_dynamic*uijRelPrPr[1]/length*abs(Jx);
      Jz = -mu_dynamic*uijRelPrPr[2]/length*abs(Jx);
    }

    uiPrPrNew[0] = uiPrPr[0]+Jx/Mi;
    uiPrPrNew[1] = uiPrPr[1]+Jy/Mi;
    uiPrPrNew[2] = uiPrPr[2]+Jz/Mi;

    uiNew[0] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,x);
    uiNew[1] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,y);
    uiNew[2] = plas_CalcVectScalarProduct(data->fp.numDim,uiPrPrNew,z);

    data->sd.coll++;

    for(idim=0; idim<data->fp.numDim; idim++){
      ent->vel[idim] = uiNew[idim];
    }
  }

  /*
  for (int r=0; r<data->fp.numDim; ++r)
    delete[] M[r];
  delete[] M;
  */
  for (idim=0; idim<data->fp.numDim; ++idim)
    free(M[idim]);
  free(M);
}

