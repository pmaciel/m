
#include "common.h"


/*
 * This function contains the material database for dispersed
 * particles, droplets and bubbles.
 */

void plas_CalcMaterialData(PLAS_DATA *data, double T, double p)
{
  if(data->ip.material==MAT_COPPER){

    //***Material data for copper***//

    data->md.rhoDisp = 8920.0;
    //data->md.muDisp = 0.0;
    data->md.cpDisp = 381.0;
    data->md.kDisp = 387.6;
    //data->md.sigDisp = 0.0;
    data->md.epsDisp = 0.03;
    //data->md.satPresDisp = 0.0;
    //data->md.vapPres = 0.0;
    //data->md.latHeatDisp = 0.0;
    //data->md.molarMassDisp = 0.0;
    //data->md.molarMassDispVap = 0.0;
    //data->md.binaryDiffCoeff = 0.0;
    //data->md.massDiffCoeff = 0.0;
    //data->md.HeDisp = 0.0;

  } else if(data->ip.material==MAT_POLY){

    //***Material data for C8H8***//

    data->md.rhoDisp = 1050.0;
    //data->md.muDisp = 0.0;
    data->md.cpDisp = 1300.0;
    data->md.kDisp = 0.08;
    //data->md.sigDisp = 0.0;
    //data->md.epsDisp = 0.0;
    //data->md.satPresDisp = 0.0;
    //data->md.vapPres = 0.0;
    //data->md.latHeatDisp = 0.0;
    //data->md.molarMassDisp = 0.0;
    //data->md.molarMassDispVap = 0.0;
    //data->md.binaryDiffCoeff = 0.0;
    //data->md.massDiffCoeff = 0.0;
    //data->md.HeDisp = 0.0;

  } else if(data->ip.material==MAT_WATER){

    //***Material data for water***//

    data->md.rhoDisp = 998.2;
    data->md.muDisp = 2.414e-5*pow(10.0,(247.8/(T-140.0)));
    data->md.cpDisp = (1.0/18.015)*(276376.0-2090.1*T+8.125*pow(T,2.0)-0.014116*pow(T,3.0)+9.3701e-6*pow(T,4.0));
    data->md.kDisp = (2.9041*(3.0+pow(1.0-(T/611.966),2.0/3.0)))/(3.0+20.0*pow(1.0-(273.0/611.966),2.0/3.0));
    data->md.sigDisp = 0.07275;
    data->md.epsDisp = 0.95;
    data->md.satPresDisp = exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T));
    //data->md.vapPres = 0.0;
    data->md.latHeatDisp = (5.2053e7/18.015)*pow((1.0-(T/611.966)),0.3199-0.212*(T/611.966)+0.25795*pow((T/611.966),2.0));;
    data->md.molarMassDisp = 18.015;
    data->md.molarMassDispVap = 29.0;
    data->md.binaryDiffCoeff = 0.220e-4;
    //data->md.massDiffCoeff = 0.0;
    //data->md.HeDisp = 0.0;

  } else if(data->ip.material==MAT_NHEPTANE){

    //***Material data for n-heptane***//

    data->md.rhoDisp = -941.03+19.9618*T-0.08612051*pow(T,2.0)+1.579494e-4*pow(T,3.0)-1.089345e-7*pow(T,4.0);
    data->md.muDisp = 0.386e-3;
    data->md.cpDisp = 799.3401062-126.5095282565*T+0.5279613848638*pow(T,2.0)-1664.890863e-8*pow(T,3.0)+644.6826474e-11*pow(T,4.0);
    data->md.kDisp = 0.25844890110-4.5549450549e-4*T;
    data->md.sigDisp = 0.059*pow(1.0-T/540.17,0.121);
    data->md.epsDisp = 1.0; // T and D dependancy still to be implemented...
    data->md.satPresDisp = 1.0e5*pow(10.0,4.02677-1258.34/(T-53.85));
    //data->md.vapPres = 0.0;
    data->md.latHeatDisp = 317.8e3*pow((540.17-T)/(540.17-371.4),0.38);
    data->md.molarMassDisp = 100.21;
    data->md.molarMassDispVap = 28.01;
    data->md.binaryDiffCoeff =
      ((1.8583e-7)/(p*pow(6.498,2.0)*(1.06036/(pow(T*399.3,0.15610)+0.19300/exp(0.47635*T*399.3))+1.76474/(exp(3.89411*T*399.3)))))
      *pow(pow(T,1.5)*(1.0/100.21+1.0/28.01),0.5);
    //data->md.massDiffCoeff = 0.0;
    //data->md.HeDisp = 0.0;

  } else if(data->ip.material==MAT_HYDROGEN){

    //***Material data for hydrogen***//

    data->md.rhoDisp = 0.0943;
    data->md.muDisp = 8.411e-6;
    data->md.cpDisp = 14283.0;
    data->md.kDisp = 0.1672;
    data->md.sigDisp = 0.07275;
    //data->md.epsDisp = 0.0;
    //data->md.satPresDisp = 0.0;
    data->md.vapPres = 2337.0;
    //data->md.latHeatDisp = 0.0;
    data->md.molarMassDisp = 2.0159;
    //data->md.molarMassDispVap = 0.0;
    //data->md.binaryDiffCoeff = 0.0;
    data->md.massDiffCoeff = 1.61e-9;
    data->md.HeDisp = 1.0/(7.099e4*101.325);

  } else if(data->ip.material==MAT_OXYGEN){

    //***Material data for oxygen***//

    data->md.rhoDisp = 1.2999;
    data->md.muDisp = 1.919e-5;
    data->md.cpDisp = 919.31;
    data->md.kDisp = 0.0246;
    //data->md.sigDisp = 0.0;
    //data->md.epsDisp = 0.0;
    //data->md.satPresDisp = 0.0;
    //data->md.vapPres = 0.0;
    //data->md.latHeatDisp = 0.0;
    data->md.molarMassDisp = 15.9994;
    //data->md.molarMassDispVap = 0.0;
    //data->md.binaryDiffCoeff = 0.0;
    //data->md.massDiffCoeff = 0.0;
    //data->md.HeDisp = 0.0;


  } else if(data->ip.material==MAT_AIR){

    //***Material data for air***//

    data->md.rhoDisp = (28.97e-3)/(Ru*T)*(p-exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T))+4.0*0.075/2e-5);//1.225;
    data->md.muDisp = 1.7894e-5;
    data->md.cpDisp = 1006.0;
    data->md.kDisp = 0.0242;
    data->md.sigDisp = 0.07275; //As bubbly flow air-water
    data->md.epsDisp = 0.0;
    data->md.satPresDisp = exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T));
    data->md.vapPres = 2337.0;
    data->md.latHeatDisp = 0.0;
    data->md.molarMassDisp = 28.97;
    data->md.molarMassDispVap = 0.0;
    data->md.binaryDiffCoeff = 0.0;
    data->md.massDiffCoeff = 2.0e-9;
    data->md.HeDisp = 2.1583e-7*exp(1342.0*(1.0/T-1.0/293.15));
  }
}

