
#include <cmath>
#include "mfactory.h"
#include "plas_material.h"


namespace aux {


// define some materials
struct air : plas_material
{
  void update(double _T, double _p) { plas_material::update(_T,_p);
    rho     = (28.97e-3)/(Ru*T)*(p-exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T))+4.0*0.075/2e-5);
    satPres = exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T));
    He      = 2.1583e-7*exp(1342.0*(1.0/T-1.0/293.15));
  }
  air() {
    flowtype = FLOW_BUBBLY;
    update(273.15,101325.);

    rho             = 1.225;
    mu              = 1.7894e-5;
    cp              = 1.006;
    k               = 0.0242;
    sig             = 0.07275;  // as bubbly flow air-water
    eps             = 0.;
  //satPres         = dependent
    vapPres         = 2337.;
    latHeat         = 0.;
    molarMass       = 28.97;
    molarMassVap    = 0.;
    binaryDiffCoeff = 0.;
    massDiffCoeff   = 2.0e-9;
  //He              = dependent
  }
};


struct copper : plas_material
{
  copper() {
    flowtype = FLOW_PARTIC;

    rho             = 8920.;
  //mu              = ?
    cp              = 381.;
    k               = 387.6;
  //sig             = ?
    eps             = 0.03;
  //satPres         = ?
  //vapPres         = ?
  //latHeat         = ?
  //molarMass       = ?
  //molarMassVap    = ?
  //binaryDiffCoeff = ?
  //massDiffCoeff   = ?
  //He              = ?
  }
};


struct hydrogen : plas_material
{
  hydrogen() {
    flowtype = FLOW_BUBBLY;

    rho             = 0.0943;
    mu              = 8.411e-6;
    cp              = 14283.;
    k               = 0.1672;
    sig             = 0.07275;
  //eps             = ?
  //satPres         = ?
    vapPres         = 2337.;
  //latHeat         = ?
    molarMass       = 2.0159;
  //molarMassVap    = ?
  //binaryDiffCoeff = ?
    massDiffCoeff   = 1.61e-9;
    He              = 1.0/(7.099e4*101.325);
  }
};


struct nheptane : plas_material
{
  void update(double _T, double _p) { plas_material::update(_T,_p);
    rho     = -941.03+19.9618*T-0.08612051*pow(T,2.0)+1.579494e-4*pow(T,3.0)-1.089345e-7*pow(T,4.0);
    cp      = 799.3401062-126.5095282565*T+0.5279613848638*pow(T,2.0)-1664.890863e-8*pow(T,3.0)+644.6826474e-11*pow(T,4.0);
    k       = 0.25844890110-4.5549450549e-4*T;
    sig     = 0.059*pow(1.0-T/540.17,0.121);
    satPres = 1.0e5*pow(10.0,4.02677-1258.34/(T-53.85));
    latHeat = 317.8e3*pow((540.17-T)/(540.17-371.4),0.38);
    binaryDiffCoeff =
      ((1.8583e-7)/(p*pow(6.498,2.0)*(1.06036/(pow(T*399.3,0.15610)+0.19300/exp(0.47635*T*399.3))+1.76474/(exp(3.89411*T*399.3)))))
      *pow(pow(T,1.5)*(1.0/100.21+1.0/28.01),0.5);
  }
  nheptane() {
    flowtype = FLOW_DROPLET;
    update(273.15,101325.);

  //rho             = dependent
    mu              = 0.386e-3;
  //cp              = dependent
  //k               = dependent
  //sig             = dependent
    eps             = 1.;  // FIXME T and D dependent
  //satPres         = dependent
  //vapPres         = ?
  //latHeat         = dependent
    molarMass       = 100.21;
    molarMassVap    = 28.01;
  //binaryDiffCoeff = dependent
  //massDiffCoeff   = ?
  //He              = ?
  }
};


struct nitrogen : plas_material
{
  void update(double _T, double _p) { plas_material::update(_T,_p);
    mu  = (6.5592e-7*pow(T,0.6081))/(1.0+54.715/T);
    cp  = (6.50+0.001*T)*4.184/(28.01e-3);
    k   = 2.5*(6.5592e-7*pow(T,0.6081))/(1.0+54.715/T)*((6.50+0.001*T)*4.184/(28.01e-3)-8.314)/28.01;
  }
  nitrogen() {
    flowtype = FLOW_BUBBLY;
    update(273.15,101325.);

    rho             = 1.25;
  //mu              = dependent
  //cp              = dependent
  //k               = dependent
  //sig             = ?
  //eps             = ?
  //satPres         = ?
  //vapPres         = ?
  //latHeat         = ?
  //molarMass       = ?
  //molarMassVap    = ?
  //binaryDiffCoeff = ?
  //massDiffCoeff   = ?
  //He              = ?
  }
};


struct oxygen : plas_material
{
  oxygen() {
    flowtype = FLOW_BUBBLY;

    rho             = 1.2999;
    mu              = 1.919e-5;
    cp              = 919.31;
    k               = 0.0246;
  //sig             = ?
  //eps             = ?
  //satPres         = ?
  //vapPres         = ?
  //latHeat         = ?
    molarMass       = 15.9994;
  //molarMassVap    = ?
  //binaryDiffCoeff = ?
  //massDiffCoeff   = ?
  //He              = ?
  }
};


struct poly : plas_material
{
  poly() {
    flowtype = FLOW_PARTIC;

    rho             = 1050.;
  //mu              = ?
    cp              = 1300.;
    k               = 0.08;
  //sig             = ?
  //eps             = ?
  //satPres         = ?
  //vapPres         = ?
  //latHeat         = ?
  //molarMass       = ?
  //molarMassVap    = ?
  //binaryDiffCoeff = ?
  //massDiffCoeff   = ?
  //He              = ?
  }
};


struct water : plas_material
{
  void update(double _T, double _p) { plas_material::update(_T,_p);
    mu = 2.414e-5*pow(10.0,(247.8/(T-140.0)));
    cp = (1.0/18.015)*(276376.0-2090.1*T+8.125*pow(T,2.0)-0.014116*pow(T,3.0)+9.3701e-6*pow(T,4.0));
    k  = (2.9041*(3.0+pow(1.0-(T/611.966),2.0/3.0)))/(3.0+20.0*pow(1.0-(273.0/611.966),2.0/3.0));
    satPres = exp((-6096.9385/T)+21.2409642-2.711193e-2*T+1.673952e-5*pow(T,2.0)+2.433502*log(T));
    latHeat = (5.2053e7/18.015)*pow((1.0-(T/611.966)),0.3199-0.212*(T/611.966)+0.25795*pow((T/611.966),2.0));;
  }
  water() {
    flowtype = FLOW_DROPLET;
    update(273.15,101325.);

    rho             = 998.2;
  //mu              = dependent
  //cp              = dependent
  //k               = dependent
    sig             = 0.07275;
    eps             = 0.95;
  //satPres         = dependent
  //vapPres         = ?
  //latHeat         = dependent
    molarMass       = 18.015;
    molarMassVap    = 29.;
    binaryDiffCoeff = 0.220e-4;
  //massDiffCoeff   = ?
  //He              = ?
  }
};


void plas_material_database()
{
  m::Register< plas_material, air      > mat01(1,          "air",      "");
  m::Register< plas_material, copper   > mat02(2,"Cu",  "","copper",   "");
  m::Register< plas_material, hydrogen > mat03(2,"H2",  "","hydrogen", "");
  m::Register< plas_material, nheptane > mat04(1,          "n-Heptane","");
  m::Register< plas_material, nitrogen > mat05(2,"N2",  "","nitrogen", "");
  m::Register< plas_material, oxygen   > mat06(2,"O2",  "","oxygen",   "");
  m::Register< plas_material, poly     > mat07(1,"C8H8",""               );
  m::Register< plas_material, water    > mat08(2,"H2O", "","water",    "");
}


}  // namespace aux

