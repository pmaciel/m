#ifndef PLAS_PLAS_MATERIAL_H
#define PLAS_PLAS_MATERIAL_H


#ifndef Ru
#define Ru 8314.472
#endif


/**
 * type of flow associated to the dispersed phase material
 * (particle, droplet, bubble)
 */
enum flowtype_t { FLOW_PARTIC=1, FLOW_DROPLET, FLOW_BUBBLY };


/**
 * structure to hold material physical properties
 */
struct plas_material
{
  // methods
  plas_material() : T(273.15/*K*/), p(101325./*Pa*/) {}
  virtual void update(double _T, double _p) { T=_T; p=_p; }

  // independent physical properties
  double
    T,
    p;
  flowtype_t flowtype;  // type of flow (particle, droplet, bubble)

  // generic (dependent) physical properties
  double
    rho,  // density
    mu,   // dynamic viscosity
    cp,   // specific heat capacity
    k;    // thermal conductivity

  // dispersed phase-specific properties
  double
    sig,              // surface tension
    eps,              // emissivity
    satPres,          // saturation pressure evaluated at surface of entities
    vapPres,          // vapour pressure
    latHeat,          // specific latent heat
    molarMass,        // molar mass
    molarMassVap,     // molar mass at vapour phase
    binaryDiffCoeff,  // binary diffusion coefficient
    massDiffCoeff,    // mass diffusivity for a gas in a liquid
    He;               // Henry law constant
};


namespace aux {


/**
 * initialize materials database
 */
void plas_material_database();


}  // namespace aux


#endif  // PLAS_PLAS_MATERIAL_H

