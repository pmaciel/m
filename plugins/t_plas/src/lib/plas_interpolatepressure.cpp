
#include "common.h"


/*
 * This function performs pressure interpolation from the
 * fluid flow solver.
 */

void plas_InterpolatePressure(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim;
  double impFac[ent->edata.numElmNodes];

  //***Compute impact factors of element nodes***//

  plas_CalcNodeImpactFactors(data,ent,impFac);

  //***Flow pressure***//

  flow->pressure = 0.0;
  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    flow->pressure += impFac[idim]*((1.0-step)*plasinterface_getPressureOld(ent->edata.elmNodes[idim])
                                        + step*plasinterface_getPressure(ent->edata.elmNodes[idim]));
  }
}

