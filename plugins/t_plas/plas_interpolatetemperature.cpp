
#include "common.h"


/*
 * This function performs temperature interpolation from the
 * fluid flow solver.
 */

void plas_InterpolateTemperature(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow, double step)
{
  int idim;
  double impFac[ent->edata.numElmNodes];

  //***Compute impact factors of element nodes***//

  plas_CalcNodeImpactFactors(data,ent,impFac);

  //***Flow temperature***//

  flow->temp = 0.0;
  for(idim=0; idim<ent->edata.numElmNodes; idim++){
    flow->temp += impFac[idim]*((1.0-step)*plasinterface_getTemperatureOld(ent->edata.elmNodes[idim])
                                    + step*plasinterface_getTemperature(ent->edata.elmNodes[idim]));
  }
}

