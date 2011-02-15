
#include "common.h"


/*
 * This file includes functions to compute flow coefficients.
 *
 * This function computes all flow coefficients.
 */

void plas_CalcEntityCoefficients(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent, LOCAL_FLOW_VARIABLES *flow)
{
  //***Compute flow coefficients***//

  ent->reynolds = plas_CalcDispReynolds(data->fp.nuCont,ent->diam,ent->normVel);
  ent->dragCoeff = plas_CalcDragCoeff(data->rp.flowType,ent->reynolds);
  ent->liftCoeff = plas_CalcLiftCoeff(data->rp.flowType);
  ent->kinRespTime = plas_CalcKinematicResponseTime(data,ent);
  ent->spalding = plas_CalcSpaldingNumber(data,flow->pressure);
  ent->prandtl = plas_CalcPrandtlNumber(data);
  ent->nusselt = plas_CalcNusseltNumber(data->ip.evapModel,ent->reynolds,ent->spalding,ent->prandtl);
  ent->thermRespTime = plas_CalcThermalResponseTime(data,ent->diam);
  ent->schmidt = plas_CalcSchmidtNumber(data);
  ent->sherwood = plas_CalcSherwoodNumber(data->ip.evapModel,ent->reynolds,ent->schmidt,ent->spalding);
  ent->massTrCoeff = plas_CalcMassTransferCoeff(data,ent->sherwood,ent->spalding);
  ent->pressBubble = plas_CalcPressBubble(data,ent->diam,flow->pressure);
  ent->rhoBubble = plas_CalcRhoBubble(data,ent->temp, ent->pressBubble);
  ent->concInterf = plas_CalcConcInterf(data, ent->pressBubble);
}

