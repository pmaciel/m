
#include "common.h"


/*
 * This file includes all functinality to perform an element
 * search for a dispersed entity.
 */

void plas_SearchDomainParallel(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent)
{
  int bfCtr = 0;
  int elmFoundAnyProc = 0;
  int idx,elmFound,procFound;
  double dist;

  //***Perform a successive neighbour search (find by chance)***//

  do{
    ent->elm = plas_RandomInteger(0,data->fp.numElm-1);
    bfCtr++;
    plas_SearchSuccessive(data,ent);
    if(ent->flag==DFLAG_ENABLED){
      elmFound = 1;
      procFound = plas_MpiGetRank();
    }else{
      elmFound = 0;
      procFound = plas_MpiGetNumProc();
    }

  } while(!elmFound && bfCtr<5);

  //***Determins if an element has been found on any process***//

  elmFoundAnyProc = plas_MpiAllMaxInt(elmFound);

  //***Brute force search (relevant elements given by flow solver)***//

  if(!elmFoundAnyProc){
    ent->elm = plasinterface_StartElementSearch(ent->pos);
    do{
      plas_SetElementGeometry(data->fp.numDim,ent);
      plas_FindMinimumElementFaceDistance(data->fp.numDim,ent,&idx,&dist);
      if(dist>(0.0-data->rp.errTol)){
        elmFound = 1;
        procFound = plas_MpiGetRank();
      }else{
        elmFound = 0;
        procFound = plas_MpiGetNumProc();
        ent->elm++;
      }
    } while(!elmFound && ent->elm<=plasinterface_EndElementSearch(ent->pos));
  }

  //***In case an element is found on more than opne process, take the minimum one***//

  procFound = plas_MpiAllMinInt(procFound);

  //***Disable entity in case of not found***//

  if(elmFound && plas_MpiGetRank()==procFound){
    ent->flag = DFLAG_ENABLED;
    ent->node = plas_FindNearestElementNode(data,ent);
  } else{
    ent->flag = DFLAG_DISABLED;
    ent->elm = -1;
    ent->node = -1;
  }
}

