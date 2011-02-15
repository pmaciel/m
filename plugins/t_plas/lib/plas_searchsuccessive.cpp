
#include "common.h"


/*
 * This file includes all functinality to perform an element
 * search for a dispersed entity.
 *
 * This function performs a successive neigbour search routine
 * following the algorithm of Loehner et al.
 */

void plas_SearchSuccessive(PLAS_DATA *data, LOCAL_ENTITY_VARIABLES *ent)
{
  int elmFound = 0;
  int leftDomain = 0;
  int lastElm = ent->elm;
  int lastNode = ent->node;
  int neighbourElm,idx;
  double dist;

  //***Successive neighbour search***//

  do{

    //***Set geometry of current element***//

    plas_SetElementGeometry(data->fp.numDim,ent);

    //***Find the face with the minimum distance to the entity***//

    plas_FindMinimumElementFaceDistance(data->fp.numDim,ent,&idx,&dist);

    if(dist>(0.0-data->rp.errTol)){

      //***In case the minimum distance is positive, the element is found***//

      elmFound = 1;
    } else{

      //***Search the neighbour element in directin of the minimum face distance***//

      neighbourElm = plasinterface_getElmNeighbour(ent->elm,idx);

      if(neighbourElm==-1){
        leftDomain = 1;
      } else{
        ent->elm = neighbourElm;
        if(lastElm!=-1){
          lastElm = neighbourElm;
        }
      }
    }
  } while(!elmFound && !leftDomain);

  //***Disable entity in case of not found***//

  if(elmFound){
    ent->flag = DFLAG_ENABLED;
    ent->node = plas_FindNearestElementNode(data,ent);
  } else{
    ent->flag = DFLAG_LEFT;
    ent->elm = lastElm;
    ent->node = lastNode;
  }
}

