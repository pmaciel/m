
#include "common.h"


/*
 * This file includes all functionality concerning entities
 * interacting with boundaries of the computational domain,
 * such as wall bounces and entities leaving through outlets.
 *
 * This function finds the index of the boundary face through
 * which an entity leaft the computational domain.
 */

void plas_FindExitFace(int numBnd, int numDim, LOCAL_ENTITY_VARIABLES *ent, int *f, int *i, int *j)
{
  int ibnd,ifac=0,found;
  double d;

  //***Loop over boundaries and faces to find the exit face***//

  found = 0;
  for(ibnd=0; ibnd<numBnd; ibnd++){
    for(ifac=0; ifac<plasinterface_getNumBndFaces(ibnd); ifac++){
      if(ent->elm==plasinterface_getBndDomElm(ibnd,ifac,ent->elm)){
        d = plas_CalcWallFaceDistance(numDim,ent->pos,ibnd,ifac);
        if(d<0){found = 1;}
      }
      if(found){break;}
    }
    if(found){break;}
  }

  //***Set the return parameters***//

  *f = found;
  *i = ibnd;
  *j = ifac;
}

