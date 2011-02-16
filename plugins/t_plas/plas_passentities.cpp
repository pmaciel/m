
#include "common.h"


/*
 * This file includes a function to pass an entity from one
 * process of a multi-processor job to another one after the
 * entity has crossed a process boundary.
 *
 * Function to pass entities from one process to another.
 */

void plas_PassEntities(PLAS_DATA *data)
{
  LOCAL_ENTITY_VARIABLES ent;
  int numProc        = plas_MpiGetNumProc();
  int irank          = plas_MpiGetRank();
  int *disps         = new int[numProc];
  int *recs          = new int[numProc];
  int leftCtrOneProc = data->sd.leftproc;
  int leftCtrAllProc = plas_MpiAllSumInt(data->sd.leftproc);
  int *foundProc     = new int[leftCtrAllProc];
  int *foundElm      = new int[leftCtrAllProc];
  int *foundNode     = new int[leftCtrAllProc];
  int idim,jdim,ient,jent,ielm,inod;
  char errMessage[100];

  double *leftDataOneProc = new double[(2*data->fp.numDim+1)*leftCtrOneProc];
#ifdef MPI
  double *leftDataAllProc = new double[(2*data->fp.numDim+1)*leftCtrAllProc];
#else
  double *leftDataAllProc = leftDataOneProc;
#endif

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(data->fp.numDim,&ent);

  //***Collect all entities that left a process***//

  jdim = 0;
  for(ient=0; ient<data->ip.numMaxEnt; ient++){
    if(data->ed[ient].flag==DFLAG_PASS){
      data->ed[ient].flag = DFLAG_DISABLED;
      for(idim=0; idim<data->fp.numDim; idim++){
        leftDataOneProc[jdim] = data->ed[ient].position[idim];
        jdim++;
      }
      for(idim=0; idim<data->fp.numDim; idim++){
        leftDataOneProc[jdim] = data->ed[ient].velocity[idim];
        jdim++;
      }
      leftDataOneProc[jdim] = data->ed[ient].diameter;
      jdim++;
    }
  }

  //***Gather all entities that left any process***//

#ifdef MPI
  MPI_Allgather(&leftCtrOneProc,1,MPI_INT,recs,1,MPI_INT,MPI_COMM_WORLD);

  disps[0] = 0;
  for (int iproc=1; iproc<numProc; iproc++){
    disps[iproc] = disps[iproc-1] + recs[iproc-1];
  }

  for (int iproc=0; iproc<numProc; iproc++){
    recs[iproc] *= (2*data->fp.numDim+1);
    disps[iproc] *= (2*data->fp.numDim+1);
  }

  MPI_Allgatherv(leftDataOneProc,(2*data->fp.numDim+1)*leftCtrOneProc,MPI_DOUBLE,leftDataAllProc,recs,disps,MPI_DOUBLE,MPI_COMM_WORLD);
#endif

  //***Search for all free entities on all processes***//

  idim = 0;
  for(ient=0; ient<leftCtrAllProc; ient++){

    for(jdim=0; jdim<data->fp.numDim; jdim++){
      ent.pos[jdim] = leftDataAllProc[idim+jdim];
    }

    //***Element search***//

    plas_SearchDomainParallel(data,&ent);

    if(ent.flag==DFLAG_ENABLED){
      foundProc[ient] = irank;
      foundElm[ient] = ent.elm;
      foundNode[ient] = ent.node;
    } else{
      foundProc[ient] = -1;
    }
    idim += 2*data->fp.numDim+1;
  }

  //***Determine a definite process on which an entity is found***//

  plas_MpiAllMaxIntArray(foundProc,leftCtrAllProc);

  //***Broadcast information of found entities***//

  idim = 0;
  jent = 0;

  for(ient=0; ient<leftCtrAllProc; ient++){

    if(foundProc[ient]>-1){
      ielm = foundElm[ient];
      plas_MpiBroadcastInt(&ielm,1,foundProc[ient]);
      inod = foundNode[ient];
      plas_MpiBroadcastInt(&inod,1,foundProc[ient]);
    }

    if(foundProc[ient]!=irank){
      idim += 2*data->fp.numDim+1;
      continue;
    }

    while(data->ed[jent].flag==DFLAG_ENABLED){
      jent++;
      if(jent==data->ip.numMaxEnt){
        sprintf(errMessage,"Memory exceeded when passing trajectory information.");
        plas_TerminateOnError(errMessage);
      }
    }

    //***Generate entity on the new process***//

    data->sd.passed++;
    data->ed[jent].flag = DFLAG_ENABLED;
    data->ed[jent].element = ielm;
    data->ed[jent].node = inod;

    for(jdim=0; jdim<data->fp.numDim; jdim++){
      data->ed[jent].position[jdim] = leftDataAllProc[idim];
      idim++;
    }

    for(jdim=0; jdim<data->fp.numDim; jdim++){
      data->ed[jent].velocity[jdim] = leftDataAllProc[idim];
      idim++;
    }

    data->ed[jent].diameter = leftDataAllProc[idim];
    idim++;
  }

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);

  delete[] foundProc;
  delete[] foundElm;
  delete[] foundNode;
  delete[] leftDataOneProc;
#ifdef MPI
  delete[] leftDataAllProc;
#endif
  delete[] disps;
  delete[] recs;
}

