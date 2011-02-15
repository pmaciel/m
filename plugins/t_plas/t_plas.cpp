
#include "plas_driver.h"


PLAS_DATA          dplasdata;
DRIVER_PARAMETERS  dparam;
DRIVER_GAMBIT_MESH dmesh;
DRIVER_FLOW_FIELD  dflow;


/*
 * Main file of the PLaS driver program, which is a pseudo
 * flow solver. It gives a steady state flow field without
 * solving the flow. It serves for testing purposes of PLaS.
 *
 * Main routine of the PLaS driver.
 */

int main(int argc, char *argv[])
{
  int i;

#ifdef MPI
  int proc,rank;
#endif

  printf("-------------------\n");
  printf("PLaS driver program\n");
  printf("-------------------\n");

  //***Read driver parameters***//

  printf("Reading driver.conf...\n");
  plasdriver_ReadDriverDataFile(&dparam);

  //***Read mesh***//

  printf("Reading Gambit neutral file \"%s.neu\"...\n",dparam.gridString);
  plasdriver_ReadGambitNeutralFile(&dmesh,&dparam);

  //***Geometry calculations***//

  printf("Performing geometry calculations...\n");
  printf("   Elements around nodes...\n");
  plasdriver_CalcElmsAroundNode(&dmesh);
  printf("   Element neighbours...\n");
  plasdriver_CalcElementNeighbours(&dmesh);
  printf("   Element normals...\n");
  plasdriver_CalcElementNormals(&dmesh);
  printf("   Element volumes...\n");
  plasdriver_CalcElementVolumes(&dmesh);
  printf("   Nodal volumes...\n");
  plasdriver_CalcNodalVolumes(&dmesh);

  dmesh.domainVolume = 0.0;
  for(i=0; i<dmesh.numElm; i++){
    dmesh.domainVolume += dmesh.elmVolumes[i];
    if(i==0){
      dmesh.minElmVolume = dmesh.elmVolumes[i];
      dmesh.maxElmVolume = dmesh.elmVolumes[i];
    } else{
      if(dmesh.elmVolumes[i]<dmesh.minElmVolume){
        dmesh.minElmVolume = dmesh.elmVolumes[i];
      }
      if(dmesh.elmVolumes[i]>dmesh.maxElmVolume){
        dmesh.maxElmVolume = dmesh.elmVolumes[i];
      }
    }
  }

  dparam.numUnk = dmesh.numDim+2;

  //***Define or read steady-state flow field***//

  printf("Initializing flow field...\n");
  dflow.p = (double*)calloc(dmesh.numNod,sizeof(double));
  dflow.T = (double*)calloc(dmesh.numNod,sizeof(double));
  dflow.u = (double**)calloc(dmesh.numNod,sizeof(double));
  for(i=0; i<dmesh.numNod; i++){
    dflow.u[i] = (double*)calloc(dmesh.numDim,sizeof(double));
  }
  plasdriver_InitFlowField(&dmesh,dparam.material,&dflow);

  //***Write Tecplot file of flow field***//

  printf("Writing out flow field file \"%s.plt\"\n",dparam.gridString);
  plasdriver_WriteTecplot(&dmesh,&dparam,&dflow);

  //***Initialize PLaS***//

  printf("--------------------\n");
  printf("Initializing PLaS...\n");
#ifdef MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&proc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  initPLaS(&dplasdata);
  printf("--------------------\n");

  //***Perform PLaS iterations***//

  for(dparam.iter=1; dparam.iter<=dparam.numIter; dparam.iter++){
    printf("Iteration %d\n",dparam.iter);
    runPLaS(&dplasdata);
  }

  //***Terminate PLaS***//

  printf("--------------------\n");
  printf("Terminating PLaS...\n");
  terminatePLaS(&dplasdata);
  plasdriver_FreeGambitMemory(&dmesh);
  for(i=0; i<dmesh.numNod; i++){
    free(dflow.u[i]);
  }
  free(dflow.p);
  free(dflow.u);
  free(dflow.T);

#ifdef MPI
  MPI_Finalize();
#endif
  printf("-------------------\n");

  return 0;
}

