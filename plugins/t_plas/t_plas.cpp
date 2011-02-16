
#ifdef MPI
#include <mpi.h>
#endif

#include "mfactory.h"
#include "t_plas.h"

using namespace m;


Register< mtransform,t_plas > mt_plas("-tplas","Particle Lagrangian Solver (PLaS)");


t_plas *pt_plas;


t_plas::t_plas() : m::mtransform() { ::pt_plas  = this; }
t_plas::~t_plas()                  { ::pt_plas  = NULL; }


void t_plas::transform(GetPot& o, mmesh& m)
{
  int i;

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
  dflow.u = (double**)calloc(dmesh.numNod,sizeof(double*));
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
  int proc, rank;
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
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the area of a triangle.
double t_plas::plasdriver_CalcAreaTriangle(double c[3][2])
{
  return c[0][0]*(c[1][1]-c[2][1]) + c[1][0]*(c[2][1]-c[0][1]) + c[2][0]*(c[0][1]-c[1][1]);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element neighbours.
void t_plas::plasdriver_CalcElementNeighbours(DRIVER_GAMBIT_MESH *dmesh)
{
  int neighbourFound,jnod,knod,lnod,ielm,jelm,kelm,lelm,ifac,faceNodes[4];

  dmesh->elmNeighbs = (int**)calloc(dmesh->numElm,sizeof(int*));
  for(ielm=0; ielm<dmesh->numElm; ielm++){
    dmesh->elmNeighbs[ielm] = (int*)calloc(dmesh->numElmFaces[ielm],sizeof(int));
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(dmesh,ielm,ifac,faceNodes);
      neighbourFound = 0;
      for(jnod=0; jnod<dmesh->numNodElms[faceNodes[0]]; jnod++){
        jelm = dmesh->nodElms[faceNodes[0]][jnod];
        for(knod=0; knod<dmesh->numNodElms[faceNodes[1]]; knod++){
          kelm = dmesh->nodElms[faceNodes[1]][knod];
          if(dmesh->numDim==2){
            if(jelm==kelm && jelm!=ielm){neighbourFound = 1;}
          }else if(dmesh->numDim==3){
            for(lnod=0; lnod<dmesh->numNodElms[faceNodes[2]]; lnod++){
              lelm = dmesh->nodElms[faceNodes[2]][lnod];
              if(jelm==kelm && jelm==lelm && jelm!=ielm){neighbourFound = 1;}
              if(neighbourFound==1){break;}
            }
          }
          if(neighbourFound==1){break;}
        }
        if(neighbourFound==1){break;}
      }
      if(neighbourFound==1){
        dmesh->elmNeighbs[ielm][ifac] = jelm;
      } else{
        dmesh->elmNeighbs[ielm][ifac] = -1;
      }
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element normals.
void t_plas::plasdriver_CalcElementNormals(DRIVER_GAMBIT_MESH *dmesh)
{
  int ielm,ifac,fnodes[4];

  for(ielm=0; ielm<dmesh->numElm; ielm++){
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(dmesh,ielm,ifac,fnodes);

      if(dmesh->elmTypes[ielm]==ELM_SIMPLEX && dmesh->numDim==2){

        dmesh->elmNorms[ielm][ifac][0] = dmesh->coords[fnodes[0]][1] - dmesh->coords[fnodes[1]][1];
        dmesh->elmNorms[ielm][ifac][1] = dmesh->coords[fnodes[1]][0] - dmesh->coords[fnodes[0]][0];

      } else if((dmesh->elmTypes[ielm]==ELM_SIMPLEX && dmesh->numDim==3)
                || (dmesh->elmTypes[ielm]==ELM_PRISM && ifac>2)
                || (dmesh->elmTypes[ielm]==ELM_PYRAMID && ifac>0)){

        dmesh->elmNorms[ielm][ifac][0] =
          0.5*((dmesh->coords[fnodes[2]][1]-dmesh->coords[fnodes[0]][1])
              *(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
              -(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
              *(dmesh->coords[fnodes[2]][2]-dmesh->coords[fnodes[0]][2]));
        dmesh->elmNorms[ielm][ifac][1] =
          0.5*((dmesh->coords[fnodes[2]][2]-dmesh->coords[fnodes[0]][2])
              *(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
              -(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
              *(dmesh->coords[fnodes[2]][0]-dmesh->coords[fnodes[0]][0]));
        dmesh->elmNorms[ielm][ifac][2] =
          0.5*((dmesh->coords[fnodes[2]][0]-dmesh->coords[fnodes[0]][0])
              *(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
              -(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
              *(dmesh->coords[fnodes[2]][1]-dmesh->coords[fnodes[0]][1]));

      } else if(dmesh->elmTypes[ielm]==ELM_QUAD){

        dmesh->elmNorms[ielm][ifac][0] = 2.0*(dmesh->coords[fnodes[0]][1] - dmesh->coords[fnodes[1]][1]);
        dmesh->elmNorms[ielm][ifac][1] = 2.0*(dmesh->coords[fnodes[1]][0] - dmesh->coords[fnodes[0]][0]);

      } else if(dmesh->elmTypes[ielm]==ELM_HEX
                || (dmesh->elmTypes[ielm]==ELM_PRISM && ifac<=2)
                || (dmesh->elmTypes[ielm]==ELM_PYRAMID && ifac==0)){

        dmesh->elmNorms[ielm][ifac][0] =
            ((dmesh->coords[fnodes[3]][1]-dmesh->coords[fnodes[0]][1])
            *(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
            -(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
            *(dmesh->coords[fnodes[3]][2]-dmesh->coords[fnodes[0]][2]));
        dmesh->elmNorms[ielm][ifac][1] =
            ((dmesh->coords[fnodes[3]][2]-dmesh->coords[fnodes[0]][2])
            *(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
            -(dmesh->coords[fnodes[1]][2]-dmesh->coords[fnodes[0]][2])
            *(dmesh->coords[fnodes[3]][0]-dmesh->coords[fnodes[0]][0]));
        dmesh->elmNorms[ielm][ifac][2] =
            ((dmesh->coords[fnodes[3]][0]-dmesh->coords[fnodes[0]][0])
            *(dmesh->coords[fnodes[1]][1]-dmesh->coords[fnodes[0]][1])
            -(dmesh->coords[fnodes[1]][0]-dmesh->coords[fnodes[0]][0])
            *(dmesh->coords[fnodes[3]][1]-dmesh->coords[fnodes[0]][1]));
      }
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element volumes.
void t_plas::plasdriver_CalcElementVolumes(DRIVER_GAMBIT_MESH *dmesh)
{
  int ielm;
  double c2[3][2],c3[4][3];

  dmesh->elmVolumes = (double*)calloc(dmesh->numElm,sizeof(double));

  for(ielm=0; ielm<dmesh->numElm; ielm++){

    if(dmesh->elmTypes[ielm]==ELM_SIMPLEX){
      if(dmesh->numDim==2){
        c2[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
        c2[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
        c2[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
        c2[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
        c2[2][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
        c2[2][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
        dmesh->elmVolumes[ielm] = plasdriver_CalcAreaTriangle(c2);
      } else if(dmesh->numDim==3){
        c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
        c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
        c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
        c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
        c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
        c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][1]][2];
        c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
        c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
        c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][2]][2];
        c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
        c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
        c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
        dmesh->elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      }
    } else if(dmesh->elmTypes[ielm]==ELM_QUAD){
      c2[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c2[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c2[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
      c2[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
      c2[2][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
      c2[2][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
      dmesh->elmVolumes[ielm] = plasdriver_CalcAreaTriangle(c2);
      c2[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c2[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c2[1][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
      c2[1][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
      c2[2][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c2[2][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      dmesh->elmVolumes[ielm] += plasdriver_CalcAreaTriangle(c2);
    } else if(dmesh->elmTypes[ielm]==ELM_HEX){
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][1]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      dmesh->elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][6]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][6]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][6]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][2]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][6]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][6]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][6]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][6]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][6]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][6]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][4]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][4]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][4]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][6]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][6]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][6]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][7]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][7]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][7]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh->elmTypes[ielm]==ELM_PRISM){
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][1]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][2]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      dmesh->elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][1]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][4]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][4]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][4]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][4]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][4]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][4]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][5]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][5]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][5]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh->elmTypes[ielm]==ELM_PYRAMID){
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][1]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][1]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][1]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][4]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][4]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][4]][2];
      dmesh->elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = dmesh->coords[dmesh->elmNodes[ielm][0]][0];
      c3[0][1] = dmesh->coords[dmesh->elmNodes[ielm][0]][1];
      c3[0][2] = dmesh->coords[dmesh->elmNodes[ielm][0]][2];
      c3[1][0] = dmesh->coords[dmesh->elmNodes[ielm][3]][0];
      c3[1][1] = dmesh->coords[dmesh->elmNodes[ielm][3]][1];
      c3[1][2] = dmesh->coords[dmesh->elmNodes[ielm][3]][2];
      c3[2][0] = dmesh->coords[dmesh->elmNodes[ielm][2]][0];
      c3[2][1] = dmesh->coords[dmesh->elmNodes[ielm][2]][1];
      c3[2][2] = dmesh->coords[dmesh->elmNodes[ielm][2]][2];
      c3[3][0] = dmesh->coords[dmesh->elmNodes[ielm][4]][0];
      c3[3][1] = dmesh->coords[dmesh->elmNodes[ielm][4]][1];
      c3[3][2] = dmesh->coords[dmesh->elmNodes[ielm][4]][2];
      dmesh->elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the elements around a node.
void t_plas::plasdriver_CalcElmsAroundNode(DRIVER_GAMBIT_MESH *dmesh)
{
#define MAXNODELMS 50
  int inod,ielm;

  dmesh->numNodElms = (int*)calloc(dmesh->numNod,sizeof(int));

  dmesh->nodElms = (int**)calloc(dmesh->numNod,sizeof(int*));
  for(inod=0; inod<dmesh->numNod; inod++){
    dmesh->nodElms[inod] = (int*)calloc(MAXNODELMS,sizeof(int));
  }

  for(ielm=0; ielm<dmesh->numElm; ielm++){
    for(inod=0; inod<dmesh->numElmNodes[ielm]; inod++){
      dmesh->nodElms[dmesh->elmNodes[ielm][inod]][dmesh->numNodElms[dmesh->elmNodes[ielm][inod]]] = ielm;
      dmesh->numNodElms[dmesh->elmNodes[ielm][inod]]++;
    }
  }
#undef MAXNODELMS
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the nodal dual cell volumes.
void t_plas::plasdriver_CalcNodalVolumes(DRIVER_GAMBIT_MESH *dmesh)
{
  int inod,ielm;

  dmesh->nodVolumes = (double*)calloc(dmesh->numNod,sizeof(double));

  for(inod=0; inod<dmesh->numNod; inod++){
    for(ielm=0; ielm<dmesh->numNodElms[inod]; ielm++){
      dmesh->nodVolumes[inod] += dmesh->elmVolumes[dmesh->nodElms[inod][ielm]]/dmesh->numNodElms[inod];
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the volume of a tetrahedron.
double t_plas::plasdriver_CalcVolumeTetra(double c[4][3])
{
  double v1[3],v2[3],v3[3],v4[3];

  v1[0] = c[1][0]-c[0][0];
  v1[1] = c[1][1]-c[0][1];
  v1[2] = c[1][2]-c[0][2];
  v2[0] = c[2][0]-c[0][0];
  v2[1] = c[2][1]-c[0][1];
  v2[2] = c[2][2]-c[0][2];

  plas_CalcCrossProduct_3D(v3,v1,v2);

  v4[0] = c[3][0]-c[0][0];
  v4[1] = c[3][1]-c[0][1];
  v4[2] = c[3][2]-c[0][2];

  return (plas_CalcVectScalarProduct(3,v3,v4)/6.0);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function frees Gambit memory.
void t_plas::plasdriver_FreeGambitMemory(DRIVER_GAMBIT_MESH *dmesh)
{
  int inod,ielm,ibnd,ifac;

  free(dmesh->elmTypes);
  free(dmesh->numElmNodes);
  free(dmesh->numBndFaces);

  for(inod=0; inod<dmesh->numNod; inod++){
    free(dmesh->coords[inod]);
    free(dmesh->nodElms[inod]);
  }
  free(dmesh->coords);
  free(dmesh->nodElms);

  for(ielm=0; ielm<dmesh->numElm; ielm++){
    free(dmesh->elmNodes[ielm]);
    free(dmesh->elmNeighbs[ielm]);
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      free(dmesh->elmNorms[ielm][ifac]);
    }
    free(dmesh->elmNorms[ielm]);
  }
  free(dmesh->elmNodes);
  free(dmesh->elmNeighbs);
  free(dmesh->elmNorms);

  for(ibnd=0; ibnd<dmesh->numBnd; ibnd++){
    free(dmesh->bndNames[ibnd]);
    free(dmesh->bndFaces[ibnd]);
    free(dmesh->bndDomElms[ibnd]);
  }
  free(dmesh->bndNames);
  free(dmesh->bndFaces);
  free(dmesh->bndDomElms);
  free(dmesh->numNodElms);
  free(dmesh->numElmFaces);
  free(dmesh->elmVolumes);
  free(dmesh->nodVolumes);
  free(dmesh->bndTypes);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function gets the nodes of a boundary face.
void t_plas::plasdriver_GetFaceNodes(DRIVER_GAMBIT_MESH *dmesh, int elm, int face, int *nodes)
{
  if(dmesh->elmTypes[elm]==ELM_SIMPLEX){
    if(dmesh->numDim==2){
      if(face==0){
        nodes[0] = dmesh->elmNodes[elm][0];
        nodes[1] = dmesh->elmNodes[elm][1];
        nodes[2] = -1;
        nodes[3] = -1;
      } else if(face==1){
        nodes[0] = dmesh->elmNodes[elm][1];
        nodes[1] = dmesh->elmNodes[elm][2];
        nodes[2] = -1;
        nodes[3] = -1;
      } else if(face==2){
        nodes[0] = dmesh->elmNodes[elm][2];
        nodes[1] = dmesh->elmNodes[elm][0];
        nodes[2] = -1;
        nodes[3] = -1;
      }
    } else if(dmesh->numDim==3){
      if(face==0){
        nodes[0] = dmesh->elmNodes[elm][1];
        nodes[1] = dmesh->elmNodes[elm][0];
        nodes[2] = dmesh->elmNodes[elm][2];
        nodes[3] = -1;
      } else if(face==1){
        nodes[0] = dmesh->elmNodes[elm][0];
        nodes[1] = dmesh->elmNodes[elm][1];
        nodes[2] = dmesh->elmNodes[elm][3];
        nodes[3] = -1;
      } else if(face==2){
        nodes[0] = dmesh->elmNodes[elm][1];
        nodes[1] = dmesh->elmNodes[elm][2];
        nodes[2] = dmesh->elmNodes[elm][3];
        nodes[3] = -1;
      } else if(face==3){
        nodes[0] = dmesh->elmNodes[elm][2];
        nodes[1] = dmesh->elmNodes[elm][0];
        nodes[2] = dmesh->elmNodes[elm][3];
        nodes[3] = -1;
      }
    }
  } else if(dmesh->elmTypes[elm]==ELM_QUAD){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][3];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = -1;
      nodes[3] = -1;
    }
  } else if(dmesh->elmTypes[elm]==ELM_HEX){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = dmesh->elmNodes[elm][5];
      nodes[3] = dmesh->elmNodes[elm][4];
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][3];
      nodes[2] = dmesh->elmNodes[elm][7];
      nodes[3] = dmesh->elmNodes[elm][5];
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][6];
      nodes[3] = dmesh->elmNodes[elm][7];
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = dmesh->elmNodes[elm][6];
    } else if(face==4){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][2];
      nodes[3] = dmesh->elmNodes[elm][3];
    } else if(face==5){
      nodes[0] = dmesh->elmNodes[elm][4];
      nodes[1] = dmesh->elmNodes[elm][5];
      nodes[2] = dmesh->elmNodes[elm][7];
      nodes[3] = dmesh->elmNodes[elm][6];
    }
  } else if(dmesh->elmTypes[elm]==ELM_PRISM){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = dmesh->elmNodes[elm][3];
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][5];
      nodes[3] = dmesh->elmNodes[elm][4];
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][3];
      nodes[3] = dmesh->elmNodes[elm][5];
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][1];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][4];
      nodes[2] = dmesh->elmNodes[elm][5];
      nodes[3] = -1;
    }
  } else if(dmesh->elmTypes[elm]==ELM_PYRAMID){
    if(face==0){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][3];
      nodes[3] = dmesh->elmNodes[elm][1];
    } else if(face==1){
      nodes[0] = dmesh->elmNodes[elm][0];
      nodes[1] = dmesh->elmNodes[elm][1];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = dmesh->elmNodes[elm][1];
      nodes[1] = dmesh->elmNodes[elm][3];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = dmesh->elmNodes[elm][3];
      nodes[1] = dmesh->elmNodes[elm][2];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = dmesh->elmNodes[elm][2];
      nodes[1] = dmesh->elmNodes[elm][0];
      nodes[2] = dmesh->elmNodes[elm][4];
      nodes[3] = -1;
    }
  }
}


// This file contains functionality to define or read a
// steady-state primary flow field.
// This function initializes a steady-state flow field.
void t_plas::plasdriver_InitFlowField(DRIVER_GAMBIT_MESH *dmesh, int material, DRIVER_FLOW_FIELD *dflow)
{
  // Define pressure, velocity and temperature
  double p = 101325.0;
  double u = 1.0;
  double v = 0.0;
  double w = 0.0;
  double T = 373.124;

  if (material==AIR) {

    // Material properties of air
    dflow->rho = 1.225;
    dflow->mu  = 1.7894e-5;
    dflow->cp  = 1.006;
    dflow->k   = 0.0242;

  }
  else if (material==WATER) {

    // Material properties of water
    dflow->rho = 998.2;
    dflow->mu  = 1.003e-3;
    dflow->cp  = 4.182;
    dflow->k   = 0.6;

  }
  else if (material==NITROGEN) {

    // Material properties of nitrogen
    dflow->rho = 1.25;
    dflow->mu  = (6.5592e-7*pow(T,0.6081))/(1.0+54.715/T);
    dflow->cp  = (6.50+0.001*T)*4.184/(28.01e-3);
    dflow->k   = 2.5*(6.5592e-7*pow(T,0.6081))/(1.0+54.715/T)*((6.50+0.001*T)*4.184/(28.01e-3)-8.314)/28.01;
  }

  // Impose flow variables
  for (int n=0; n<dmesh->numNod; ++n) {
    dflow->p[n]    = p;
    dflow->u[n][0] = u;
    dflow->u[n][1] = v;
    dflow->u[n][2] = w;
    dflow->T[n]    = T;
  }
}


// This file contains read and write routines of the PLaS
// driver program.
// This function reads the driver data from an input file.
void t_plas::plasdriver_ReadDriverDataFile(DRIVER_PARAMETERS *dparam)
{
  int i,ignore_i;
  char fileString[100],errMessage[100],text[100],*ignore_cp;
  FILE *inpFile;

  //***Check if file exists***//

  sprintf(fileString,"./plas.driver");
  if(fopen(fileString,"r")==0){
    sprintf(errMessage,"File \"%s\" was not found.",fileString);
    plas_TerminateOnError(errMessage);
  }

  //***Openm file***//

  inpFile = fopen(fileString,"r");

  //***Read case name***//

  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%s",dparam->gridString);
  printf("   Case name is \"%s\"\n",dparam->gridString);
  ignore_cp = fgets(text,100,inpFile);

  //***Read boundaries***//

  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&dparam->numBnd);
  printf("   Number of boundaries is %d\n",dparam->numBnd);
  dparam->bnd = (int*)calloc(dparam->numBnd,sizeof(int));
  ignore_cp  = fgets(text,100,inpFile);
  for(i=0; i<dparam->numBnd; i++){
    ignore_i = fscanf(inpFile,"%d",&dparam->bnd[i]);
  }
  ignore_cp  = fgets(text,100,inpFile);

  //***Read number of iterations***//

  dparam->numIter = plas_ReadIntParam(inpFile);
  if(dparam->numIter<0){
    sprintf(errMessage,"Bad value for numer of iterations.");
    plas_TerminateOnError(errMessage);
  }
  printf("   Number of iterations is %d\n",dparam->numIter);

  //***Read Eulerian time step size***//

  dparam->dtEul = plas_ReadDoubleParam(inpFile);
  if(dparam->dtEul<=0.0){
    sprintf(errMessage,"Bad value for Eulerian time step size.");
    plas_TerminateOnError(errMessage);
  }
  printf("   Eulerian time step is %f\n",dparam->dtEul);

  //***Read flow medium***//

  dparam->material = plas_ReadDoubleParam(inpFile);
  if(dparam->material!=AIR && dparam->material!=WATER && dparam->material!=NITROGEN){
    sprintf(errMessage,"Bad value for material identifier.");
    plas_TerminateOnError(errMessage);
  }
  if(dparam->material==AIR){
    printf("   Flow medium is AIR\n");
  } else if(dparam->material==WATER){
    printf("   Flow medium is WATER\n");
  } else if(dparam->material==NITROGEN){
    printf("   Flow medium is NITROGEN\n");
  }

  fclose(inpFile);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function reads a Gambit mesh.
void t_plas::plasdriver_ReadGambitNeutralFile(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam)
{
  int g,i,j,idim,inod,ielm,ifac,ibnd,jbnd,dummy,ignore_i;
  int ctr_tri,ctr_tet,ctr_qua,ctr_hex,ctr_pri,ctr_pyr;
  char fileString[100],dummyString[100],errMessage[100],*ignore_cp;
  FILE *inpFile;

  //***Open file***//

  sprintf(fileString,"%s.neu",dparam->gridString);
  if(fopen(fileString,"r")==0){
    printf("Neutral file reader error: File \"%s\" was not found.\n",fileString);
    exit(-1);
  }

  inpFile = fopen(fileString,"r");

  for(i=0; i<6; i++){ignore_cp = fgets(dmesh->text[i],100,inpFile);}

  //***Scan header data***//

  ignore_i  = fscanf(inpFile,"%d",&dmesh->numNod);
  ignore_i  = fscanf(inpFile,"%d",&dmesh->numElm);
  ignore_i  = fscanf(inpFile,"%d",&dmesh->numGrps);
  ignore_i  = fscanf(inpFile,"%d",&dmesh->numBnd);
  ignore_i  = fscanf(inpFile,"%d",&dmesh->numDim);
  ignore_i  = fscanf(inpFile,"%d",&dummy);
  ignore_cp = fgets(dummyString,100,inpFile);

  //***Read nodes***//

  printf("   Number of dimensions is %d\n",dmesh->numDim);
  printf("   Number of nodes is %d\n",dmesh->numNod);
  dmesh->coords = (double**)calloc(dmesh->numNod,sizeof(double*));
  for(inod=0; inod<dmesh->numNod; inod++){
    dmesh->coords[inod] = (double*)calloc(dmesh->numDim,sizeof(double));
  }

  for(i=6; i<8; i++){ignore_cp = fgets(dmesh->text[i],100,inpFile);}

  for(inod=0; inod<dmesh->numNod; inod++){
    ignore_i  = fscanf(inpFile,"%d",&dummy);
    for(idim=0; idim<dmesh->numDim; idim++){
      ignore_i  = fscanf(inpFile,"%lf",&dmesh->coords[inod][idim]);
    }
    ignore_cp = fgets(dummyString,100,inpFile);
  }

  //***Read elements***//

  printf("   Number of elements is %d\n",dmesh->numElm);
  ctr_tri = ctr_tet = ctr_qua = ctr_hex = ctr_pri = ctr_pyr = 0;
  dmesh->elmTypes = (int*)calloc(dmesh->numElm,sizeof(int));
  dmesh->numElmNodes = (int*)calloc(dmesh->numElm,sizeof(int));
  dmesh->elmNodes = (int**)calloc(dmesh->numElm,sizeof(int*));
  dmesh->numElmFaces = (int*)calloc(dmesh->numElm,sizeof(int));
  dmesh->elmNorms = (double***)calloc(dmesh->numElm,sizeof(double**));

  for(i=8; i<10; i++){ignore_cp = fgets(dmesh->text[i],100,inpFile);}

  for(ielm=0; ielm<dmesh->numElm; ielm++){

    ignore_i  = fscanf(inpFile,"%d",&dummy);
    ignore_i  = fscanf(inpFile,"%d",&dummy);

    if(dummy==1 || dummy==3 || dummy==6){
      dmesh->elmTypes[ielm] = ELM_SIMPLEX;
      if(dmesh->numDim==2){
        dmesh->numElmFaces[ielm] = 3;
        ctr_tri++;
      } else if(dmesh->numDim==3){
        dmesh->numElmFaces[ielm] = 4;
        ctr_tet++;
      }
    } else if(dummy==2){
      dmesh->elmTypes[ielm] = ELM_QUAD;
      dmesh->numElmFaces[ielm] = 4;
      ctr_qua++;
    } else if(dummy==4){
      dmesh->elmTypes[ielm] = ELM_HEX;
      dmesh->numElmFaces[ielm] = 6;
      ctr_hex++;
    } else if(dummy==5){
      dmesh->elmTypes[ielm] = ELM_PRISM;
      dmesh->numElmFaces[ielm] = 5;
      ctr_pri++;
    } else if(dummy==7){
      dmesh->elmTypes[ielm] = ELM_PYRAMID;
      dmesh->numElmFaces[ielm] = 5;
      ctr_pyr++;
    }

    dmesh->elmNorms[ielm] = (double**)calloc(dmesh->numElmFaces[ielm],sizeof(double*));
    for(ifac=0; ifac<dmesh->numElmFaces[ielm]; ifac++){
      dmesh->elmNorms[ielm][ifac] = (double*)calloc(dmesh->numDim,sizeof(double));
    }

    ignore_i  = fscanf(inpFile,"%d",&dmesh->numElmNodes[ielm]);
    dmesh->elmNodes[ielm] = (int*)calloc(dmesh->numElmNodes[ielm],sizeof(int));
    for(idim=0; idim<dmesh->numElmNodes[ielm]; idim++){
      ignore_i  = fscanf(inpFile,"%d",&dmesh->elmNodes[ielm][idim]);
      dmesh->elmNodes[ielm][idim]--;
    }
    ignore_cp = fgets(dummyString,100,inpFile);
  }

  if(dmesh->numDim==2){
    printf("      Triangles: %d\n",ctr_tri);
    printf("      Quadrilaterals: %d\n",ctr_qua);
  } else if(dmesh->numDim==3){
    printf("      Tetrahedra: %d\n",ctr_tet);
    printf("      Hexahedra: %d\n",ctr_hex);
    printf("      Prisms: %d\n",ctr_pri);
    printf("      Pyramids: %d\n",ctr_pyr);
  }

  //***Read boundary information***//

  printf("   Number of boundaries is %d\n",dmesh->numBnd);
  dmesh->bndNames = (char**)calloc(dmesh->numBnd,sizeof(char*));
  dmesh->numBndFaces = (int*)calloc(dmesh->numBnd,sizeof(int));
  dmesh->bndFaces = (int**)calloc(dmesh->numBnd,sizeof(int*));
  dmesh->bndDomElms = (int**)calloc(dmesh->numBnd,sizeof(int*));
  dmesh->bndTypes = (int*)calloc(dmesh->numBnd,sizeof(int*));
  if(dparam->numBnd!=dmesh->numBnd){
    sprintf(errMessage,"Number of boundary mismatch between driver.conf and Gambit mesh file.");
    plas_TerminateOnError(errMessage);
  }
  for(i=0; i<dmesh->numBnd; i++){
    dmesh->bndTypes[i] = dparam->bnd[i];
  }
  free(dparam->bnd);

  for(i=10; i<14; i++){ignore_cp = fgets(dmesh->text[i],100,inpFile);}

  g = dmesh->numGrps;
  do{
    do{
      ignore_cp = fgets(dummyString,100,inpFile);
    }while(strncmp(dummyString,"ENDOFSECTION",12)!=0);
  }while(--g>0);

  for(ibnd=0; ibnd<dmesh->numBnd; ibnd++){
    ignore_cp = fgets(dmesh->text[i],100,inpFile);
    i++;
    dmesh->bndNames[ibnd] = (char*)calloc(33,sizeof(char));
    for(j=0; j<32; j++){
      ignore_i  = fscanf(inpFile,"%c",&dmesh->bndNames[ibnd][j]);
    }

    //***Align boundary name***//

    i = 0;
    while(dmesh->bndNames[ibnd][i]==32){
      i++;
    }
    j = 0;
    while(i<32){
      dmesh->bndNames[ibnd][j] = dmesh->bndNames[ibnd][i];
      i++;
      j++;
    }
    dmesh->bndNames[ibnd][j] = '\0';

    ignore_i  = fscanf(inpFile,"%d",&dummy);

    if(dummy!=1){
      printf("Neutral file reader error: Boundary ITYPE is not element/cell.\n");
      exit(-1);
    }

    ignore_i  = fscanf(inpFile,"%d",&dmesh->numBndFaces[ibnd]);
    ignore_i  = fscanf(inpFile,"%d",&dummy);
    ignore_i  = fscanf(inpFile,"%d",&dummy);
    ignore_cp = fgets(dummyString,100,inpFile);
    dmesh->bndFaces[ibnd] = (int*)calloc(dmesh->numBndFaces[ibnd],sizeof(int));
    dmesh->bndDomElms[ibnd] = (int*)calloc(dmesh->numBndFaces[ibnd],sizeof(int));
    for(jbnd=0; jbnd<dmesh->numBndFaces[ibnd]; jbnd++){
      ignore_i  = fscanf(inpFile,"%d",&dmesh->bndDomElms[ibnd][jbnd]);
      dmesh->bndDomElms[ibnd][jbnd]--;
      ignore_i  = fscanf(inpFile,"%d",&dummy);
      ignore_i  = fscanf(inpFile,"%d",&dmesh->bndFaces[ibnd][jbnd]);
      dmesh->bndFaces[ibnd][jbnd]--;
      ignore_cp = fgets(dummyString,100,inpFile);
    }
    ignore_cp = fgets(dmesh->text[i],100,inpFile);
    i++;
  }

  fclose(inpFile);
}


// This file contains read and write routines of the PLaS
// driver program.
// This function writes a Tecplot output of the flow field.
void t_plas::plasdriver_WriteTecplot(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam, DRIVER_FLOW_FIELD *dflow)
{
  int inod,idim,ielm;
  char fileString[100];
  FILE *outpFile;

  sprintf(fileString,"%s.plt",dparam->gridString);
  outpFile = fopen(fileString,"w");

  fprintf(outpFile,"TITLE = \"PLaS driver output for Tecplot\"\n");

  if(dmesh->numDim==2){
    fprintf(outpFile,"VARIABLES = x, y, p, u, v, T, theta\n");
  } else if(dmesh->numDim==3){
    fprintf(outpFile,"VARIABLES = x, y, z, p, u, v, w, T, theta\n");
  }

  if(dmesh->numDim==2){
    fprintf(outpFile,"ZONE T=\"Flow Domain\", N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",dmesh->numNod,dmesh->numElm);
  } else if(dmesh->numDim==3){
    fprintf(outpFile,"ZONE T=\"Flow Domain\", N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",dmesh->numNod,dmesh->numElm);
  }

  for(inod=0; inod<dmesh->numNod; inod++){
    for(idim=0; idim<dmesh->numDim; idim++){
      fprintf(outpFile,"%f ",dmesh->coords[inod][idim]);
    }
    fprintf(outpFile,"%f ",dflow->p[inod]);
    for(idim=0; idim<dmesh->numDim; idim++){
      fprintf(outpFile,"%f ",dflow->u[inod][idim]);
    }
    fprintf(outpFile,"%f 0.\n",dflow->T[inod]);
  }

  for(ielm=0; ielm<dmesh->numElm; ielm++){

    if(dmesh->elmTypes[ielm]==ELM_SIMPLEX){
      if(dmesh->numDim==2){
        fprintf(outpFile,"%d %d %d %d\n",
          dmesh->elmNodes[ielm][0]+1,dmesh->elmNodes[ielm][1]+1,
          dmesh->elmNodes[ielm][2]+1,dmesh->elmNodes[ielm][2]+1);
      } else if(dmesh->numDim==3){
        fprintf(outpFile,"%d %d %d %d %d %d %d %d\n",
          dmesh->elmNodes[ielm][0]+1,dmesh->elmNodes[ielm][1]+1,
          dmesh->elmNodes[ielm][2]+1,dmesh->elmNodes[ielm][2]+1,
          dmesh->elmNodes[ielm][3]+1,dmesh->elmNodes[ielm][3]+1,
          dmesh->elmNodes[ielm][3]+1,dmesh->elmNodes[ielm][3]+1);
      }
    } else if(dmesh->elmTypes[ielm]==ELM_QUAD){
      fprintf(outpFile,"%d %d %d %d\n",
        dmesh->elmNodes[ielm][0]+1,dmesh->elmNodes[ielm][1]+1,
        dmesh->elmNodes[ielm][2]+1,dmesh->elmNodes[ielm][3]+1);
    } else if(dmesh->elmTypes[ielm]==ELM_HEX){
      fprintf(outpFile,"%d %d %d %d %d %d %d %d\n",
        dmesh->elmNodes[ielm][0]+1,dmesh->elmNodes[ielm][1]+1,
        dmesh->elmNodes[ielm][3]+1,dmesh->elmNodes[ielm][2]+1,
        dmesh->elmNodes[ielm][4]+1,dmesh->elmNodes[ielm][5]+1,
        dmesh->elmNodes[ielm][7]+1,dmesh->elmNodes[ielm][6]+1);
    } else if(dmesh->elmTypes[ielm]==ELM_PRISM){
      fprintf(outpFile,"%d %d %d %d %d %d %d %d\n",
        dmesh->elmNodes[ielm][0]+1,dmesh->elmNodes[ielm][2]+1,
        dmesh->elmNodes[ielm][1]+1,dmesh->elmNodes[ielm][1]+1,
        dmesh->elmNodes[ielm][3]+1,dmesh->elmNodes[ielm][5]+1,
        dmesh->elmNodes[ielm][4]+1,dmesh->elmNodes[ielm][4]+1);
    } else if(dmesh->elmTypes[ielm]==ELM_PYRAMID){
      fprintf(outpFile,"%d %d %d %d %d %d %d %d\n",
        dmesh->elmNodes[ielm][0]+1,dmesh->elmNodes[ielm][1]+1,
        dmesh->elmNodes[ielm][3]+1,dmesh->elmNodes[ielm][2]+1,
        dmesh->elmNodes[ielm][4]+1,dmesh->elmNodes[ielm][4]+1,
        dmesh->elmNodes[ielm][4]+1,dmesh->elmNodes[ielm][4]+1);
    }
  }

  fclose(outpFile);
}


/*
 * Below are interfacing functions to PLaS (called during PLaS runtime)
 */


void plasinterface_setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->flowSolver   = FLOWSOLVER_DRIVER;
  fp->numDim       = pt_plas->dmesh.numDim;
  fp->numUnk       = pt_plas->dparam.numUnk;
  fp->numNod       = pt_plas->dmesh.numNod;
  fp->numElm       = pt_plas->dmesh.numElm;
  fp->numBnd       = pt_plas->dmesh.numBnd;
  fp->rhoCont      = pt_plas->dflow.rho;
  fp->muCont       = pt_plas->dflow.mu;
  fp->nuCont       = pt_plas->dflow.mu/pt_plas->dflow.rho;
  fp->cpCont       = pt_plas->dflow.cp;
  fp->kCont        = pt_plas->dflow.k;
  fp->dtEul        = pt_plas->dparam.dtEul;
  fp->domainVolume = pt_plas->dmesh.domainVolume;
  fp->minElmVolume = pt_plas->dmesh.minElmVolume;
  fp->maxElmVolume = pt_plas->dmesh.maxElmVolume;
  fp->writeOutput  = 1;

  fp->time = 0.;
  fp->iter = 0;
}


void plasinterface_setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time += pt_plas->dparam.dtEul;
  fp->iter =  pt_plas->dparam.iter;
}


double plasinterface_getBndFaceRefCoord(int bnd, int bface, int dim)
{
  int faceNodes[4];
  ::pt_plas->plasdriver_GetFaceNodes(&pt_plas->dmesh,pt_plas->dmesh.bndDomElms[bnd][bface],pt_plas->dmesh.bndFaces[bnd][bface],faceNodes);

  return pt_plas->dmesh.coords[faceNodes[0]][dim];
}


double plasinterface_getElmFaceMiddlePoint(int elm, int eface, int dim)
{
  int faceNodes[4];
  ::pt_plas->plasdriver_GetFaceNodes(&pt_plas->dmesh,elm,eface,faceNodes);

  int ctr = 0;
  double coord = 0.;
  for (int ifac=0; ifac<4; ++ifac)
    if (faceNodes[ifac]!=-1) {
      coord += pt_plas->dmesh.coords[faceNodes[ifac]][dim];
      ++ctr;
    }
  coord /= ctr;

  return coord;
}


double plasinterface_getBndFaceNormComp           (int bnd, int face, int dim)    { return pt_plas->dmesh.elmNorms[pt_plas->dmesh.bndDomElms[bnd][face]][pt_plas->dmesh.bndFaces[bnd][face]][dim]; }
double plasinterface_getElmNormComp               (int elm, int eface, int dim)   { return pt_plas->dmesh.elmNorms[elm][eface][dim]; }
double plasinterface_getElmVol                    (int elm)                       { return pt_plas->dmesh.elmVolumes[elm]; }
double plasinterface_getEulerianTimeScale         (int nod)                       { return 0.; }
double plasinterface_getNodCoord                  (int nod, int dim)              { return pt_plas->dmesh.coords[nod][dim]; }
double plasinterface_getNodVol                    (int nod)                       { return pt_plas->dmesh.nodVolumes[nod]; }
double plasinterface_getPerBndOffset              (int bnd, int dim)              { return 0.; }  // TODO: implement
double plasinterface_getPressure                  (int nod)                       { return pt_plas->dflow.p[nod]; }
double plasinterface_getPressureOld               (int nod)                       { return pt_plas->dflow.p[nod]; }
double plasinterface_getTemperature               (int nod)                       { return pt_plas->dflow.T[nod]; }
double plasinterface_getTemperatureOld            (int nod)                       { return pt_plas->dflow.T[nod]; }
double plasinterface_getVelocityComp              (int nod, int dim)              { return pt_plas->dflow.u[nod][dim]; }
double plasinterface_getVelocityCompOld           (int nod, int dim)              { return pt_plas->dflow.u[nod][dim]; }
double plasinterface_getVelocityDerivativeComp    (int nod, int idim, int jdim)   { return 0.; }
double plasinterface_getVelocityDerivativeCompOld (int nod, int idim, int jdim)   { return 0.; }
int    plasinterface_getBndDomElm                 (int bnd, int bface, int dummy) { return pt_plas->dmesh.bndDomElms[bnd][bface]; }
int    plasinterface_getElementType               (int elm)                       { return pt_plas->dmesh.elmTypes[elm]; }
int    plasinterface_getElmNeighbour              (int elm, int eface)            { return pt_plas->dmesh.elmNeighbs[elm][eface]; }
int    plasinterface_getElmNode                   (int elm, int enod)             { return pt_plas->dmesh.elmNodes[elm][enod]; }
int    plasinterface_getNumBndFaces               (int bnd)                       { return pt_plas->dmesh.numBndFaces[bnd]; }
int    plasinterface_getPerBndFlag                (int bnd)                       { return pt_plas->dmesh.bndTypes[bnd]<0?  1:0; }
int    plasinterface_getWallBndFlag               (int bnd)                       { return pt_plas->dmesh.bndTypes[bnd]==1? 1:0; }
int    plasinterface_EndElementSearch             (double *pos)                   { return pt_plas->dmesh.numElm-1; }
int    plasinterface_StartElementSearch           (double *pos)                   { return 0; }
void   plasinterface_setPartitioningData          (PLAS_PART_DATA *part)          { part->numNodPairs = 0; }
void   plasinterface_screenOutput                 (char *text)                    { printf("INFO: %s",text); }
void   plasinterface_screenWarning                (char *text)                    { printf("WARN: %s",text); }

