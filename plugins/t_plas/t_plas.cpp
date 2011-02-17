
#include "mfactory.h"
#include "t_plas.h"

using namespace m;


Register< mtransform,t_plas > mt_plas("-tplas","Particle Lagrangian Solver (PLaS)");


void t_plas::transform(GetPot& o, mmesh& m)
{
  screenOutput("read driver parameters (driver.conf)");
  plasdriver_ReadDriverDataFile();


  screenOutput("converting m::mmesh to dmesh...");
  M = &m;
  plasdriver_ReadGambitNeutralFile();


  screenOutput("performing geometry calculations: elements around nodes...");
  plasdriver_CalcElmsAroundNode();
  screenOutput("performing geometry calculations: element neighbours...");
  plasdriver_CalcElementNeighbours();
  screenOutput("performing geometry calculations: element normals...");
  plasdriver_CalcElementNormals();
  screenOutput("performing geometry calculations: element volumes...");
  plasdriver_CalcElementVolumes();
  screenOutput("performing geometry calculations: nodal volumes...");
  plasdriver_CalcNodalVolumes();


  for (int i=0; i<dmesh.numElm; ++i){
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


  dparam.numUnk = M->d()+2;


  screenOutput("initializing flow field...");
  dflow.p = new double [M->n()];
  dflow.T = new double [M->n()];
  dflow.u = new double*[M->n()];
  for (unsigned i=0; i<M->n(); ++i)
    dflow.u[i] = new double[M->d()];
  plasdriver_InitFlowField(dparam.material);


  screenOutput("initializing PLaS...");
  plas::initialize();
  screenOutput("initializing PLaS.");


  screenOutput("perform PLaS iterations...");
  for (dparam.iter=1; dparam.iter<=dparam.numIter; ++dparam.iter) {
    screenOutput("iterate...");
    plas::run();
  }
  screenOutput("perform PLaS iterations.");


  screenOutput("terminating PLaS...");
  ///plasdriver_FreeGambitMemory();
  for (unsigned i=0; i<M->n(); ++i)
    delete[] dflow.u[i];
  delete[] dflow.p;
  delete[] dflow.u;
  delete[] dflow.T;
  screenOutput("terminating PLaS.");
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
void t_plas::plasdriver_CalcElementNeighbours()
{
  int neighbourFound,jnod,knod,lnod,ielm,jelm,kelm,lelm,ifac,faceNodes[4];

  dmesh.elmNeighbs = new int*[dmesh.numElm];
  for(ielm=0; ielm<dmesh.numElm; ielm++){
    dmesh.elmNeighbs[ielm] = new int[dmesh.numElmFaces[ielm]];
    for(ifac=0; ifac<dmesh.numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(ielm,ifac,faceNodes);
      neighbourFound = 0;
      for(jnod=0; jnod<dmesh.numNodElms[faceNodes[0]]; jnod++){
        jelm = dmesh.nodElms[faceNodes[0]][jnod];
        for(knod=0; knod<dmesh.numNodElms[faceNodes[1]]; knod++){
          kelm = dmesh.nodElms[faceNodes[1]][knod];
          if(M->d()==2){
            if(jelm==kelm && jelm!=ielm){neighbourFound = 1;}
          }else if(M->d()==3){
            for(lnod=0; lnod<dmesh.numNodElms[faceNodes[2]]; lnod++){
              lelm = dmesh.nodElms[faceNodes[2]][lnod];
              if(jelm==kelm && jelm==lelm && jelm!=ielm){neighbourFound = 1;}
              if(neighbourFound==1){break;}
            }
          }
          if(neighbourFound==1){break;}
        }
        if(neighbourFound==1){break;}
      }
      if(neighbourFound==1){
        dmesh.elmNeighbs[ielm][ifac] = jelm;
      } else{
        dmesh.elmNeighbs[ielm][ifac] = -1;
      }
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element normals.
void t_plas::plasdriver_CalcElementNormals()
{
  int ielm,ifac,fnodes[4];

  for(ielm=0; ielm<dmesh.numElm; ielm++){
    for(ifac=0; ifac<dmesh.numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(ielm,ifac,fnodes);

      if(dmesh.elmTypes[ielm]==ELM_SIMPLEX && M->d()==2){

        dmesh.elmNorms[ielm][ifac][0] = M->vv[1][fnodes[0]] - M->vv[1][fnodes[1]];
        dmesh.elmNorms[ielm][ifac][1] = M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]];

      } else if((dmesh.elmTypes[ielm]==ELM_SIMPLEX && M->d()==3)
                || (dmesh.elmTypes[ielm]==ELM_PRISM && ifac>2)
                || (dmesh.elmTypes[ielm]==ELM_PYRAMID && ifac>0)){

        dmesh.elmNorms[ielm][ifac][0] =
          0.5*((M->vv[1][fnodes[2]] - M->vv[1][fnodes[0]])
              *(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
              -(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
              *(M->vv[2][fnodes[2]] - M->vv[2][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][1] =
          0.5*((M->vv[2][fnodes[2]] - M->vv[2][fnodes[0]])
              *(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
              -(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
              *(M->vv[0][fnodes[2]] - M->vv[0][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][2] =
          0.5*((M->vv[0][fnodes[2]] - M->vv[0][fnodes[0]])
              *(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
              -(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
              *(M->vv[1][fnodes[2]] - M->vv[1][fnodes[0]]));

      } else if(dmesh.elmTypes[ielm]==ELM_QUAD){

        dmesh.elmNorms[ielm][ifac][0] = 2.0*(M->vv[1][fnodes[0]] - M->vv[1][fnodes[1]]);
        dmesh.elmNorms[ielm][ifac][1] = 2.0*(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]]);

      } else if(dmesh.elmTypes[ielm]==ELM_HEX
                || (dmesh.elmTypes[ielm]==ELM_PRISM && ifac<=2)
                || (dmesh.elmTypes[ielm]==ELM_PYRAMID && ifac==0)){

        dmesh.elmNorms[ielm][ifac][0] =
            ((M->vv[1][fnodes[3]] - M->vv[1][fnodes[0]])
            *(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
            -(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
            *(M->vv[2][fnodes[3]] - M->vv[2][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][1] =
            ((M->vv[2][fnodes[3]] - M->vv[2][fnodes[0]])
            *(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
            -(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
            *(M->vv[0][fnodes[3]] - M->vv[0][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][2] =
            ((M->vv[0][fnodes[3]] - M->vv[0][fnodes[0]])
            *(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
            -(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
            *(M->vv[1][fnodes[3]] - M->vv[1][fnodes[0]]));
      }
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element volumes.
void t_plas::plasdriver_CalcElementVolumes()
{
  int ielm;
  double c2[3][2],c3[4][3];

  dmesh.elmVolumes = new double[dmesh.numElm];

  for(ielm=0; ielm<dmesh.numElm; ielm++){

    if(dmesh.elmTypes[ielm]==ELM_SIMPLEX){
      if(M->d()==2){
        c2[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
        c2[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
        c2[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
        c2[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
        c2[2][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
        c2[2][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
        dmesh.elmVolumes[ielm] = plasdriver_CalcAreaTriangle(c2);
      } else if(M->d()==3){
        c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
        c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
        c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
        c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
        c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
        c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][1]];
        c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
        c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
        c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][2]];
        c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
        c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
        c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
        dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      }
    } else if(dmesh.elmTypes[ielm]==ELM_QUAD){
      c2[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c2[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c2[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
      c2[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
      c2[2][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
      c2[2][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
      dmesh.elmVolumes[ielm] = plasdriver_CalcAreaTriangle(c2);
      c2[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c2[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c2[1][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
      c2[1][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
      c2[2][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c2[2][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcAreaTriangle(c2);
    } else if(dmesh.elmTypes[ielm]==ELM_HEX){
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][1]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][6]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][6]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][6]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][2]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][6]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][6]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][6]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][6]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][6]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][6]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][4]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][4]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][4]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][6]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][6]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][6]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][7]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][7]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][7]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh.elmTypes[ielm]==ELM_PRISM){
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][1]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][2]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][1]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][4]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][4]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][4]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][4]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][4]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][4]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][5]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][5]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][5]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh.elmTypes[ielm]==ELM_PYRAMID){
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][1]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][1]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][1]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][4]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][4]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][4]];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][dmesh.elmNodes[ielm][0]];
      c3[0][1] = M->vv[1][dmesh.elmNodes[ielm][0]];
      c3[0][2] = M->vv[2][dmesh.elmNodes[ielm][0]];
      c3[1][0] = M->vv[0][dmesh.elmNodes[ielm][3]];
      c3[1][1] = M->vv[1][dmesh.elmNodes[ielm][3]];
      c3[1][2] = M->vv[2][dmesh.elmNodes[ielm][3]];
      c3[2][0] = M->vv[0][dmesh.elmNodes[ielm][2]];
      c3[2][1] = M->vv[1][dmesh.elmNodes[ielm][2]];
      c3[2][2] = M->vv[2][dmesh.elmNodes[ielm][2]];
      c3[3][0] = M->vv[0][dmesh.elmNodes[ielm][4]];
      c3[3][1] = M->vv[1][dmesh.elmNodes[ielm][4]];
      c3[3][2] = M->vv[2][dmesh.elmNodes[ielm][4]];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the elements around a node.
void t_plas::plasdriver_CalcElmsAroundNode()
{
#define MAXNODELMS 50
  int inod,ielm;

  dmesh.numNodElms = new int[M->n()];

  dmesh.nodElms = new int*[M->n()];
  for (unsigned inod=0; inod<M->n(); inod++) {
    dmesh.nodElms[inod] = new int[MAXNODELMS];
  }

  for(ielm=0; ielm<dmesh.numElm; ielm++){
    for(inod=0; inod<dmesh.numElmNodes[ielm]; inod++){
      dmesh.nodElms[dmesh.elmNodes[ielm][inod]][dmesh.numNodElms[dmesh.elmNodes[ielm][inod]]] = ielm;
      dmesh.numNodElms[dmesh.elmNodes[ielm][inod]]++;
    }
  }
#undef MAXNODELMS
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the nodal dual cell volumes.
void t_plas::plasdriver_CalcNodalVolumes()
{
  int ielm;

  dmesh.nodVolumes = new double[M->n()];

  for (unsigned inod=0; inod<M->n(); inod++) {
    for(ielm=0; ielm<dmesh.numNodElms[inod]; ielm++){
      dmesh.nodVolumes[inod] += dmesh.elmVolumes[dmesh.nodElms[inod][ielm]]/dmesh.numNodElms[inod];
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
// This function gets the nodes of a boundary face.
void t_plas::plasdriver_GetFaceNodes(int elm, int face, int *nodes)
{
  if(dmesh.elmTypes[elm]==ELM_SIMPLEX){
    if(M->d()==2){
      if(face==0){
        nodes[0] = dmesh.elmNodes[elm][0];
        nodes[1] = dmesh.elmNodes[elm][1];
        nodes[2] = -1;
        nodes[3] = -1;
      } else if(face==1){
        nodes[0] = dmesh.elmNodes[elm][1];
        nodes[1] = dmesh.elmNodes[elm][2];
        nodes[2] = -1;
        nodes[3] = -1;
      } else if(face==2){
        nodes[0] = dmesh.elmNodes[elm][2];
        nodes[1] = dmesh.elmNodes[elm][0];
        nodes[2] = -1;
        nodes[3] = -1;
      }
    } else if(M->d()==3){
      if(face==0){
        nodes[0] = dmesh.elmNodes[elm][1];
        nodes[1] = dmesh.elmNodes[elm][0];
        nodes[2] = dmesh.elmNodes[elm][2];
        nodes[3] = -1;
      } else if(face==1){
        nodes[0] = dmesh.elmNodes[elm][0];
        nodes[1] = dmesh.elmNodes[elm][1];
        nodes[2] = dmesh.elmNodes[elm][3];
        nodes[3] = -1;
      } else if(face==2){
        nodes[0] = dmesh.elmNodes[elm][1];
        nodes[1] = dmesh.elmNodes[elm][2];
        nodes[2] = dmesh.elmNodes[elm][3];
        nodes[3] = -1;
      } else if(face==3){
        nodes[0] = dmesh.elmNodes[elm][2];
        nodes[1] = dmesh.elmNodes[elm][0];
        nodes[2] = dmesh.elmNodes[elm][3];
        nodes[3] = -1;
      }
    }
  } else if(dmesh.elmTypes[elm]==ELM_QUAD){
    if(face==0){
      nodes[0] = dmesh.elmNodes[elm][0];
      nodes[1] = dmesh.elmNodes[elm][1];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==1){
      nodes[0] = dmesh.elmNodes[elm][1];
      nodes[1] = dmesh.elmNodes[elm][2];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = dmesh.elmNodes[elm][2];
      nodes[1] = dmesh.elmNodes[elm][3];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = dmesh.elmNodes[elm][3];
      nodes[1] = dmesh.elmNodes[elm][0];
      nodes[2] = -1;
      nodes[3] = -1;
    }
  } else if(dmesh.elmTypes[elm]==ELM_HEX){
    if(face==0){
      nodes[0] = dmesh.elmNodes[elm][0];
      nodes[1] = dmesh.elmNodes[elm][1];
      nodes[2] = dmesh.elmNodes[elm][5];
      nodes[3] = dmesh.elmNodes[elm][4];
    } else if(face==1){
      nodes[0] = dmesh.elmNodes[elm][1];
      nodes[1] = dmesh.elmNodes[elm][3];
      nodes[2] = dmesh.elmNodes[elm][7];
      nodes[3] = dmesh.elmNodes[elm][5];
    } else if(face==2){
      nodes[0] = dmesh.elmNodes[elm][3];
      nodes[1] = dmesh.elmNodes[elm][2];
      nodes[2] = dmesh.elmNodes[elm][6];
      nodes[3] = dmesh.elmNodes[elm][7];
    } else if(face==3){
      nodes[0] = dmesh.elmNodes[elm][2];
      nodes[1] = dmesh.elmNodes[elm][0];
      nodes[2] = dmesh.elmNodes[elm][4];
      nodes[3] = dmesh.elmNodes[elm][6];
    } else if(face==4){
      nodes[0] = dmesh.elmNodes[elm][1];
      nodes[1] = dmesh.elmNodes[elm][0];
      nodes[2] = dmesh.elmNodes[elm][2];
      nodes[3] = dmesh.elmNodes[elm][3];
    } else if(face==5){
      nodes[0] = dmesh.elmNodes[elm][4];
      nodes[1] = dmesh.elmNodes[elm][5];
      nodes[2] = dmesh.elmNodes[elm][7];
      nodes[3] = dmesh.elmNodes[elm][6];
    }
  } else if(dmesh.elmTypes[elm]==ELM_PRISM){
    if(face==0){
      nodes[0] = dmesh.elmNodes[elm][0];
      nodes[1] = dmesh.elmNodes[elm][1];
      nodes[2] = dmesh.elmNodes[elm][4];
      nodes[3] = dmesh.elmNodes[elm][3];
    } else if(face==1){
      nodes[0] = dmesh.elmNodes[elm][1];
      nodes[1] = dmesh.elmNodes[elm][2];
      nodes[2] = dmesh.elmNodes[elm][5];
      nodes[3] = dmesh.elmNodes[elm][4];
    } else if(face==2){
      nodes[0] = dmesh.elmNodes[elm][2];
      nodes[1] = dmesh.elmNodes[elm][0];
      nodes[2] = dmesh.elmNodes[elm][3];
      nodes[3] = dmesh.elmNodes[elm][5];
    } else if(face==3){
      nodes[0] = dmesh.elmNodes[elm][0];
      nodes[1] = dmesh.elmNodes[elm][2];
      nodes[2] = dmesh.elmNodes[elm][1];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = dmesh.elmNodes[elm][3];
      nodes[1] = dmesh.elmNodes[elm][4];
      nodes[2] = dmesh.elmNodes[elm][5];
      nodes[3] = -1;
    }
  } else if(dmesh.elmTypes[elm]==ELM_PYRAMID){
    if(face==0){
      nodes[0] = dmesh.elmNodes[elm][0];
      nodes[1] = dmesh.elmNodes[elm][2];
      nodes[2] = dmesh.elmNodes[elm][3];
      nodes[3] = dmesh.elmNodes[elm][1];
    } else if(face==1){
      nodes[0] = dmesh.elmNodes[elm][0];
      nodes[1] = dmesh.elmNodes[elm][1];
      nodes[2] = dmesh.elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = dmesh.elmNodes[elm][1];
      nodes[1] = dmesh.elmNodes[elm][3];
      nodes[2] = dmesh.elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = dmesh.elmNodes[elm][3];
      nodes[1] = dmesh.elmNodes[elm][2];
      nodes[2] = dmesh.elmNodes[elm][4];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = dmesh.elmNodes[elm][2];
      nodes[1] = dmesh.elmNodes[elm][0];
      nodes[2] = dmesh.elmNodes[elm][4];
      nodes[3] = -1;
    }
  }
}


// This file contains functionality to define or read a
// steady-state primary flow field.
// This function initializes a steady-state flow field.
void t_plas::plasdriver_InitFlowField(int material)
{
  // Define pressure, velocity and temperature
  double p = 101325.0;
  double u = 1.0;
  double v = 0.0;
  double w = 0.0;
  double T = 373.124;

  if (material==AIR) {

    // Material properties of air
    dflow.rho = 1.225;
    dflow.mu  = 1.7894e-5;
    dflow.cp  = 1.006;
    dflow.k   = 0.0242;

  }
  else if (material==WATER) {

    // Material properties of water
    dflow.rho = 998.2;
    dflow.mu  = 1.003e-3;
    dflow.cp  = 4.182;
    dflow.k   = 0.6;

  }
  else if (material==NITROGEN) {

    // Material properties of nitrogen
    dflow.rho = 1.25;
    dflow.mu  = (6.5592e-7*pow(T,0.6081))/(1.0+54.715/T);
    dflow.cp  = (6.50+0.001*T)*4.184/(28.01e-3);
    dflow.k   = 2.5*(6.5592e-7*pow(T,0.6081))/(1.0+54.715/T)*((6.50+0.001*T)*4.184/(28.01e-3)-8.314)/28.01;
  }

  // Impose flow variables
  for (unsigned n=0; n<M->n(); ++n) {
    dflow.p[n]    = p;
    dflow.u[n][0] = u;
    dflow.u[n][1] = v;
    dflow.u[n][2] = w;
    dflow.T[n]    = T;
  }
}


// This file contains read and write routines of the PLaS
// driver program.
// This function reads the driver data from an input file.
void t_plas::plasdriver_ReadDriverDataFile()
{
  int ignore_i;
  char text[100],*ignore_cp;
  FILE *inpFile;


  screenOutput("check if file \"plas.driver\" exists...");
  if (!fopen("plas.driver","r"))
    plas_TerminateOnError("file \"plas.driver\" was not found");


  screenOutput("open file...");
  inpFile = fopen("plas.driver","r");


  screenOutput("read case name...");
  ignore_cp = fgets(text,100,inpFile);
  ignore_cp = fgets(text,100,inpFile); //ignore_i  = fscanf(inpFile,"%s",dparam.gridString);
  ignore_cp = fgets(text,100,inpFile);


  screenOutput("read boundaries...");
  ignore_cp = fgets(text,100,inpFile);
  ignore_i  = fscanf(inpFile,"%d",&dparam.numBnd);
  dparam.bnd = new int[dparam.numBnd];
  ignore_cp  = fgets(text,100,inpFile);
  for (int i=0; i<dparam.numBnd; ++i)
    ignore_i = fscanf(inpFile,"%d",&dparam.bnd[i]);
  ignore_cp = fgets(text,100,inpFile);


  screenOutput("read number of iterations...");
  dparam.numIter = plas_ReadIntParam(inpFile);
  if (dparam.numIter<0)
    plas_TerminateOnError("bad value for numer of iterations.");


  screenOutput("read Eulerian time step size...");
  dparam.dtEul = plas_ReadDoubleParam(inpFile);
  if (dparam.dtEul<=0.)
    plas_TerminateOnError("bad value for Eulerian time step size");


  screenOutput("read flow medium...");
  dparam.material = plas_ReadDoubleParam(inpFile);
  if (dparam.material!=AIR && dparam.material!=WATER && dparam.material!=NITROGEN)
    plas_TerminateOnError("bad value for material identifier.");


  screenOutput("close file...");
  fclose(inpFile);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function reads a Gambit mesh.
void t_plas::plasdriver_ReadGambitNeutralFile()
{
#if 0
  m_zinner.resize(M->z());
  m_zbound.resize(M->z());
  for (unsigned i=0; i<M->z(); ++i) {
    if      (M->d(i)==M->d())   { m_zinner[i].nelem = M->e(i); m_zinner[i].etype = M->vz[i].t; }
    else if (M->d(i)==M->d()-1) { m_zbound[i].nelem = M->e(i); m_zbound[i].etype = M->vz[i].t; }
  }
#endif


  // general information
  dmesh.numElm = 0;
  for (unsigned i=0; i<M->z(); ++i)
    dmesh.numElm += M->d(i)==M->d()?   M->e(i):0;

  // set volume elements
  dmesh.elmTypes    = new int     [dmesh.numElm];
  dmesh.numElmNodes = new int     [dmesh.numElm];
  dmesh.elmNodes    = new int*    [dmesh.numElm];
  dmesh.numElmFaces = new int     [dmesh.numElm];
  dmesh.elmNorms    = new double**[dmesh.numElm];
  for (unsigned iz=0, ie=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()) {
      for (unsigned i=ie; i<M->e(iz)+ie; ++i) {
        switch(M->vz[iz].t) {
          case (FETRIANGLE):      dmesh.elmTypes[i] = ELM_SIMPLEX; dmesh.numElmFaces[i] = 3; dmesh.numElmNodes[i] = 3; break;
          case (FEQUADRILATERAL): dmesh.elmTypes[i] = ELM_QUAD;    dmesh.numElmFaces[i] = 4; dmesh.numElmNodes[i] = 4; break;
          case (FETETRAHEDRON):   dmesh.elmTypes[i] = ELM_SIMPLEX; dmesh.numElmFaces[i] = 4; dmesh.numElmNodes[i] = 4; break;
          case (FEBRICK):         dmesh.elmTypes[i] = ELM_HEX;     dmesh.numElmFaces[i] = 6; dmesh.numElmNodes[i] = 8; break;
          case (PRISM3):          dmesh.elmTypes[i] = ELM_PRISM;   dmesh.numElmFaces[i] = 5; dmesh.numElmNodes[i] = 5; break;
          case (PYRAMID4):        dmesh.elmTypes[i] = ELM_PYRAMID; dmesh.numElmFaces[i] = 5; dmesh.numElmNodes[i] = 5; break;
          default:                dmesh.elmTypes[i] = 0;
        }
        if (dmesh.elmTypes[i]) {
          dmesh.elmNorms[i] = new double*[dmesh.numElmFaces[i]]; for (int j=0; j<dmesh.numElmFaces[i]; ++j) dmesh.elmNorms[i][j] = new double[M->d()];
          dmesh.elmNodes[i] = new int    [dmesh.numElmNodes[i]]; for (int j=0; j<dmesh.numElmNodes[i]; ++j) dmesh.elmNodes[i][j] = M->vz[iz].e2n[i-ie].n[j];
        }
      }
      ie += M->e(iz);
    }
  }


  // mark boundary nodes
  std::vector< bool > nod_isbnd(M->n(),false);
  for (unsigned iz=0; iz<M->z(); ++iz) if (M->d(iz)==M->d()-1) {
    for (unsigned ie=0; ie<M->e(iz); ++ie) {
      std::vector< unsigned >::const_iterator in;
      for (in = M->vz[iz].e2n[ie].n.begin(); in!=M->vz[iz].e2n[ie].n.end(); ++in)
        nod_isbnd[*in] = true;
    }
  }


#if 0
  // mark boundary elements
  std::vector< std::vector< bool > > elm_isbnd(dmesh.numElm,std::vector< bool >(M->n(),false));
  for (unsigned iz=0, ib=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()-1) {
      for (unsigned jz=0, ie=0; jz<M->z(); ++jz)
        if (M->d(jz)==M->d()) {
          for (unsigned j=0; j<M->e(jz); ++j)
            for (std::vector< unsigned >::const_iterator n=M->vz[jz].e2n[j].n.begin(); n!=M->vz[jz].e2n[j].n.end(); ++n)
              elm_isbnd[ib][ie+j] = elm_isbnd[ib][ie+j] || nod_isbnd[ib][*n];
          ie += M->e(jz);
        }
      ++ib;
    }
  }

  // set number of boundary faces
  dmesh.numBndFaces = new int[dmesh.numBnd];
  for (unsigned iz=0, ib=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()-1) {
    }
  }
#endif



  // set boundary elements
  dmesh.bndFaces    = new int*[M->z()];
  dmesh.bndDomElms  = new int*[M->z()];
  dmesh.bndTypes    = new int [M->z()];
  for (unsigned i=0; i<M->z(); i++)
    dmesh.bndTypes[i] = dparam.bnd[i];
  delete[] dparam.bnd;

  for (unsigned iz=0, ib=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()-1) {

      // set number of boundary faces
      dmesh.numBndFaces[ib] = M->e(iz);

      // (iterate over this boundary's elements)
      typedef std::vector< unsigned > t_element_nodes;
      for (unsigned be=0; be<M->e(iz); ++be) {
        //t_element_nodes &elm_bnd = M->vz[iz].e2n[be].n;



      }

      // set boundary faces "inner" element
      dmesh.bndDomElms[ib] = new int[M->e(iz)];
      for (unsigned jz=0, ie=0; jz<M->z(); ++jz) {
        if (M->d(jz)==M->d()) {

#if 0
          t_element_nodes &a = M->vz[0].e2n[0].n;
#endif
              //        dmesh.bndDomElms[ib][i] = ?;

          ie += M->e(jz);
        }
      }

      // set boundary faces boundary face index (0-based)
      dmesh.bndFaces[ib]   = new int[M->e(iz)];
      for (unsigned i=0; i<M->e(iz); ++i) {
// FIXME     dmesh.bndFaces  [ib][i] = ?;
      }

#if 0
/**
 * comparison of boundary element with "inner cells" connectivity, returning the
 * "inner cell" index and respective opposite node global/local indices.
 * note: it only checks cells which touch a boundary
 * @param[in] e2n "inner" connectivity
 * @param[in] e2n_isbnd marker for cells sitting on the boundary
 * @param[in] ben boundary element to compare
 * @param[out] bcell "inner" corresponding cell (0-based)
 * @param[out] bnode ... cell opposite node, global index (0-based)
 * @param[out] binc ... cell opposite node, local index (0-based)
 */
{
  const std::vector< m::melem > e2n;
  const std::vector< unsigned > ben;
  int *bcell, *bnode, *binc;
{
  const int Nvtfce = 9999;
  bcell = bnode = binc = NULL;

  bool match = false;
  for (unsigned c=0; c<e2n.size() && !match; ++c) {
    if (!elm_isbnd[ib][c])
      continue;
    const std::vector< unsigned >& cell = e2n[c].n;

    // check for "inner" element matching all nodes from the boundary element
    int n_match = 0;
    for (unsigned i=0; i<ben.size(); ++i)
      n_match += std::count(cell.begin(),cell.end(),ben[i])? 1:0;
    match = n_match==Nvtfce;

    // if found, locate which node index (local/global) is not contained
    if (match) {
      *bcell = (int) c;
      for (unsigned i=0; i<cell.size(); ++i)
        if (!std::count(ben.begin(),ben.end(),cell[i])) {
          *bnode = (int) cell[i];
          *binc  = (int)      i;
          break;
        }
    }
  }
  if (!match)
    plas_TerminateOnError("boundary element not found");
}
}
#endif

      ++ib;
    }
  }
}


void t_plas::setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->numDim       = M->d();
  fp->numUnk       = dparam.numUnk;
  fp->numNod       = M->n();
  fp->numElm       = dmesh.numElm;
  fp->numBnd       = M->z();
  fp->rhoCont      = dflow.rho;
  fp->muCont       = dflow.mu;
  fp->nuCont       = dflow.mu/dflow.rho;
  fp->cpCont       = dflow.cp;
  fp->kCont        = dflow.k;
  fp->dtEul        = dparam.dtEul;
  fp->minElmVolume = dmesh.minElmVolume;
  fp->maxElmVolume = dmesh.maxElmVolume;
  fp->writeOutput  = 1;

  fp->time = 0.;
  fp->iter = 0;
}


void t_plas::setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time += dparam.dtEul;
  fp->iter =  dparam.iter;
}


double t_plas::getBndFaceRefCoord(int bnd, int bface, int dim)
{
  int faceNodes[4];
  plasdriver_GetFaceNodes(dmesh.bndDomElms[bnd][bface],dmesh.bndFaces[bnd][bface],faceNodes);
  return M->vv[dim][faceNodes[0]];
}


double t_plas::getElmFaceMiddlePoint(int elm, int eface, int dim)
{
  int faceNodes[4];
  plasdriver_GetFaceNodes(elm,eface,faceNodes);

  int ctr = 0;
  double coord = 0.;
  for (int ifac=0; ifac<4; ++ifac)
    if (faceNodes[ifac]!=-1) {
      coord += M->vv[dim][faceNodes[ifac]];
      ++ctr;
    }
  coord /= ctr;

  return coord;
}
