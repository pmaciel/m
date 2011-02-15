
#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function computes element volumes.
 */

void plasdriver_CalcElementVolumes(DRIVER_GAMBIT_MESH *dmesh)
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

