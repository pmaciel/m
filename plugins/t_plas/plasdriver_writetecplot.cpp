

#include "plas_driver.h"


/*
 * This file contains read and write routines of the PLaS
 * driver program.
 *
 * This function writes a Tecplot output of the flow field.
 */

void plasdriver_WriteTecplot(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam, DRIVER_FLOW_FIELD *dflow)
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

