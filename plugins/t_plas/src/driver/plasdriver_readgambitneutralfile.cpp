
#include <cstring>

#include "plas_driver.h"


/*
 * This file contains routines to compute the geometry of the
 * mesh used for the steady-state flow solution.
 *
 * This function reads a Gambit mesh.
 */

void plasdriver_ReadGambitNeutralFile(DRIVER_GAMBIT_MESH *dmesh, DRIVER_PARAMETERS *dparam)
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

