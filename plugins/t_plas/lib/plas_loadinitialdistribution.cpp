
#include "common.h"


/*
 * This file contains routines for the generation of entities,
 * be it to form an initial set of entities or to create
 * entites due to secondary phase boundary conditions.
 *
 * This function loads a distribution of entities from a file.
 */

void plas_LoadInitialDistribution(PLAS_DATA *data, char *inpString)
{
  LOCAL_ENTITY_VARIABLES ent;
  int ient,jent,idim,numIni;
  char buffer[200],*search;
  long fpos;
  FILE *inpFile;

  //***Read data file***//

  inpFile = fopen(inpString,"r");

  /* locate zone starting position ("^ZONE.+(dispersed)")) */
  /* set numIni (I) and data->fp.time (SOLUTIONTIME) */
  /* back to correct position */

  fpos = 0L;
  data->fp.time = 0.;
  while (fgets(buffer,200,inpFile)!=NULL) {
    if (strstr(buffer,"ZONE")==buffer && strstr(buffer,"(dispersed)")!=NULL) {
      if ((search = strstr(buffer," I="))!=NULL)
        sscanf(search+3,"%d",&numIni);
      if ((search = strstr(buffer," SOLUTIONTIME="))!=NULL)
        sscanf(search+14,"%lf",&data->fp.time);
      fpos = ftell(inpFile);
    }
  }
  numIni = (numIni<data->ip.numMaxEnt? numIni:data->ip.numMaxEnt);
  fseek(inpFile,fpos,SEEK_SET);

  //***Allocation of local data structure***//

  plas_AllocateLocalEntityVar(data->fp.numDim,&ent);

  //***Loop over entities to generate***//

  jent = 0;
  for (ient=0; ient<numIni && fgets(buffer,200,inpFile)!=NULL; ++ient) {

    //***Read entity data from data file***//

    if (data->fp.numDim==2) {
      sscanf(buffer,"%lf %lf %lf %lf %lf %lf",
        &ent.pos[0],&ent.pos[1],
        &ent.vel[0],&ent.vel[1],
        &ent.temp,&ent.diam);
    }
    else if (data->fp.numDim==3) {
      sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",
        &ent.pos[0],&ent.pos[1],&ent.pos[2],
        &ent.vel[0],&ent.vel[1],&ent.vel[2],
        &ent.temp,&ent.diam);
    }

    //***Perform element search***//

    plas_SearchDomainParallel(data,&ent);

    //***Initialize entity***//

    if (ent.flag==DFLAG_ENABLED) {
      data->ed[jent].flag = DFLAG_ENABLED;
      data->ed[jent].element = ent.elm;
      data->ed[jent].node = plas_FindNearestElementNode(data,&ent);
      data->ed[jent].diameter = ent.diam;
      data->ed[jent].temperature = ent.temp;
      for(idim=0; idim<data->fp.numDim; idim++){
        data->ed[jent].position[idim] = ent.pos[idim];
        data->ed[jent].velocity[idim] = ent.vel[idim];
      }
      jent++;
    }
  }

  fclose(inpFile);

  //***De-allocation of local data structure***//

  plas_DeallocateLocalEntityVar(&ent);
}

