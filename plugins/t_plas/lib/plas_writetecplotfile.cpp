
#include "common.h"


/*
 * This file contains all write functionality, to the screen,
 * to output files and to Tecplot.
 *
 * This function appends Tecplot output of the dispersed
 * phase entities into the previously created file
 */

void plas_WriteTecplotFile(PLAS_DATA *data, char *outpString, int iter, double time)
{
  FILE *outpFile;
  int ient, iproc, idim;
  const int numProc = plas_MpiGetNumProc();
  const int irank   = plas_MpiGetRank();

  for (iproc=0; iproc<numProc; ++iproc) {
    if (iproc==irank) {
      outpFile = fopen(outpString,"a");
      if (irank==0) {

        // write header (if no entities are active, write dummy)
        fprintf(outpFile,"ZONE T=\"Entities (dispersed)\" I=%d AUXDATA ITER=\"%d\" SOLUTIONTIME=%12.6f DATAPACKING=POINT\n",data->sd.enabled? data->sd.enabled:1,iter,time);
        if (!data->sd.enabled) {
          for (idim=0; idim<data->fp.numDim*2+3; ++idim)
            fprintf(outpFile,"0 ");
          fprintf(outpFile,"\n");
        }

      }
      for (ient=0; ient<data->ip.numMaxEnt && data->sd.enabled; ++ient) {
        if (data->ed[ient].flag==DFLAG_ENABLED || data->ed[ient].flag==DFLAG_CREATED) {

          PLAS_ENTITY_DATA *e = &(data->ed[ient]);
          for (idim=0; idim<data->fp.numDim; ++idim)
            fprintf(outpFile,"%.12f ",e->position[idim]);
          for (idim=0; idim<data->fp.numDim; ++idim)
            fprintf(outpFile,"%.12f ",e->velocity[idim]);
          fprintf(outpFile,"%.12f %.12f %.12f\n",e->diameter,e->temperature,0.);

        }
      }
      fclose(outpFile);
    }
    plas_MpiBarrier();
  }
}

