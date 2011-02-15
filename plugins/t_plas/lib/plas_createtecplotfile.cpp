
#include "common.h"


/*
 * This file contains all write functionality, to the screen,
 * to output files and to Tecplot.
 *
 * This function creates a Tecplot file of the dispersed
 * phase entities.
 */

void plas_CreateTecplotFile(PLAS_DATA *data, char *outpString)
{
  if (!plas_MpiGetRank()) {

    FILE* outpFile;
    outpFile = fopen(outpString,"w");
    if (data->fp.numDim==2) {
      fprintf(outpFile,"VARIABLES = \"X\" \"Y\" \"U\" \"V\" \"d\" \"T\" \"theta\"\n");
    } else if(data->fp.numDim==3) {
      fprintf(outpFile,"VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"d\" \"T\" \"theta\"\n");
    }
    fclose(outpFile);

  }
}

