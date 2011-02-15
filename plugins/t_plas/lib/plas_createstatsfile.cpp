
#include "common.h"


/*
 * This file contains all write functionality, to the screen,
 * to output files and to Tecplot.
 *
 * This function writes the PLaS statistics to a file.
 */

void plas_CreateStatsFile(PLAS_DATA *data, char *outpString)
{
  FILE *outpFile;
  if (!plas_MpiGetRank()) {
    outpFile = fopen(outpString,"w");
    if (outpFile==NULL)
      return;

    fprintf(outpFile,"Iter\tTime\tEnt\tIn\tOut\tBounc\tColl\tPer\tPass\tLost\tRe_disp\tNu_disp\tdt_Lagr\tsubit\n");
    fclose(outpFile);
  }
}

