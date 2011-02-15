
#include "common.h"


/*
 * This file contains all write functionality, to the screen,
 * to output files and to Tecplot.
 *
 * This function writes the PLaS statistics to a file.
 */

void plas_WriteStatsFile(PLAS_DATA *data, char *outpString, int iter, double time)
{
  FILE *outpFile;
  if (!plas_MpiGetRank()) {
    outpFile = fopen(outpString,"a");
    if (outpFile==NULL)
      return;

    fprintf(outpFile,"%5d\t%6.2f\t%6d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%11.4e\t%11.4e\t%11.4e\t%6.2f\n",
    iter,time,
    data->sd.enabled,data->sd.in,data->sd.out,data->sd.bounce,data->sd.coll,
    data->sd.periodic,data->sd.passed,data->sd.lost,data->sd.reynoldsAvg,data->sd.nusseltAvg,
    data->sd.dtLagrAvg,data->sd.subIterAvg);
    fclose(outpFile);
  }
}

