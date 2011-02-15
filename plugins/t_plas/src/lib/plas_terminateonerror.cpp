
#include "common.h"


/*
 * This file contains all write functionality, to the screen,
 * to output files and to Tecplot.
 *
 * This routine terminates PLaS due to a fatal error. It
 * writes out an error message.
 */

void plas_TerminateOnError(char *errMessage)
{
  char screenMessage[120];

  plasinterface_screenWarning((char*) "*** PLAS FATAL ERROR!\n");
  sprintf(screenMessage,"*** %s\n",errMessage);
  plasinterface_screenWarning(screenMessage);
  plasinterface_screenWarning((char*) "*** Terminating program.\n");

  exit(-1);
}

