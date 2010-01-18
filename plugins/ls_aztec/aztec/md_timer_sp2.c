/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_sp2.c,v $
 *
 * $Author: glhenni $
 *
 * $Date: 1996/02/27 15:23:13 $
 *
 * $Revision: 1.4 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char *cvs_timersp2_id =
  "$Id: md_timer_sp2.c,v 1.4 1996/02/27 15:23:13 glhenni Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <sys/time.h>

double second()

{
  double time;
  struct timestruc_t itime;

  gettimer(TIMEOFDAY, &itime);
  time  = (double) itime.tv_sec;
  time += (double) itime.tv_nsec / (double) NS_PER_SEC;
  return time;
}
