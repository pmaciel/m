/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_generic.c,v $
 *
 * $Author: glhenni $
 *
 * $Date: 1996/02/27 15:23:05 $
 *
 * $Revision: 3.2 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char *cvs_timergen_id =
  "$Id: md_timer_generic.c,v 3.2 1996/02/27 15:23:05 glhenni Exp $";
#endif

/*
     Clock rolls negative when shifts to 32nd bit
     wrap time is about  = 2147 seconds = about 36 minutes
     clockp = double(1<<30)*4.*1.e-6 is the proper increment
     PARAMETER (CLOCKP = 4294.967296)
*/

#include <time.h>

double second()
/*
 * Returns system cpu and wall clock time in seconds
 */
{
  time_t  itime;
  clock_t num_ticks, clock();
  double  cpu, wall;

  num_ticks=clock();
  if(num_ticks==-1) cpu=-1.;
  else cpu = ((double) num_ticks)/((double) CLOCKS_PER_SEC);
  /* wall= ((double) time(&itime)); */

  return(cpu);
}
