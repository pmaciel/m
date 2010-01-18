/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_intel.c,v $
 *
 * $Author: glhenni $
 *
 * $Date: 1996/02/27 15:23:08 $
 *
 * $Revision: 3.5 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char *cvs_timerint_id =
  "$Id: md_timer_intel.c,v 3.5 1996/02/27 15:23:08 glhenni Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/* INTEL timer */

double second()

{

  double dclock();

  return (dclock());

} /* second */
