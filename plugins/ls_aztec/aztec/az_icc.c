/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_icc.c,v $
 *
 * $Author: tuminaro $
 *
 * $Date: 1996/04/26 20:35:37 $
 *
 * $Revision: 1.8 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_icc.c,v 1.8 1996/04/26 20:35:37 tuminaro Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void icc()

/*******************************************************************************

  C driver for incomplete choleski factor of sparse matrix.

  This routine takes a VBR format sparse matrix and determines the incomplete
  Choleski decomposition. Then it performs the solve step on the traingular
  systems.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Order of linear system.

  NExt:            Number of external variables.

  M:               Number of nonzeros in sparse matrix A.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  x:               On input, contains the current solution to the linear system.

*******************************************************************************/

{

  /* local variables */

  /**************************** execution begins ******************************/

  (void) fprintf(stderr, "WARNING: icc preconditioning not implemented\n");
  (void) fprintf(stderr, "         Try ilu instead.\n");
  (void) fprintf(stderr, "Returning with no preconditioning...\n");

} /* icc */
