/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_icc_msr.c,v $
 *
 * $Author: sahutch $
 *
 * $Date: 1996/01/22 18:04:54 $
 *
 * $Revision: 1.6 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_icc_msr.c,v 1.6 1996/01/22 18:04:54 sahutch Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdio.h>
#include <math.h>

#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void icc_pcgpak(int N, double val[], int bindx[], double l[], int ijl[],
                double u[], int iju[], int itemp[], double dtemp[])

/*******************************************************************************

  Routine to calculate incomplete sparse lower and upper factors Choleski
  (l,ijl), (u,iju) for the sparse VBR matrix (val, bindx).

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Order of linear system.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  bindx:           Index array used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  l:               On output, the sparse lower incomplete triagular factor of A.

  ijl:             Pointer array for sparse lower factor l.

  u:               On output sparse upper incomplete triagular factor of A.

  iju:             Pointer array for sparse upper factor u.

  itemp:           Integer work array of size N.

  dtemp:           Double work array of size N.

*******************************************************************************/

{

  /* local variables */

  int    nlast;
  int    i, j, k, *j_ptrs, nzero_ptr;
  double lij, lik, lii, *a_j_ptrs;
  int    nptrs, nj, nk;
  int    sort_flag;

  /**************************** execution begins ******************************/

  /* set pointers to work space */

  j_ptrs   = &itemp[0];
  a_j_ptrs = &dtemp[0];

  ijl[0] = N + 1;

  /* first element */

  l[0]   = sqrt(val[0]);
  ijl[1] = ijl[0];

  for (i = 1; i < N; i++) {     /* loop through rows */
    AZ_gather_ptrs(i, N, val, bindx, &nptrs, j_ptrs, a_j_ptrs);
    sort_flag = 1;

    if (sort_flag == 1 && nptrs > 1) {
      AZ_sort_ptrs(nptrs, j_ptrs, a_j_ptrs);
    }

    ijl[i+1]  = ijl[i];
    nzero_ptr = ijl[i];

    /* calculate and load nonzero off diagonals */

    for (nj = 0; nj < nptrs; nj++) {

      j = j_ptrs[nj];           /* operate only on nonzeros of a */

      if (j < i) {              /* look in lower section only */

        /* start inner sum */

        lij = a_j_ptrs[nj];

        nlast = ijl[j+1] - ijl[j]; /* nonzeros in row j */

        for (nk = 0; nk < nlast; nk++) {

          k = ijl[ijl[j] + nk]; /* form product only */

          if (k < j) {          /* for nonzeros of l */
            lij -= AZ_get_elm(i, k, l, ijl) * AZ_get_elm(j, k, l, ijl);
          }
        }

        ijl[i+1]++;                /* count nonzero */
        ijl[nzero_ptr] = j;        /* pointer to nonzero */
        l[nzero_ptr++] = lij/l[j]; /* value of nonzero */
      }
    }

    /* calculate and load diagonal */

    lii = val[i];

    nlast = ijl[i+1] - ijl[i];  /* nonzeros in row i */

    for (nk = 0; nk < nlast; nk++) {
      k = ijl[ijl[i] + nk];     /* form product only */
      if (k < i) {              /* for nonzeros of l */
        lik  = AZ_get_elm(i, k, l, ijl);
        lii -= lik*lik;
      }
    }

    l[i] = sqrt(lii);
  }

  AZ_transpose(N, l, ijl, u, iju, itemp);

} /* icc_pcgpak */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_transpose(int N, double l[], int ijl[], double lt[], int ijlt[],
                  int row_counter[])

/*******************************************************************************

  Routine to produce the sparse transpose matrix (lt, ijlt) corresponding to the
  sparse matrix (l, ijl).

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Order of linear system.

  l:               On output, the sparse lower incomplete triagular factor of A.

  ijl:             Pointer array for sparse lower factor l.

  lt:              On output sparse transpose of l.

  ijlt:            Pointer array for sparse transpose lt.

  row_counter:     Integer work array of size N.

*******************************************************************************/

{

  /* local variables */

  int irow, ijl_col;
  int ijlt_col, ijlt_row, ijlt_ptr;
  int nzeros, k , ijl_ptr;

  /**************************** execution begins ******************************/

  /* set diagonal elements and initialize ijlt column pointers*/

  for (irow = 0; irow < N; irow++) {
    lt[irow]          = l[irow];
    ijlt[irow]        = row_counter[irow] = 0;
  }

  /* determine number of nonzeros per row in lt */

  for (irow = 0; irow < N; irow++) {
    ijl_col = ijl[irow];
    nzeros  = ijl[irow+1] - ijl[irow];

    /* count each nonzero column entry for the appropriate row of lt */

    for (k = 0; k < nzeros; k++) {
      ijlt_row = ijl[ijl_col + k];
      ijlt[ijlt_row + 1]++;
    }
  }

  /* now produce pointers to columns of lt */

  ijlt[0] = N+1;

  for (irow = 0; irow < N; irow++) ijlt[irow+1] += ijlt[irow];

  /* fill in nonzero off diagonal elements of lt */

  for (irow = 0; irow < N; irow++) {
    nzeros = ijl[irow+1] - ijl[irow];

    for (k = 0; k < nzeros; k++) {
      ijl_ptr  = ijl[irow] + k;

      ijlt_col = irow;          /* switch columns and rows */
      ijlt_row = ijl[ijl_ptr];

      ijlt_ptr = ijlt[ijlt_row]  + row_counter[ijlt_row];

      ijlt[ijlt_ptr] = ijlt_col;
      lt[ijlt_ptr]   = l[ijl_ptr];

      row_counter[ijlt_row]++;  /* count nonzeros already
                                   loaded in (lt,ijlt) */
    }
  }

} /* transpose */
