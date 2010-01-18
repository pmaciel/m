/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_scaling.c,v $
 *
 * $Author: tuminaro $
 *
 * $Date: 1996/04/26 20:25:30 $
 *
 * $Revision: 1.16 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_scaling.c,v 1.16 1996/04/26 20:25:30 tuminaro Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "az_aztec.h"

/* static functions */

static void calc_blk_diag_Chol(double *val, int *indx, int *bindx, int *rpntr,
                               int *cpntr, int *bpntr, double *L, int *d_indx,
                               int *d_bindx, int *d_rpntr, int *d_bpntr,
                               int *data_org);

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_scale_f(double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
                int bpntr[], double b[], double x[], int options[],
                int data_org[], int proc_config[], int action)

/*******************************************************************************

  Scale the matrix, rhs and initial guess as specified by the user.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Flag which determines whether to scale the matrix or to
                   rescale the solution.

*******************************************************************************/

{

  /* local variables */

  char *yo = "AZ_scale: ";

  /**************************** execution begins ******************************/

  switch (options[AZ_scaling]) {

  case AZ_none:

    /* no scaling , just return */

    break;

  case AZ_Jacobi:

    /* block Jacobi scaling */

    if (action == AZ_SCALE_MAT)
      AZ_block_diagonal_scaling(val, indx, bindx, rpntr, cpntr, bpntr, b,
                                options, data_org, proc_config);
    break;

  case AZ_BJacobi:

    /* block Jacobi scaling */

    if (action == AZ_SCALE_MAT)
      AZ_block_diagonal_scaling(val, indx, bindx, rpntr, cpntr, bpntr, b,
                                options, data_org, proc_config);
    break;

  case AZ_row_sum:

    /* left row-sum scaling */

    if (action == AZ_SCALE_MAT)
      AZ_row_sum_scaling(val, indx, bindx, rpntr, cpntr, bpntr, b, data_org,
                         options);
    break;

  case AZ_sym_diag:
    if (action == AZ_SCALE_MAT)
      AZ_sym_diagonal_scaling(val, bindx, b, data_org, options, x,
      indx, bpntr, rpntr, cpntr);
    else if (action == AZ_RESCALE_SOL)
      AZ_sym_rescale_sl(x, data_org);
    break;

  case AZ_sym_row_sum:
    if (action == AZ_SCALE_MAT) {
      AZ_sym_row_sum_scaling_sl(val, bindx, b, data_org, options, x);
    }
    else if (action == AZ_RESCALE_SOL)
      AZ_sym_rescale_sl(x, data_org);
    break;

  case AZ_sym_BJacobi:

    /* symmetric block Jacobi scaling */
    if (action == AZ_SCALE_MAT)
      AZ_sym_block_diagonal_scaling(val, indx, bindx, rpntr, cpntr, bpntr, b,
                                    options, data_org, proc_config);
#ifdef next_version
    else if (action == AZ_RESCALE_SOL)
      AZ_sym_rescale_vbr(x, data_org);
#endif
    break;

  default:
    (void) fprintf(stderr, "%sERROR: invalid scaling option.\n"
                   "          options[AZ_scaling] = %d\n", yo,
                   options[AZ_scaling]);
  exit(-1);
  }

} /* AZ_scale */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_block_diagonal_scaling(double val[], int indx[], int bindx[],
                               int rpntr[], int cpntr[], int bpntr[],
                               double b[], int options[], int data_org[],
                               int proc_config[])

/*******************************************************************************

  Routine to block Jacobi scale the sparse matrix problem.  Note: this scales
  the matrix and the right-hand side and the resulting matrix is non-symmetric.

  If the matrix is in MSR format, it is treated as point entry (block size 1)
  and standard Jacobi scaling is performed.  Else, if the matrix is in VBR
  format, block Jacobi scaling is performed.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Flag which determines whether to scale the matrix or to
                   rescale the solution.

*******************************************************************************/

{

  /* local variables */

  register int   iblk_row, i, j, k, ival, d_ival, jblk, ib;
  int            bpoff, idoff;
  int            d_bpoff, d_idoff;
  int            m1, n1, ib1, ib2;
  int            ione = 1, itemp;
  int            m, Proc, N;
  int            tsize;
  int            j_last, bindx_row;
  int            max_blk;
  static int    *d3_indx, *d3_bindx, *d3_rpntr, *d3_bpntr;
  static double *d3_inv, *sc_vec;
  double         zero = 0.0, one = 1.0;
  double        *c, *work;
  char           None[2];
  char          *yo = "AZ_block_diagonal_scaling: ";

  /**************************** execution begins ******************************/

  /* initialize */

  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
  m    = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  Proc = proc_config[AZ_node];

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

    /***** MSR Matrix => point Jacobi scaling *****/

    sc_vec = (double *) AZ_manage_memory(N*sizeof(double), AZ_ALLOC,
                                         data_org[AZ_name], "sc_vec",
                                         &itemp);

    if ((options[AZ_pre_calc] >= AZ_reuse) && (itemp == AZ_NEW_ADDRESS)) {
      (void) fprintf(stderr, "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo,
                     data_org[AZ_name]);
      exit(-1);
    }

    if (options[AZ_pre_calc] <= AZ_recalc) {
      for (ib = 0; ib < N; ib++) {
        if (fabs(val[ib]) > DBL_MIN) sc_vec[ib] = 1.0 / val[ib];
        else                         sc_vec[ib] = 1.0;

        val[ib] = 1.0;
        j_last  = bindx[ib+1] - bindx[ib];
        bindx_row = bindx[ib];

        for (j = 0; j < j_last; j++) {
          k       = bindx_row + j;
          val[k] *= sc_vec[ib];
        }
      }
    }

    for (ib = 0; ib < N; ib++)
      b[ib] *= sc_vec[ib];
  }

  else {

    /***** VBR Matrix => block Jacobi scaling *****/

    strcpy(None, "N");

    /* First, compute the block-diagonal inverse (if not already computed) */

    tsize = 0;
    for (i = 0; i < m; i++)
      tsize += (rpntr[i+1] - rpntr[i]) * (cpntr[i+1] - cpntr[i]);

    d3_indx  = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name],
                                           "d3_indx", &itemp);
    d3_bindx = (int *)    AZ_manage_memory(m*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name],
                                           "d3_bindx", &itemp);
    d3_rpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name],
                                           "d3_rpntr", &itemp);
    d3_bpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name],
                                           "d3_bpntr", &itemp);
    d3_inv   = (double *) AZ_manage_memory(tsize*sizeof(double), AZ_ALLOC,
                                           data_org[AZ_name],
                                           "d3_inv", &itemp);

    if ((options[AZ_pre_calc] >= AZ_reuse) && (itemp == AZ_NEW_ADDRESS)) {
      (void) fprintf(stderr, "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo,
                     data_org[AZ_name]);
      exit(-1);
    }

    if (options[AZ_pre_calc] <= AZ_recalc) {
      AZ_calc_blk_diag_inv(val, indx, bindx, rpntr, cpntr, bpntr, d3_inv,
                           d3_indx, d3_bindx, d3_rpntr, d3_bpntr, data_org);

      /* offset of the first block */

      bpoff = *bpntr;
      idoff = *indx;

      d_bpoff = *d3_bpntr;
      d_idoff = *d3_indx;

      /* scale the matrix 'A' */

      max_blk = 0;
      for (i = 0; i < m + data_org[AZ_N_ext_blk] ; i++) {
        if ( cpntr[i+1]-cpntr[i] > max_blk )
          max_blk = cpntr[i+1] - cpntr[i];
      }

      work = (double *) malloc(max_blk*max_blk*sizeof(double));
      if (work == NULL) {
        (void) fprintf(stderr, "%sERROR: not enough memory for diagonal\n"
                       "      scaling. Not able to allocate work\n"
                       "      array. Must run a smaller problem\n", yo);
        exit(-1);
      }

      /* loop over the block rows */

      for (iblk_row = 0; iblk_row < m; iblk_row++) {

        /* find out how many rows are in this block row */

        m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

        /* starting index of current row block */

        ival = indx[bpntr[iblk_row] - bpoff] - idoff;

        /* starting index of current block row for diagonal scaling blocks */

        d_ival  = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;

        /* loop over the (block) columns in the current (block) row */

        for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++){
          jblk = bindx[j];

          /* the starting point column index of the current block */

          ib1 = cpntr[jblk];

          /* ending point column index of the current block */

          ib2 = cpntr[jblk+1];

          /* number of columns in the current block */

          n1    = ib2 - ib1;
          itemp = m1*n1;

          if (jblk == iblk_row) {

            /* diagonal block => set to identity */

            if (m1 != n1) {
              if (Proc == 0) {
                (void) fprintf(stderr, "%sERROR: diagonal blocks are not\n"
                               "       square inside scaling\n", yo);
              }
              exit(-1);
            }

            for (i = 0; i < m1; i++)
              for (k = 0; k < n1; k++)
                if (i == k)
                  val[ival + i + m1*k] = 1.0;
                else
                  val[ival + i + m1*k] = 0.0;
          }

          else {
            if (itemp > max_blk*max_blk) {
              (void) fprintf(stderr, "%sERROR: block size (%d) is too big =>\n",
                             yo, itemp);
              exit(-1);
            }

            /* Dense matrix-matrix multiplication */

#if defined (hp)
            blas_$dgemm(None, None, &m1, &n1, &m1, &one, &d3_inv[d_ival],
                        &m1, &val[ival], &m1, &zero, work, &m1, strlen(None),
                        strlen(None));
            vec_$dcopy( work, val+ival, &itemp);
#else
            dgemm_(None, None, &m1, &n1, &m1, &one, &d3_inv[d_ival], &m1,
                   &val[ival], &m1, &zero, work, &m1, strlen(None),
                   strlen(None));
            dcopy_(&itemp, work, &ione, val+ival, &ione);
#endif
          }
          ival += itemp;
        }
      }

      free((void *) work);
    }

    /* lastly, scale the rhs */

    c = (double *) calloc((unsigned) N, sizeof(double));
    if (c == NULL) {
      (void) fprintf(stderr, "%sERROR: not enough memory for block diagonal\n"
                     "       scaling. Not able to allocate c\n"
                     "       array. Must run a smaller problem\n", yo);
      exit(-1);
    }

    AZ_matvec_mult(d3_inv, d3_indx, d3_bindx, d3_rpntr, d3_rpntr, d3_bpntr, b,
                   c, 1, data_org);

#if defined (hp)
    vec_$dcopy(c, b, &N);
#else
    dcopy_(&N, c, &ione, b, &ione);
#endif

    free((void *) c);
  }

} /* AZ_block_diagonal_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_block_diagonal_scaling(double val[], int indx[], int bindx[],
                                   int rpntr[], int cpntr[], int bpntr[],
                                   double b[], int options[], int data_org[],
                                   int proc_config[])

/*******************************************************************************

  Routine to symmetric block Jacobi scale the sparse matrix problem.  Note: this
  scales the matrix, the solution and the right-hand side.

  If the matrix is in MSR format, it is treated as point entry (block size 1)
  and standard Jacobi scaling is performed.  Else, if the matrix is in VBR
  format, block symmetric Jacobi scaling is performed.

  For the DVBR format, the following system is created from Ax = b:

       (trans(inv(L)) A inv(L) (L x) = trans(inv(L)) b

  where L is the Cholesky factor of the block diagonal portion of A.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Flag which determines whether to scale the matrix or to
                   rescale the solution.

*******************************************************************************/

{

  /* local variables */

  register int   iblk_row, i, j, k, ival, jblk;
  register int   it, jt, icount, iwork, ilow, ib = 0;
  int            bpoff, idoff;
  int            d_bpoff, d_idoff;
  int            m1, n1, ib1, ib2, idL;
  int            ione = 1, itemp, itemp2;
  int            m, Proc;
  int            tsize;
  int            max_blk;
  static int    *d3_indx, *d3_bindx, *d3_rpntr, *d3_bpntr;
  static double *L;
  double         done = 1.0;
  double        *work;
  char           None[2];
  char          *side = "L", *uplo = "L", *transa = "N", *diag = "N";
  char          *yo = "sym_AZ_block_diagonal_scaling: ";

  /**************************** execution begins ******************************/

  printf("Error: Symmetric block scaling not implemented\n");
  exit(-1);

  /* initialize */

  m    = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  Proc = proc_config[AZ_node];

  strcpy(None, "N");

  /*
   * First, compute the block-diagonal Cholesky factorization, its inverse and
   * its inverse's transpose (if not already computed).
   */

  /* malloc temporary space */

  tsize = 0;
  for (i = 0; i < m; i++)
    tsize += (rpntr[i+1] - rpntr[i]) * (cpntr[i+1] - cpntr[i]);

  d3_indx  = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name],
                                         "d3_indx", &itemp);
  d3_bindx = (int *)    AZ_manage_memory(m*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name],
                                         "d3_bindx", &itemp);
  d3_rpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name],
                                         "d3_rpntr", &itemp);
  d3_bpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name],
                                         "d3_bpntr", &itemp);
  L        = (double *) AZ_manage_memory(tsize*sizeof(double), AZ_ALLOC,
                                         data_org[AZ_name], "L",
                                         &itemp);

  if ((options[AZ_pre_calc] >= AZ_reuse) && (itemp == AZ_NEW_ADDRESS)) {
    (void) fprintf(stderr, "%sERROR: Previous scaling not found for matrix "
                   "with\ndata_org[AZ_name] = %d. Check\n"
                   "options[AZ_pre_calc]\n", yo, data_org[AZ_name]);
    exit(-1);
  }

  /*
   * If necessary, calculate the Cholesky factors (L) of the diagonal blocks and
   * store in L and the d3_ pointer vectors.
   */

  if (options[AZ_pre_calc] <= AZ_recalc) {
    calc_blk_diag_Chol(val, indx, bindx, rpntr, cpntr, bpntr, L, d3_indx,
                       d3_bindx, d3_rpntr, d3_bpntr, data_org);

    /* offset of the first block */

    bpoff = *bpntr;
    idoff = *indx;

    d_bpoff = *d3_bpntr;
    d_idoff = *d3_indx;

    /* symmetrically scale the matrix 'A' using the Cholesky factors */

    max_blk = 0;
    for (i = 0; i < m + data_org[AZ_N_ext_blk] ; i++) {
      if ( cpntr[i+1]-cpntr[i] > max_blk )
        max_blk = cpntr[i+1] - cpntr[i];
    }
    work = (double *) malloc(max_blk*max_blk*sizeof(double));
    if (work == NULL) {
      (void) fprintf(stderr, "%sERROR: not enough memory for diagonal\n"
                     "      scaling. Not able to allocate work\n"
                     "      array. Must run a smaller problem\n", yo);
      exit(-1);
    }

    /* loop over the block rows */

    for (iblk_row = 0; iblk_row < m; iblk_row++) {

      /* find out how many rows are in this block row */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row] - bpoff] - idoff;

      /* loop over the (block) columns in the current (block) row */

      for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++){
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = cpntr[jblk];

        /* ending point column index of the current block */

        ib2 = cpntr[jblk+1];

        /* number of columns in the current block */

        n1    = ib2 - ib1;
        itemp = m1*n1;

        if (jblk == iblk_row) {

          /* diagonal block => set to identity */

          if (m1 != n1) {
            if (Proc == 0) {
              (void) fprintf(stderr, "%sERROR: diagonal blocks are not\n"
                             "       square inside scaling\n", yo);
            }
            exit(-1);
          }

          for (i = 0; i < m1; i++) {
            itemp2 = ival + i;
            for (k = 0; k < n1; k++)
              if (i == k)
                val[itemp2 + m1*k] = 1.0;
              else
                val[itemp2 + m1*k] = 0.0;
          }
        }

        else if (jblk < iblk_row) {

          /* lower off-diagonal block */

          if (itemp > max_blk*max_blk) {
            (void) fprintf(stderr, "%sERROR: block size (%d) is too big =>\n",
                           yo, itemp);
            exit(-1);
          }

          /*
           * Fill the work array with the proper block in "A" and take its
           * transpose.
           */

          icount = 0;
          for (it = 0; it < m1; it++) {
            for (jt = 0; jt < n1; jt++) {
              work[icount] = val[ival+icount];
              icount++;
            }
          }

          AZ_dtrans(&m1, &n1, work);

          /* starting index for the diagonal L block for this first operation */

          idL = d3_indx[d3_bpntr[jblk] - d_bpoff] - d_idoff;

          /* perform a backsolve on L*work' = A' to get 'work' array */

          dtrsm_(side, uplo, transa, diag, &m1, &n1, &done, L+idL, &m1, work,
                 &m1, strlen(side), strlen(uplo), strlen(transa), strlen(diag));

          /* need the transpose of the work array */

          AZ_dtrans(&m1, &n1, work);

          /* starting index for the diagonal 'L' block for the next operation */

          idL = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;

          /* perform a backsolve on L*work2 = work */

          dtrsm_(side, uplo, transa, diag, &m1, &n1, &done, L+idL, &m1, work,
                 &m1, strlen(side), strlen(uplo), strlen(transa), strlen(diag));

          /* copy the transpose of this result into the proper block in 'val' */

          iwork  = AZ_get_sym_indx(iblk_row, jblk, indx, bindx, bpntr);
          icount = 0;
          for (i = 0; i < n1; i++)
            for (k = 0; k < m1; k++)
              *(val + iwork + i + k*n1) = *(work + icount++);
        }

        ival += itemp;
      }

      /* scale the RHS */

      idL = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;
      dtrsm_(side, uplo, transa, diag, &m1, &ione, &done, L+idL, &m1, b+ib,
             &m1, strlen(side), strlen(uplo), strlen(transa), strlen(diag));
      ib += m1;
    }

    free((void *) work);

    /*
     * Copy the lower block to their corresponding upper blocks for symmetry.
     */

    /* loop over the block rows */

    for (iblk_row = 0; iblk_row < m; iblk_row++) {

      /* find out how many rows are in this block row */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row] - bpoff] - idoff;

      /* loop over the (block) columns in the current (block) row */

      for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++){
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = cpntr[jblk];

        /* ending point column index of the current block */

        ib2 = cpntr[jblk+1];

        /* number of columns in the current block */

        n1    = ib2 - ib1;
        itemp = m1*n1;

        if (jblk > iblk_row) {

          /*
           * Lower off-diagonal block => copy to corresponding upper block.
           * NOTE: we do not pass in cpntr as the matrix should be symmetric.
           */

          ilow = AZ_get_sym_indx(iblk_row, jblk, indx, bindx, bpntr);

          icount = 0;
          for (i = 0; i < n1; i++)
            for (k = 0; k < m1; k++)
              val[ilow + i + k*n1] = val[ival + icount++];
        }

        ival += itemp;
      }
    }
  }

} /* sym_AZ_block_diagonal_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_row_sum_scaling(double val[], int indx[], int bindx[], int rpntr[],
                        int cpntr[], int bpntr[], double b[], int data_org[],
                        int options[])

/*******************************************************************************

  Routine to left row-sum scale the sparse matrix problem; Note: this scales
  the entire matrix problem Ax = b

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int indx_ptr, irow, iblk_row, jcol, jblk_col, icol_blk, iblk = 0;
  int          num_blk_rows, num_col_blks, num_blk_cols, I, J, t_index;
  int          m, k, N;
  int          j_last, bindx_row;
  double       row_sum = 0.0, sign = 0.0, inv_row_sum, *sc_vec;
  char        *yo = "AZ_row_sum_scaling: ";

  /**************************** execution begins ******************************/

  N      = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sc_vec = (double *) AZ_manage_memory(N*sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], "sc_vec", &k);

  if ((options[AZ_pre_calc] >= AZ_reuse) && (k == AZ_NEW_ADDRESS)) {
    (void) fprintf(stderr, "%sERROR: Previous scaling not found for matrix "
                   "with\ndata_org[AZ_name] = %d. Check\n"
                   "options[AZ_pre_calc]\n", yo,
                   data_org[AZ_name]);
    exit(-1);
  }

  if (options[AZ_pre_calc] <= AZ_recalc) {
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      for (irow = 0; irow < N; irow++) {

        /* get row sum */

        j_last  = bindx[irow+1] - bindx[irow];
        bindx_row = bindx[irow];
        row_sum = fabs(val[irow]);

        for(jcol = 0; jcol < j_last; jcol++) {
          k        = bindx_row + jcol;
          row_sum += fabs(val[k]);
        }

        row_sum = row_sum * sgn(val[irow]);

        if (fabs(row_sum) < DBL_MIN) {
          (void) fprintf(stderr, "%sERROR: Row %d is all zero's\n", yo, irow);
          exit(-1);
        }

        sc_vec[irow] = 1.0 / row_sum;

        /* scale matrix row */

        val[irow] *= sc_vec[irow];
        for (jcol = 0; jcol < j_last; jcol++) {
          k       = bindx_row + jcol;
          val[k] *= sc_vec[irow];
        }
      }
    }

    else {
      m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

      /* loop over the block rows */

      for (iblk_row = 0; iblk_row < m; iblk_row++) {

        /* find out how many rows are in this block row */

        num_blk_rows = rpntr[iblk_row+1] - rpntr[iblk_row];

        /* find out how many block columns are in this block row */

        num_col_blks = bpntr[iblk_row+1] - bpntr[iblk_row];

        /* loop over all the rows in this block row */

        for (irow = 0; irow < num_blk_rows; irow++) {
          I = rpntr[iblk_row] + irow;     /* true matrix row number */

          /* loop over all the column blocks in this block row */

          for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

            /* find out which column block this is */

            icol_blk = bindx[iblk];
            indx_ptr = indx[iblk++];

            /* find out how many columns are in this block */

            num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

            /* loop over all the columns in this block */

            for (jcol = 0; jcol < num_blk_cols; jcol++) {
              J     = cpntr[icol_blk] + jcol;     /* true matrix column */
              t_index = indx_ptr + jcol*num_blk_rows + irow;

              /* diagonal entry => get sign */

              if (I == J) sign = sgn(val[t_index]);

              row_sum += fabs(val[t_index]);
            }
          }

          /* reset the block counter */

          iblk -= num_col_blks;

          if ( fabs(sign) < (1.0 - sqrt(DBL_EPSILON)) ) {
            (void) fprintf(stderr, "%sERROR: sign not set => no diagonal "
                           "entry.\n         inside row_sum.\n", yo);
            exit(-1);
          }

          if (fabs(row_sum) == 0.0) {
            (void) fprintf(stderr,"%sERROR: row %d is all 0's.\n", yo, I);
            exit(-1);
          }

          inv_row_sum = sign / row_sum;
          sc_vec[I]   = inv_row_sum;
          row_sum     = sign = 0.0;

          /* scale the matrix */

          for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

            /* find out which column block this is */

            icol_blk = bindx[iblk];
            indx_ptr = indx[iblk++];

            /* find out how many columns are in this block */

            num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

            /* loop over all the columns in this block */

            for (jcol = 0; jcol < num_blk_cols; jcol++) {
              J           = cpntr[icol_blk] + jcol;     /* true matrix column */
              t_index       = indx_ptr + jcol * num_blk_rows + irow;
              val[t_index] *= inv_row_sum;
            }
          }

          /* reset the block counter */

          iblk -= num_col_blks;
        }

        /* last row in this block row => offset the correction above */

        iblk += num_col_blks;
      }
    }
  }

  for (irow = 0; irow < N; irow++)
    b[irow] *= sc_vec[irow];

} /* AZ_row_sum_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void sym_row_sum_scaling(int proc_config[])

/*******************************************************************************

  Routine to symmetricaly row-sum scale sparse matrix problem; Note: this scales
  the entire matrix problem Ax = b, the routine sym_rescale must be used to
  transform solution back to recover soution to original problem.

  Author:          , SNL,
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  char *yo = "sym_row_sum_scaling: ";

  /**************************** execution begins ******************************/

  if (proc_config[AZ_node]  == 0) {
    (void) fprintf(stderr, "%sWARNING: sym_row_sum_scaling not implemented"
                   "           for VBR matrices\n"
                   "No preconditioning performed\n", yo);
  }

} /* sym_row_sum_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_diagonal_scaling(double val[], int bindx[], double b[],
                        int data_org[], int options[], double x[],
                        int indx[], int bpntr[], int rpntr[], int cpntr[])

/*******************************************************************************

  Routine to symmetrically diagonally scale sparse matrix problem; Note: this
  scales the entire matrix problem Ax = b, the routine sym_rescale must be used
  to transform solution back to recover soution to original problem.

  Author:          John N. Shadid, SNL, 1421 (MSR format)
  =======          Lydie Prevost, SNL 1422 (VBR format) 

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  bindx:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  x:               Current solution vector.

*******************************************************************************/

{

  /* local variables */

  register int j, k, irow, icol;
  int          N, m;
  int          j_last, bindx_row, i;
  double       *sc_vec;
  int count;

  char        *yo = "AZ_sym_diagonal_scaling: ";


  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];
  m    = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];



  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], "sc_vec", &i);
  if ((options[AZ_pre_calc] >= AZ_reuse) && (i == AZ_NEW_ADDRESS)) {
    (void) fprintf(stderr, "%sERROR: Previous scaling not found for matrix "
                   "with\ndata_org[AZ_name] = %d. Check\n"
                   "options[AZ_pre_calc]\n\n", yo,
                   data_org[AZ_name]);
    exit(-1);
  }



  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

    /***** MSR Matrix *****/


    if (options[AZ_pre_calc] <= AZ_recalc) {
      for (irow = 0; irow < N; irow++) {

        /* scale matrix */

        j_last  = bindx[irow+1] - bindx[irow];
        bindx_row = bindx[irow];

        if (fabs(val[irow]) < DBL_MIN) {
          (void) fprintf(stderr, "%sERROR: diagonal of row %d is zero\n", yo,
                         irow);
          exit(-1);
        }

        sc_vec[irow] = 1.0 / sqrt(fabs(val[irow]));

        for (j = 0; j < j_last; j++) {
          k       = bindx_row + j;
          val[k] *= sc_vec[irow];
        }
        val[irow] *= sc_vec[irow];
      }

      /* do right diagonal scaling */

      AZ_exchange_bdry(sc_vec, data_org);

      /* index through rows of matrix */

      for (irow = 0; irow < N; irow++) {
        val[irow] *= sc_vec[irow];

        j_last     = bindx[irow+1] - bindx[irow];
        bindx_row    = bindx[irow];

        for (j = 0; j < j_last; j++) {
          k       = bindx_row + j;
          val[k] *= sc_vec[bindx[k]];
        }
      }
    }

  }

  else {

    /***** VBR Matrix *****/
   
    /* find the diagonal terms */

    if (options[AZ_pre_calc] <= AZ_recalc) {

     /* index through block rows of matrix */
   
         for (irow = 0; irow < m; irow++) {
   
         /* for all the blocks in a row */
   
         for (k = bpntr[irow]; k < bpntr[irow+1]; k++) {
            icol = bindx[k];
   
            count = 0;
   
            for (i = rpntr[irow]; i < rpntr[irow+1]; i++) {
                for (j = cpntr[icol]; j < cpntr[icol+1]; j++) {
                  if ( icol == irow && i ==j ){
                     sc_vec[i] =  1.0 / sqrt(fabs(val[indx[k]+count]));
                  }
                  count++;
               }
            }
   
         }
      }


     /* do left and right diagonal scaling */
   
      AZ_exchange_bdry(sc_vec, data_org);
   
     /* index through rows of matrix */
   
      for (irow = 0; irow < m; irow++) {
   
         /* for all the blocks in a row */
   
         for (k = bpntr[irow]; k < bpntr[irow+1]; k++) {
            icol = bindx[k];
   
            count = 0;
   
            for (i = rpntr[irow]; i < rpntr[irow+1]; i++) {
                for (j = cpntr[icol]; j < cpntr[icol+1]; j++) {
                  val[indx[k]+count] *=  sc_vec[i] * sc_vec[j] ;
                  count++;
               }
            }
   
         }
      }
   
    }
  }

  /* rescale right hand side and solution */

  for (irow = 0; irow < N; irow++) b[irow] *= sc_vec[irow];
  for (irow = 0; irow < N; irow++) x[irow] /= sc_vec[irow];

} /* sym_diagonal_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_row_sum_scaling_sl(double val[], int bindx[], double b[],
                               int data_org[], int options[], double x[])

/*******************************************************************************

  Routine to symmetrically diagonally scale sparse matrix problem; Note: this
  scales the entire matrix problem Ax = b, the routine sym_rescale must be used
  to transform solution back to recover soution to original problem.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  bindx:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  b:               Right hand side of linear system.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  x:               Current solution vector.

*******************************************************************************/

{

  /* local variables */

  register int j, k, irow;
  int          N;
  int          i, jcol, j_last, bindx_row;
  double       row_sum, *sc_vec;
  char        *yo = "AZ_sym_row_sum_scaling_sl: ";

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], "sc_vec", &i);

  if ((options[AZ_pre_calc] >= AZ_reuse) && (i == AZ_NEW_ADDRESS)) {
    (void) fprintf(stderr, "%sERROR: Previous scaling not found for matrix "
                   "with\ndata_org[AZ_name] = %d. Check\n"
                   "options[AZ_pre_calc]\n", yo,
                   data_org[AZ_name]);
    exit(-1);
  }

  if (options[AZ_pre_calc] <= AZ_recalc) {
    for(irow = 0; irow < N; irow++) {

      /* get row sum */

      j_last  = bindx[irow+1] - bindx[irow];
      bindx_row = bindx[irow];
      row_sum = fabs(val[irow]);

      for (jcol = 0; jcol < j_last; jcol++) {
        k        = bindx_row + jcol;
        row_sum += fabs(val[k]);
      }

      if (fabs(row_sum) < DBL_MIN) {
        (void) fprintf(stderr, "%sERROR: Row %d is all zero's\n", yo, irow);
        exit(-1);
      }

      sc_vec[irow] = 1.0 / sqrt(fabs(row_sum));

      /* scale matrix row */

      val[irow] *= sc_vec[irow];

      for (jcol = 0; jcol < j_last; jcol++) {
        k       = bindx_row + jcol;
        val[k] *= sc_vec[irow];
      }
    }

    AZ_exchange_bdry(sc_vec, data_org);

    /* do right diagonal scaling */
    /* index through rows of matrix */

    for (irow = 0; irow < N; irow++) {
      val[irow] *= sc_vec[irow];
      j_last     = bindx[irow+1] - bindx[irow];
      bindx_row    = bindx[irow];

      for (j = 0; j < j_last; j++) {
        k       = bindx_row + j;
        val[k] *= sc_vec[bindx[k]];
      }
    }
  }

  for (irow = 0; irow < N; irow++) b[irow] *= sc_vec[irow];
  for (irow = 0; irow < N; irow++) x[irow] /= sc_vec[irow];

} /* sym_row_sum_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_rescale_sl(double x[], int data_org[])

/*******************************************************************************

  Routine to symmetrically diagonally rescale the sparse matrix problem;
  Note: this rescales the entire matrix problem Ax = b.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Current solution vector.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N, k;
  double      *sc_vec;
  char        *yo = "AZ_sym_rescale_sl: ";

  /**************************** execution begins ******************************/

  if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
      (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) ) {
     (void) fprintf(stderr,"%sWARNING: Matrix type is neither MSR nor VBR\n");
     return;
  }



  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], "sc_vec", &k);

  if (k == AZ_NEW_ADDRESS) {
    (void) fprintf(stderr, "%sWARNING: Scaling vector not found: "
                   "Not rescaling solution\n", yo);
    return;
  }

  for (i = 0; i < N; i++) x[i] = x[i] * sc_vec[i];

  AZ_exchange_bdry(x, data_org);

} /* AZ_sym_rescale_sl */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void calc_blk_diag_Chol(double *val, int *indx, int *bindx, int *rpntr,
                               int *cpntr, int *bpntr, double *L, int *d_indx,
                               int *d_bindx, int *d_rpntr, int *d_bpntr,
                               int *data_org)

/*******************************************************************************

  Routine to calculate the Cholesky factors of the block-diagonal portion of the
  sparse matrix in 'val' and the associated integer pointer vectors. This is
  used for scaling and/or preconditioning.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  L:               Vector containing the upper Cholesky factors of the diagonal
                   blocks.

  d_indx:          The 'indx' array corresponding to the inverse-block
                   diagonals.

  d_bindx:         The 'bindx' array corresponding to the inverse-block
                   diagonals.

  d_rpntr:         The 'rpntr' array corresponding to the inverse-block
                   diagonals.

  d_bpntr:         The 'bpntr' array corresponding to the inverse-block
                   diagonals.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  register int i, j, iblk_row, jblk, icount = 0, iblk_count = 0, ival;
  int          m1, n1, itemp, iL;
  int          m;
  int          bpoff, idoff;
  int          info;
  char        *yo = "calc_blk_diag_Chol: ", *uplo = "L";

  /**************************** execution begins ******************************/

  m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (m == 0) return;

  /* offset of the first block */

  bpoff = *bpntr;
  idoff = *indx;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    /* number of rows in the current row block */

    m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* starting index of current row block */

    ival = indx[bpntr[iblk_row] - bpoff] - idoff;

    /* loop over column block numbers, looking for the diagonal block */

    for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++) {
      jblk = bindx[j];

      /* determine the number of columns in this block */

      n1 = cpntr[jblk+1] - cpntr[jblk];

      itemp = m1*n1;

      if (jblk == iblk_row) {   /* diagonal block */

        /* error check */

        if (n1 != m1) {
          (void) fprintf(stderr, "%sERROR: diagonal blocks are not square.\n",
                         yo);
          exit(-1);
        }
        else {

          /* fill the vectors */

          d_indx[iblk_count]  = icount;
          d_rpntr[iblk_count] = rpntr[iblk_row];
          d_bpntr[iblk_count] = iblk_row;
          d_bindx[iblk_count] = iblk_row;

          for (i = 0; i < itemp; i++) L[icount++] = val[ival + i];

          /* Compute the Cholesky factors for this block */

          iL = d_indx[d_bpntr[iblk_row] - *d_bpntr] - *d_indx;
          dpotrf_(uplo, &m1, L+iL, &m1, &info, strlen(uplo));

          if (info < 0) {
            (void) fprintf(stderr, "%sERROR: argument %d is illegal\n", yo,
                           -info);
            exit(-1);
          }
          else if (info > 0) {
            (void) fprintf(stderr, "%sERROR: the factorization has produced a "
                           "singular L with L[%d][%d] being exactly zero\n",
                           yo, info, info);
            exit(-1);
          }

          iblk_count++;
        }
        break;
      }
      else
        ival += itemp;
    }
  }

  d_indx[iblk_count]  = icount;
  d_rpntr[iblk_count] = rpntr[iblk_row];
  d_bpntr[iblk_count] = iblk_row;

} /* calc_blk_diag_Chol */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#ifdef next_version
void AZ_sym_rescale_vbr(double x[], int data_org[])

/*******************************************************************************

  Routine to symmetrically block-diagonally rescale the sparse matrix problem;
  Note: this rescales the entire matrix problem Ax = b.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Current solution vector.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  int          N, k;
  double      *sc_vec;
  char        *yo = "AZ_sym_rescale_vbr: ";

  /**************************** execution begins ******************************/

  if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) return;

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], "sc_vec", &k);

  if (k == AZ_NEW_ADDRESS) {
    (void) fprintf(stderr, "%sWARNING: Scaling vector not found - "
                   "not rescaling solution\n", yo);
    return;
  }
  /*
    for (i = 0; i < N; i++) x[i] = x[i] * sc_vec[i];

    AZ_exchange_bdry(x, data_org);
    */
} /* AZ_sym_rescale_vbr */
#endif
