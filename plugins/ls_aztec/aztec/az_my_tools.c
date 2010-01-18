/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_my_tools.c,v $
 *
 * $Author: Edwin van der Weide $
 *
 * $Date: 1997/07/02 19:30:20 $
 *
 * $Revision: Modified version of AZ_transform (val is not needed anymore) $
 *            Possibility to include block milu preconditioner (more efficient 
 *            in terms of memory, not performance)
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_my_tools.c,v 1.33 1996/06/24 19:30:20 tuminaro Exp $";
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
#include <malloc.h>
#include <float.h>
#include <string.h>
#include "az_aztec.h"

extern int AZ_sys_msg_type;
extern int AZ_using_fortran;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_reorder_matrix_mod(int N_update, int bindx[],
                           int update_index[], int extern_index[], int indx[],
                           int rnptr[], int bnptr[], int N_external,
                           int cnptr[], int option, int mat_type)

/*******************************************************************************

  Reorder the matrix so that it corresponds to the new ordering given by
  'update_index' and 'extern_index'.

  IMPORTANT: This routine assumes that update_index[] contains two sequencies of
  numbers that are ordered but intertwined. For example,

  update_index:  4 5 0 6 1 2 3 7

  seq 1 =>    0   1 2 3

  seq 2 =>4 5   6 7

  Author:          Ray S. Tuminaro, SNL, 1422
  =======
  Modified:        Edwin van der Weide, VKI
                   val does not occur any more in the subroutine

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  rnptr,
  bindx,
  indx,
  bnptr,
  cnptr:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  update_index:    update_index[i] gives the local numbering of global point
                   'update[i]'.

  extern_index:    extern_index[i] gives the local numbering of global point
                   'external[i]'.

  N_external:      Number of external elements on this processor.

  option:

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

*******************************************************************************/

{

  /* local variables */

  int   start, end;
  int   indx_length;
  int  *temp;
  int   i, j;
  char *yo = "AZ_reorder_matrix_mod: ";

  /**************************** execution begins ******************************/

  if (mat_type == AZ_MSR_MATRIX) {
    start = N_update+1;      /* first nonzero offdiag */
    end   = bindx[N_update]; /* last nonzero          */
  }
  else if (mat_type == AZ_VBR_MATRIX) {
    start = 0;               /* first nonzero */
    end   = bnptr[N_update]; /* last nonzero  */

    /* reorder cnptr[] */

    /* 1) convert the cnptr array to give the blk size */

    AZ_convert_ptrs_to_values(cnptr, N_update + N_external);

    /* 2) order the internal blocks.  NOTE: AZ_sortqlists() can only be used for
     * the internal elements as it expects the list to correspond to 2 ordered
     * sequences that are intermixed.
     */

    AZ_sortqlists((char *) cnptr, 0, update_index, N_update, sizeof(int),
                  N_update);

    /* 3) order the external blocks */

    temp = (int *) calloc((unsigned)(N_external + 1), sizeof(int));
    if (temp == NULL) {
      (void) fprintf(stderr,
                     "%sERROR: not enough memory to malloc temporary space\n",
                     yo);
      exit(-1);
    }

    for (i = 0; i < N_external; i++)
      temp[extern_index[i] - N_update] = cnptr[i + N_update];

    for (i = 0; i < N_external; i++) cnptr[i + N_update] = temp[i];
    free((char *) temp);

    /* 4) reconvert cnptr to give pointer information */

    AZ_convert_values_to_ptrs(cnptr, N_update + N_external, 0);
  }
  else {
    (void) fprintf(stderr, "%sERROR: matrix is not MSR or VBR\n", yo);
    exit(-1);
  }

  /*
   * Change column indices (bindx) to reflect new ordering depending
   * on whether or not a point is internal or external.
   */

  for (i = start; i < end; i++) {
    if (bindx[i] < N_update) bindx[i] = update_index[bindx[i]];
    else                     bindx[i] = extern_index[bindx[i] - N_update];
  }

  if (option == AZ_EXTERNS) return;

  /* reorder rows */

  if (mat_type == AZ_MSR_MATRIX) {

    /* We move the rows in four steps:
     *  1) sort the first N_update values of 'val'.
     *  2) sort the first N_update values of 'bindx'.
     *     We do this by first converting the ptrs to values
     *     representing the number of nonzero off diagonals.
     *  3) sort the off diagonal column indices.
     *  4) sort the off diagonal matrix nozeros.
     */

    j = bindx[0];
    AZ_convert_ptrs_to_values(bindx, N_update);

    AZ_sortqlists((char *) &(bindx[N_update + 1]), bindx, update_index,
                  end - N_update - 1, sizeof(int), N_update);
    AZ_sortqlists((char *) bindx, 0, update_index, N_update,
                  sizeof(int), N_update);
    AZ_convert_values_to_ptrs(bindx, N_update, j);
  }
  else {
    indx_length = bnptr[N_update];
    AZ_convert_ptrs_to_values(indx, indx_length);
    AZ_convert_ptrs_to_values(bnptr, N_update);
    AZ_convert_ptrs_to_values(rnptr, N_update);

    /* move indx */

    AZ_sortqlists((char *) indx, bnptr, update_index, indx_length, sizeof(int),
                  N_update);

    /* move bindx */

    AZ_sortqlists((char *) bindx, bnptr, update_index, indx_length, sizeof(int),
                  N_update);

    /* move bnptr */

    AZ_sortqlists((char *) bnptr, 0, update_index, N_update, sizeof(int),
                  N_update);

    /* move rnptr */

    AZ_sortqlists((char *) rnptr, 0, update_index, N_update, sizeof(int),
                  N_update);

    AZ_convert_values_to_ptrs(rnptr, N_update, 0);
    AZ_convert_values_to_ptrs(bnptr, N_update, 0);
    AZ_convert_values_to_ptrs(indx, indx_length, 0);
  }

} /* AZ_reorder_matrix_mod */


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void AZ_transform_mod(int proc_config[], int *external[], int bindx[], 
                      int update[], int *update_index[], int *extern_index[],
                      int *data_org[], int N_update, int indx[], int bnptr[],
                      int rnptr[], int *cnptr[], int mat_type)

/*******************************************************************************

  Convert a global distributed matrix to a parallel local distributed matrix.
  This includes the following steps:
      1) reorder matrix rows so that all the rows corresponding to internal
         points preceed all the rows corresponding to border points.
      2) replace global indicies by local indicies.
      3) make a list of the external unknowns and store them in external[].
      4) make a list of processors which update each external unknown and store
         this list in extern_proc where extern_proc[i] is the processor that
         updates external[i].
      5) make 2 arrays (update_index[], extern_index[]) which define the mapping
         between global and local indicies. In particular, global index
         update[i] corresponds to the locally numbered variable update_index[i]
         and global index external[i] corresponds to the locally numbered
         variable extern_index[i].
      6) Initialize all the quanities in data_org[] to their appropriate values
         (so that communication can properly occur).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======
  Modified:        Edwin van der Weide, VKI
                   val does not occur anymore in the subroutine

  Return code:     void
  ============

  Parameter list:
  ===============

  proc_config:     Processor information corresponding to:
                      proc_config[AZ_node] = name of this processor
                      proc_config[AZ_N_procs] = # of processors used

  external:        On output, list of external blocks

  update:          On input, list of pts to be updated on this node

  update_index,    On output, ordering of update and external locally on this
  extern_index:    processor. For example  'update_index[i]' gives the index
                   location of the block which has the global index 'update[i]'.

  data_org:        On output, indicates how the data is set out on this node.
                   For example, data_org[] contains information on how many
                   unknowns are internal, external and border unknowns as well
                   as which points need to be communicated. See Aztec User's
                   guide for more details.

  N_update         Number of points to be updated on this node.

  bindx            On input, global distributed matrix (MSR or VBR) arrays
  indx, bnptr,     holding matrix values. On output, local renumbered matrix
  rnptr, cnptr:    (DMSR or DMVBR).

  mat_type:        Type of matrix (AZ_MSR_MATRIX or AZ_VBR_MATRIX).

*******************************************************************************/

{
  int        i, ii, j;
  static int mat_name = 1;

  int         N_extern;   /* Number of pts needed by this processor for
                             matrix-vector multiply but not updated by this
                             processor.  */
  int         N_internal, /* Number of pts which can be updated without
                             communication */
              N_border;   /* Number of pts to be updated requiring communication
                           */
  int        *extern_proc;
  int        *tcnptr = NULL;

  /*
   * Compute the external points and change the global indices to
   * local indices. That is,
   *   On input:                        On output:
   *      bindx[k] = update[j]      ==>   bindx[k] = j
   *      bindx[k] = external[j]    ==>   bindx[k] = j + N_update
   */

  AZ_find_local_indices(N_update, bindx, update, external, &N_extern, mat_type,
                        bnptr);

  /* compute the cnptr array for VBR matrices */

  if (mat_type == AZ_VBR_MATRIX) {
    if (!AZ_using_fortran)
      *cnptr = (int *) calloc(N_update + N_extern + 1, sizeof(int));

    tcnptr = *cnptr;

    for (i = 0; i < N_update; i++) tcnptr[i] = rnptr[i+1] - rnptr[i];

    for (i = 0; i < N_update; i++) {
      for (j = bnptr[i]; j < bnptr[i+1]; j++) {
        ii = bindx[j];

        if ((ii >= N_update) && ( tcnptr[ii] == 0)) {
          tcnptr[ii] = (indx[j+1]-indx[j]) / (rnptr[i+1]-rnptr[i]);
        }
      }
    }

    AZ_convert_values_to_ptrs(tcnptr, N_update + N_extern, 0);
  }

  /*
   * Read or compute (and sort) the processor numbers of the processors which
   * update the external points.
   */

  i                = AZ_using_fortran;
  AZ_using_fortran = AZ_FALSE;

  AZ_find_procs_for_externs(N_update, update, *external, N_extern, proc_config,
                            &extern_proc);
  AZ_using_fortran = i;

  /*
   * Determine a new ordering for the points:
   *    a) lowest numbers for internal points,
   *    b) next lowest numbers for border points
   *    c) highest nubers for the external points
   *       NOTE: external points updated by the same processor are consecutively
   *             ordered.
   */

  if (!AZ_using_fortran) {
    *update_index = (int *) calloc(N_update + 1, sizeof(int));
    *extern_index = (int *) calloc(N_extern + 1, sizeof(int));
  }

  if (*extern_index == NULL)  {
    (void) fprintf(stderr,
                   "Error: Not enough space in main() for extern_index[]\n");
    exit(1);
  }

  AZ_order_ele(*update_index, *extern_index, &N_internal, &N_border, N_update,
               bnptr, bindx, extern_proc, N_extern, AZ_ALL, mat_type);

  /*
   * Permute the matrix using the new ordering.  IMPORTANT: This routine assumes
   * that update_index[] contains 2 sequencies that are ordered but
   * intertwined. See AZ_reorder_matrix().
   */

  AZ_reorder_matrix_mod(N_update, bindx, *update_index, *extern_index,
                        indx, rnptr, bnptr, N_extern, tcnptr, AZ_ALL,mat_type);

  /*
   * Initialize 'data_org' so that local information can be exchanged to update
   * the external points.
   */

  AZ_set_message_info(N_extern, *extern_index, N_update, *external, extern_proc,
                      update, *update_index, proc_config, tcnptr, data_org,
                      mat_type);

  (*data_org)[AZ_name]       = mat_name;
  (*data_org)[AZ_N_int_blk]  = N_internal;
  (*data_org)[AZ_N_bord_blk] = N_border;
  (*data_org)[AZ_N_ext_blk]  = N_extern;

  if (mat_type == AZ_VBR_MATRIX) {
    (*data_org)[AZ_N_internal] = rnptr[N_internal];
    (*data_org)[AZ_N_external] = tcnptr[N_update + N_extern] - rnptr[N_update];
    (*data_org)[AZ_N_border]   = rnptr[N_update] - rnptr[N_internal];
  }

  else {
    (*data_org)[AZ_N_internal] = N_internal;
    (*data_org)[AZ_N_external] = N_extern;
    (*data_org)[AZ_N_border]   = N_border;
  }

  mat_name++;
  free(extern_proc);

} /* AZ_transform_mod */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_block_milu(double val[], int indx[], int bindx[],
                   int cpntr[], int bpntr[], int new_blks, int new_N,
                   int length, double buffer[], double x[], int options[],
                   int data_org[])

/*******************************************************************************

  Routine preconditions the vector 'x' using a modified incomplete block 
  factorization.
  If the factorization has not already been computed, this routine calculates
  the necessary block diagonal matrix.

  Note: At the moment only the option options[AZ_overlap] == AZ_none works properly.
        If desirable this could be changed in the future (not by me I think).

  Author:          Edwin van der Weide, VKI
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the entries of the matrix. The matrix is
                   stored block-row-by-block-row. Each block entry is dense and
                   stored by columns (VBR).

  indx,
  bindx,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  new_blks:        Number of blocks rows in expanded system corresponding to
                   domain decomposition.

  new_N:           Number of rows in expanded system corresponding to domain
                   decomposition.

  length:          Length of buffer.

  buffer:          Buffer holding extra rows received from neighbors.

  x:               On input, contains the current solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  int     i, j, k, l ;
  int     valsize;
  int    *ipvt;
  int     job;
  int     N, blks, nrow, ncol;
  double *Mtmp, *z, *x_help ;
  double  rcond, det[2], alpha, beta ;
  char    T = 'N' ;
  char    None[2] ;

  /* these arrays hold the block diagonal matrix */

  static int *indx2, *diag_block;

  /*  diag_block == int work array of size n points to each diagonal block */

  static double *val2;

  int            st;
  static int     previous_factors = -1;
  static int     max_row = 0 ;
  static int     first_call = 1 ;

  /* external functions */

  extern void AZ_lower_bmilu_solve(int new_blks, double *val, int *indx, 
                                   int *bindx, int *cpntr, int *bpntr,
                                   double *val2, int *indx2, double *x, 
                                   double *x_help) ;
  extern void AZ_upper_bmilu_solve(int new_blks, double *val, int *indx,
                                   int *bindx, int *cpntr, int *bpntr,
                                   double *val2, int *indx2, double *x, 
                                   double *x_help) ;

  /**************************** execution begins ******************************/

  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
  blks = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (N == 0) return;

  /***************************************************************************/
  /* At the moment only options[AZ_overlap] == AZ_none works. To save memory  */
  /* in the allocation new_N and new_blks are set to N and blks resp.         */
  /****************************************************************************/

  new_blks = blks ;
  new_N    = N ;

  if (options[AZ_pre_calc] == AZ_reuse) {

    /*
     * Previous block diagonal is being used. Find pointers corresponding to this
     * previous block diagonal.
     */

    indx2 = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                      data_org[AZ_name], "indx2", &st);
    if (st == AZ_NEW_ADDRESS) {
      (void) fprintf(stderr, "Error: options[AZ_pre_calc] == AZ_reuse and "
                     "previous factors\n       not found. "
                     "Check data_org[AZ_name].\n");
      exit(-1);
    }

    diag_block = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], "diag_block", &st);
    valsize = indx2[new_blks] + 1 ;
    val2    = (double *) AZ_manage_memory((valsize+1)*sizeof(double),
                                           AZ_ALLOC, data_org[AZ_name], "val2",
                                           &st);
  }

  else if (options[AZ_pre_calc] <= AZ_recalc) {

    /* compute a new block diagonal matrix */

    if (options[AZ_overlap] != AZ_none) {
      (void) fprintf(stderr, "Error: options[AZ_overlap] != AZ_none for "
                     "bmilu preconditioner.\n"
                     "Options not yet implemented.\n");
      exit(-1);
    }

    /* Manage the memory. */

    indx2 = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                      data_org[AZ_name], "indx2", &st);

    if( (!first_call) && (st==AZ_NEW_ADDRESS) ) {
      (void) fprintf(stderr, "This is not the first call and it cannot "
                     "find previous \n allocated memory for indx2.\n"
                     "Check data_org[AZ_name].\n");
      exit(-1);
    }

    diag_block = (int *) AZ_manage_memory((new_blks+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], "diag_block", &st);

    if( (!first_call) && (st==AZ_NEW_ADDRESS) ) {
      (void) fprintf(stderr, "This is not the first call and it cannot "
                     "find previous \n allocated memory for diag_block.\n"
                     "Check data_org[AZ_name].\n");
      exit(-1);
    }

    /* If this is the first call, figure out some values. */

    if( first_call ) {

      /* Set values to indx2 and determine what max_row is. */

      indx2[0] = 0 ;
      max_row  = 0 ;
      for(i=0; i<new_blks; i++)
      {
        indx2[i+1] = indx2[i] + (cpntr[i+1]-cpntr[i])*(cpntr[i+1]-cpntr[i]) ;
        max_row = max(max_row,(cpntr[i+1]-cpntr[i])) ;
      }
      max_row += 1 ;    /* Just to be sure it allocates something later on. */

      /* Determine the values of diag_block. */

      for(i=0; i<new_blks; i++)
      {
        diag_block[i] = -1 ;
        for(j=bpntr[i]; j<bpntr[i+1]; j++)
        {
          if(bindx[j]==i)
	  {
            diag_block[i] = j ;
            break ;
	  }
        }
        if(diag_block[i] == -1)
	{
          (void) fprintf(stderr, "No diagonal entry found for block %d.\n", i) ;
          (void) fprintf(stderr, "This implementation of the bmilu preconditioner "
                                 "cannot be used in this case.\n") ;
          exit(-1) ;
	}
      }
    /* NPW 16/7/00: forced always to think it's the first call */
    /*first_call = 0 ;*/
    }

    /* Manage the memory for val2. */

    valsize = indx2[new_blks] + 1 ;
    val2    = (double *) AZ_manage_memory((valsize+1)*sizeof(double),
                                           AZ_ALLOC, data_org[AZ_name], "val2",
                                           &st);

    /* Allocate some memory for help variables. */

    Mtmp = (double *) calloc(sizeof(double), max_row*max_row) ;
    z    = (double *) calloc(sizeof(double), max_row) ;
    ipvt = (int *) calloc(sizeof(int), max_row) ;

    /* Give the string None the value "N". */

    strcpy(None,"N") ;

    /* Copy the block diagonal of the full matrix into val2. */

    for(i=0; i<new_blks; i++)
      for(j=0; j<(indx2[i+1]-indx2[i]); j++)
        val2[indx2[i]+j] = val[indx[diag_block[i]]+j] ;

    /* Construct (the inverse of) the block diagonal of the preconditioner. */

    for(i=0; i<new_blks; i++)
    {
      /* Compute the inverse of this central block. */

      nrow = cpntr[i+1] - cpntr[i] ;
      dgeco_(val2+indx2[i], &nrow, &nrow, ipvt, &rcond, z) ;
      job = 01 ;
      dgedi_(val2+indx2[i], &nrow, &nrow, ipvt, det, z, &job) ;

      /* Loop over the entries for this point. Only modify succesors. */
      /* Also check whether the succesor is part of the update set.   */
      /* If not the entry is ignored (No overlap in preconditioner).  */

      for(j=bpntr[i]; j<bpntr[i+1]; j++)
      {
        l = bindx[j] ;
        if( (l>i) && (l<new_blks) )
	{
          /* Both a succesor and in update set. */

          ncol = cpntr[l+1] - cpntr[l] ;

          /* First matrix multiplication. */

          alpha = 1.0 ;
          beta  = 0.0 ;
          dgemm_(&T, &T, &nrow, &ncol, &nrow, &alpha, val2+indx2[i], &nrow,
                 val+indx[j], &nrow, &beta, Mtmp, &nrow, strlen(None), 
                 strlen(None)) ;

          /* Find the entry for the down jacobian. */

          for(k=bpntr[l]; k<bpntr[l+1]; k++)
            if(bindx[k] == i) break ;
          if(k == bpntr[l+1])
	  {
            (void) fprintf(stderr, "Down jacobian not found in"
                     "bmilu preconditioner.\n"
                     "This should not be possible.\n");
            exit(-1);
	  }

          /* Second matrix multiplication and substraction. */

          alpha = -1.0 ;
          beta  =  1.0 ;
          dgemm_(&T, &T, &ncol, &ncol, &nrow, &alpha, val+indx[k], &ncol,
                 Mtmp, &nrow, &beta, val2+indx2[l], &ncol, strlen(None),
                 strlen(None)) ;
	}
      }
    }

    /* Release some memory. */

    free((void *) Mtmp) ;
    free((void *) z) ;
    free((void *) ipvt) ;
  }

  else if (previous_factors != data_org[AZ_name] ) {
    (void) fprintf(stderr, "Warning: Using a previous factorization as a "
                   "preconditioner\neven though matrix "
                   "(data_org[AZ_name]) has changed\n");
  }

  previous_factors = data_org[AZ_name];

  if( !max_row )     /* Determine the value of max_row. */
  {
    for(i=0; i<new_blks; i++)
      max_row = max(max_row,(cpntr[i+1]-cpntr[i])) ;
    max_row += 1 ;  /* Just to be sure. */
  }

  /* Allocate memory for x_help. */

  x_help = (double *) calloc(sizeof(double), max_row) ;

  /* Solution of the system M.x_new = x_old. M is the preconditioner. */
  /* Solution is obtained by a lower and an upper sweep.              */

  AZ_lower_bmilu_solve(new_blks, val, indx, bindx, cpntr, bpntr,
                       val2, indx2, x, x_help) ;
  AZ_upper_bmilu_solve(new_blks, val, indx, bindx, cpntr, bpntr,
                       val2, indx2, x, x_help) ;

  free((void *) x_help) ;

}  /* AZ_block_milu */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_lower_bmilu_solve(int new_blks, double *val, int *indx, int *bindx,
                          int *cpntr, int *bpntr, double *val2, int *indx2,
                          double *x, double *x_help)

/*******************************************************************************

  Lower triangular solver for L.x_new = x_old with L partially stored in val and 
  partially in val2. As this is for bmilu preconditioner, val contains the
  original matrix and val2 the (inverse) of the block diagonal of the 
  preconditioner.


  Author:          Edwin van der Weide, VKI
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  new_blks:        Number of (block) rows in the matrix.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  val2:            Array containing the nonzero entries of the block diagonal
                   of the bmilu preconditioner.

  indx2:           pointer to numbers in val2.

  x:               On input:  Right hand side of the system.
                   On output: Solution of the system.

  x_help:          help array used in the sweep

*******************************************************************************/

{

  /* local variables */

  int     i, j, l ;
  int     nrow, ncol ;
  int     ione = 1;
  double *x_pntr ;
  double  plus_one = 1.0, minus_one = -1.0, my_zero = 0.0;
  char   *T = "N";

  /* external functions */

  extern void AZ_my_dgemv3(int m, int n, double *a, double *x, double *y);

  /**************************** execution begins ******************************/

  /* Loop over the block rows. */

  for(i=0; i<new_blks; i++) {

    /* Determine the row dimension and set the pointer. */

    nrow   = cpntr[i+1] - cpntr[i] ;
    x_pntr = x + cpntr[i] ;

    /* Loop over the lower blocks in the current block row.        */
    /* This automatically implies this block is in the update set. */

    for(j=bpntr[i]; j<bpntr[i+1]; j++) {

      l = bindx[j] ;

      /* Check whether l is smaller than i. If this is the case the value of  */
      /* the product between the "down" jacobian and x[l] must be substracted */
      /* from x[i].                                                           */

      if(l < i) {

	/* Determine the column dimension. */

        ncol = cpntr[l+1] - cpntr[l] ;

        /* Dense matrix vector multiplication. */

        if (nrow == 1 && ncol == 1)
          *x_pntr -= val[indx[j]]*x[cpntr[l]] ;
        else if(nrow < 10)
          AZ_dgemv3(nrow, ncol, val+indx[j], x+cpntr[l], x_pntr);
        else
          dgemv_(T, &nrow, &ncol, &minus_one, val+indx[j], &nrow, x+cpntr[l], 
                 &ione, &plus_one, x_pntr, &ione, strlen(T));
      }
    }

    /* Now the multiplication with the inverse of the diagonal needs to be done */
    /* to obtain the desired result. The array x_help is used to avoid trouble. */

    if (nrow == 1)
      x_help[0] = val2[indx2[i]]*x[cpntr[i]] ;
    else if(nrow < 10)
      AZ_my_dgemv3(nrow, nrow, val2+indx2[i], x_pntr, x_help);
    else
      dgemv_(T, &nrow, &nrow, &plus_one, val2+indx2[i], &nrow, x_pntr,
             &ione, &my_zero, x_help, &ione, strlen(T));

    for(j=0; j<nrow; j++)
      x[cpntr[i]+j] = x_help[j] ;
  }
}  /* AZ_lower_bmilu_solve */



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_upper_bmilu_solve(int new_blks, double *val, int *indx, int *bindx,
                          int *cpntr, int *bpntr, double *val2, int *indx2,
                          double *x, double *x_help)

/*******************************************************************************

  Upper triangular solver for P.U.x_new = x_old with U partially stored in val
  and partially in val2. The inverse of the block diagonal matrix P is stored 
  in val2. For the bmilu preconditioner the diagonal of U is the inverse of P.
  The resulting algorithm is a loop over the succesors in which 2 matrix-vector
  products are performed.

  Note: this routine is not meant for a standard upper triangular solve.


  Author:          Edwin van der Weide, VKI
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  new_blks:        Number of (block) rows in the matrix.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  val2:            Array containing the nonzero entries of the block diagonal
                   of the bmilu preconditioner.

  indx2:           pointer to numbers in val2.

  x:               On input:  Right hand side of the system.
                   On output: Solution of the system.

  x_help:          help array used in the sweep

*******************************************************************************/

{

  /* local variables */

  int     i, j, l ;
  int     nrow, ncol ;
  int     ione = 1;
  double *x_pntr ;
  double  plus_one = 1.0, minus_one = -1.0;
  char   *T = "N";

  /* external functions */

  extern void AZ_dgemv3_plus(int m, int n, double *a, double *x, double *y);

  /**************************** execution begins ******************************/

  /* Loop over the block rows, starting at highest number, */
  /* for it is a backward sweep.                           */

  for(i=(new_blks-1); i>=0; i--) {

    /* Determine the row dimension and set the pointer. */

    nrow   = cpntr[i+1] - cpntr[i] ;
    x_pntr = x + cpntr[i] ;

    /* Set the elements of x_help to zero. */

    for(j=0; j<nrow; j++)
      x_help[j] = 0.0 ;

    /* Loop over the upper blocks in the current block row.             */
    /* Check explicitly whether this block is in the update set or not. */

    for(j=bpntr[i]; j<bpntr[i+1]; j++) {

      l = bindx[j] ;

      /* Check whether l is bigger than i and in the update set. If this */
      /* is the case a vector, which is the result of a matrix vector    */
      /* multiplication, is added to x_help[i].                          */

      if((l>i) && (l<new_blks)) {

	/* Determine the column dimension. */

        ncol = cpntr[l+1] - cpntr[l] ;

        /* Dense matrix vector multiplications. */
        /* x_help is used as temporal storage.  */

        if (nrow == 1 && ncol == 1)
          x_help[0] += val[indx[j]]*x[cpntr[l]] ;
        else if(nrow < 10)
          AZ_dgemv3_plus(nrow, ncol, val+indx[j], x+cpntr[l], x_help);
        else
          dgemv_(T, &nrow, &ncol, &plus_one, val+indx[j], &nrow, x+cpntr[l],
                 &ione, &plus_one, x_help, &ione, strlen(T));
      }
    }

    /* Now the result of the central jacobian times x_help must be */
    /* substracted from the x itself.                              */

    if (nrow == 1 && ncol == 1)
      *x_pntr -= val2[indx2[i]]*x_help[0] ;
    else if(nrow < 10)
      AZ_dgemv3(nrow, nrow, val2+indx2[i], x_help, x_pntr);
    else
      dgemv_(T, &nrow, &nrow, &minus_one, val2+indx2[i], &nrow, x_help,
                 &ione, &plus_one, x_pntr, &ione, strlen(T));
  }
}  /* AZ_upper_bmilu_solve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_my_dgemv3(int m, int n, double *a, double *x, double *y)

/*******************************************************************************

  Perform the matrix-vector operation y = Ax where x and y are vectors and
  A is an m-by-n matrix.  This modeled on the blas 2 routine dgemv and used for
  small matrix sizes.  Also, there is a special case for m = n = 5, (3d fluids).

  Note: this routine is almost identical to AZ_dgemv3

  Author:          Edwin van der Weide, VKI
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of rows of the matrix a.

  n:               Number of columns of the matrix a.

  a:               The leading m-by-n part of the array a must contain the
                   matrix of coefficients.

  x:               The vector of length at least n to multiply by A.

  y:               The vector of length at least m which, on output, will
                   contain the result of the Ax operation summed into it.

*******************************************************************************/

{

  /* local variables */

  register int i = 0;

  /**************************** execution begins ******************************/

 switch (m) {

  case 5:

    *y = *(y+1) = *(y+2) = *(y+3) = *(y+4) = 0.0 ;
    while (i++ < n) {
      *y     += *a++ * *x;
      *(y+1) += *a++ * *x;
      *(y+2) += *a++ * *x;
      *(y+3) += *a++ * *x;
      *(y+4) += *a++ * *x;

      x++;

    }

    break;

  default:

     for(i=0; i<m; i++) *(y+i) = 0.0 ;
     i = 0 ;
     while (i++ < n) {
      register int j = 0;
      while (j++ < m)
        *y++ += *a++ * *x;
      y -= m; x++;
    }
  }
} /* AZ_my_dgemv3 */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


void AZ_dgemv3_plus(int m, int n, double *a, double *x, double *y)

/*******************************************************************************

  Perform the matrix-vector operation y = y + Ax where x and y are vectors and
  A is an m-by-n matrix.  This modeled on the blas 2 routine dgemv and used for
  small matrix sizes.  Also, there is a special case for m = n = 5, (3d fluids).

  Note: this routine is almost identical to AZ_dgemv3

  Author:          Edwin van der Weide, VKI
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of rows of the matrix a.

  n:               Number of columns of the matrix a.

  a:               The leading m-by-n part of the array a must contain the
                   matrix of coefficients.

  x:               The vector of length at least n to multiply by A.

  y:               The vector of length at least m which, on output, will
                   contain the result of the y + Ax operation summed into it.

*******************************************************************************/

{

  /* local variables */

  register int i = 0;

  /**************************** execution begins ******************************/

 switch (m) {
  case 5:
    while (i++ < n) {
      *y     += *a++ * *x;
      *(y+1) += *a++ * *x;
      *(y+2) += *a++ * *x;
      *(y+3) += *a++ * *x;
      *(y+4) += *a++ * *x;

      x++;

    }

    break;

  default:

     while (i++ < n) {
      register int j = 0;
      while (j++ < m)
        *y++ += *a++ * *x;
      y -= m; x++;
    }
  }

} /* AZ_dgemv3_plus */
