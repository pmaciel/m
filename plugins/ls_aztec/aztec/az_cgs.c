/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_cgs.c,v $
 *
 * $Author: tuminaro $
 *
 * $Date: 1996/07/30 20:44:56 $
 *
 * $Revision: 1.16 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_cgs.c,v 1.16 1996/07/30 20:44:56 tuminaro Exp $";
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
#include <float.h>
#include "az_aztec.h"

void AZ_pcgs(double val[], int indx[], int bindx[], int rpntr[], int cpntr[],
         int bpntr[], double b[], double x[], double weight[], int options[],
         double params[], int data_org[], int proc_config[], double status[])

/*******************************************************************************

  This routine uses Sonneveld's(1984, 1989) Conjugate Gradient Squared algorithm
  to solve the nonsymmetric matrix problem Ax = b.

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

  weight:          Vector of weights for convergence norm #4.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N, NN, converged, one = 1, iter = 1, r_avail = AZ_TRUE, j;
  int          precond_flag, print_freq, proc, brkdown_will_occur = AZ_FALSE;
  double       alpha, beta, nalpha, true_scaled_r;
  double       *u, *v, *r, *rtilda, *p, *q, *w, init_time;
  double       rhonm1 = 1.0, rhon, rec_residual, actual_residual = -1.0;
  double       sigma, scaled_r_norm, epsilon, brkdown_tol = DBL_EPSILON;

  /**************************** execution begins ******************************/

  /* pull needed values out of parameter arrays */

  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  precond_flag = options[AZ_precond];
  epsilon      = params[AZ_tol];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];

  /* allocate memory for required vectors */

  NN     = (N + data_org[AZ_N_external])*sizeof(double) + 1;
  /* +1: make sure everybody allocates something */

  w      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "w in cgs", &j);
  u      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "u in cgs", &j);
  p      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "p in cgs", &j);
  r      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "r in cgs", &j);
  q      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "q in cgs", &j);
  rtilda = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "rtilda in cgs",
                                       &j);
  v      = (double *) AZ_manage_memory(NN, AZ_ALLOC, AZ_SYS, "v in cgs", &j);

  AZ_compute_residual(val, indx, bindx, rpntr, cpntr, bpntr, b, x, r, data_org);

  /* q, p <- 0 */

  for (i = 0; i < N; i++) q[i] = p[i] = 0.0;

  /* set rtilda */

  if (options[AZ_aux_vec] == AZ_resid) dcopy_(&N, r, &one, rtilda, &one);
  else AZ_random_vector(rtilda, data_org, proc_config);

  /*
   * Compute a few global scalars:
   * 1) ||r|| corresponding to options[AZ_conv]
   * 2) scaled ||r|| corresponding to options[AZ_conv]
   * 3) rhon = <rtilda, r>
   */

  AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, r,
                            weight, &rec_residual, &scaled_r_norm, options,
                            data_org, proc_config,&r_avail,r,rtilda,&rhon,
                            AZ_FIRST_TIME);
  true_scaled_r = scaled_r_norm;

  if ((options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) && (proc == 0))
    (void) fprintf(stdout, "\t\titer:    0\t\tresidual = %e\n", scaled_r_norm);

  converged = scaled_r_norm < epsilon;

  for (iter = 1; iter <= options[AZ_max_iter] && !converged; iter++) {
    if (brkdown_will_occur) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config);

      AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                params, true_scaled_r, actual_residual, options,
                                proc_config);
      return;
    }

    if (fabs(rhon) < brkdown_tol) {   /* possible breakdown problem */
      if (AZ_breakdown_f(N, r, rtilda, rhon, proc_config))
        brkdown_will_occur = AZ_TRUE;
      else
        brkdown_tol = 0.1*fabs(rhon);
    }

    beta   = rhon / rhonm1;
    rhonm1 = rhon;

    /* u = r + beta*q            */
    /* p = u + beta*(q + beta*p) */
    /* w = M^-1 p                */
    /* v = A w                   */

    for (i = 0; i < N; i++) u[i] = r[i] + beta * q[i];
    daxpy_(&N, &beta, p, &one, q, &one);

    for (i = 0; i < N; i++) p[i] = u[i] + beta * q[i];
    dcopy_(&N, p, &one, w, &one);

    if (iter==1) init_time = AZ_second();

    if (precond_flag) AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, w,
                                      options, data_org, proc_config, params);

    if (iter==1) status[AZ_first_precond] = AZ_second() - init_time;

    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, w, v, 1, data_org);

    sigma = AZ_gdot(N, rtilda, v, proc_config);
    if (fabs(sigma) < brkdown_tol) { /* possible problem */

      if (AZ_breakdown_f(N, rtilda, v, sigma, proc_config)) {

        /* break down */

        AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v,
                               weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config);

        AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                  params, true_scaled_r, actual_residual,
                                  options, proc_config);
        return;
      }
      else brkdown_tol = 0.1 * fabs(sigma);
    }

    alpha  = rhon / sigma;
    nalpha = -alpha;

    /* q = u - alpha*v  */
    /* w = M^-1 (u + q) */
    /* x = x + alpha*w  */
    /* v = Aw           */
    /* r = r - alpha*v  */

    for (i = 0; i < N; i++) q[i] = u[i] - alpha * v[i];
    for (i = 0; i < N; i++) w[i] = u[i] + q[i];

    if (precond_flag) AZ_precondition(val, indx, bindx, rpntr, cpntr, bpntr, w,
                                      options, data_org, proc_config, params);

    daxpy_(&N, &alpha, w, &one, x, &one);
    AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, w, v, 1, data_org);
    daxpy_(&N, &nalpha, v, &one, r, &one);

    /*
     * compute a few global scalars:
     *     1) ||r||                corresponding to options[AZ_conv]
     *     2) scaled ||r||         corresponding to options[AZ_conv]
     *     3) rhon = <rtilda, r>
     */

    AZ_compute_global_scalars(val, indx, bindx, rpntr, cpntr, bpntr, x, b, r,
                              weight, &rec_residual, &scaled_r_norm, options,
                              data_org, proc_config, &r_avail, r, rtilda, &rhon,
                              AZ_NOT_FIRST);

    if ( (iter%print_freq == 0) && proc == 0 )
      (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                     scaled_r_norm);

    /* convergence tests */

    if (scaled_r_norm < epsilon) {
      AZ_scale_true_residual(val, indx, bindx, rpntr, cpntr, bpntr, x, b, v,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config);

      converged = true_scaled_r < params[AZ_tol];

      /*
       * Note: epsilon and params[AZ_tol] may not be equal due to a previous
       * call to AZ_get_new_eps().
       */

      if (!converged &&
          (AZ_get_new_eps(&epsilon, scaled_r_norm, true_scaled_r,
                          proc_config) == AZ_QUIT)) {

        /*
         * Computed residual has converged, actual residual has not converged,
         * AZ_get_new_eps() has decided that it is time to quit.
         */

        AZ_terminate_status_print(AZ_loss, iter, status, rec_residual, params,
                                  true_scaled_r, actual_residual, options,
                                  proc_config);
        return;
      }
    }
  }

  iter--;

  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings))
    (void) fprintf(stdout, "\t\titer: %4d\t\tresidual = %e\n", iter,
                   scaled_r_norm);

  /* check if we exceeded maximum number of iterations */

  if (converged) {
    i             = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* pcgs */
