
#include "common.h"


/*
 * This file contains routines to generate random numbers.
 *
 * This functions generates random numbers according to a
 * Gaussian distribution with mean m and standard deviation s.
 *
 * (c) Copyright 1994, Everett F. Carter Jr.
 * Permission is granted by the author to use this software
 * for any application provided this copyright notice
 * is preserved.
 */

double plas_RandomGaussian(float m, float s)
{
  double xx1, xx2, w, yy1;
  static double yy2;
  static int use_last = 0;

  if (use_last){
    yy1 = yy2;
    use_last = 0;
  } else{
    do {
      xx1 = 2.0 * plas_RandomDouble() - 1.0;
      xx2 = 2.0 * plas_RandomDouble() - 1.0;
      w = xx1 * xx1 + xx2 * xx2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    yy1 = xx1 * w;
    yy2 = xx2 * w;
    use_last = 1;
  }

  return( m + yy1 * s );
}

