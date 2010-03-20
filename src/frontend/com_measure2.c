/* New routines to evaluate the .measure cards.
   Entry point is function get_measure2(), called by fcn do_measure()
   from measure.c.
   Patches by Bill Swartz from 2009-05-18 and 2009-08-21 are included.

   $Id$
*/
#include <config.h>
#include <ngspice.h>
#include <memory.h>

#include <fteext.h>
#include <wordlist.h>

#include "vectors.h"
#include <math.h>
#include "dotcards.h"
#include "com_measure2.h"

typedef enum {
   MEASUREMENT_OK = 0,
   MEASUREMENT_FAILURE = 1
 } MEASURE_VAL_T ;

#define MEASURE_DEFAULT -1
#define MEASURE_LAST_TRANSITION  -2

typedef struct measure
{
  char *result;

  char *m_vec;      // name of the output variable which determines the beginning of the measurement
  char *m_vec2;     // second output variable to measure if applicable
  char *m_analysis;     // analysis type (tran, dc or ac)
  char m_vectype;       // type of vector m_vec (vm, vi, vr, vp, vdb)
  char m_vectype2;      // type of vector m_vec2 (vm, vi, vr, vp, vdb)
  int m_rise;       // count number of rise events
  int m_fall;       // count number of fall events
  int m_cross;      // count number of rise/fall aka cross  events
  double m_val;     // value of the m_ver at which the counter for crossing, rises or falls is incremented by one
  double m_td;      // amount of delay before the measurement should start
  double m_from;    // measure only in a time window - starting time of window
  double m_to;      // measurement window - ending time
  double m_at;      // measure at the specified time
  double m_measured;    // what we measured
  double m_measured_at; // what we measured at the given time

} MEASURE, *MEASUREPTR ;

typedef enum AnalysisType {
   AT_UNKNOWN, AT_DELAY, AT_TRIG,
   AT_FIND, AT_WHEN,
   AT_AVG, AT_MIN, AT_MAX, AT_RMS, AT_PP,
   AT_INTEG, AT_DERIV,
   AT_ERR, AT_ERR1, AT_ERR2, AT_ERR3, AT_MIN_AT, AT_MAX_AT
} ANALYSIS_TYPE_T ;

/** return precision (either 5 or value of environment variable NGSPICE_MEAS_PRECISION) */
int
measure_get_precision(void)
{
   char *env_ptr;
   int  precision = 5;

   if ( ( env_ptr = getenv("NGSPICE_MEAS_PRECISION") ) ) {
      precision = atoi(env_ptr);
   }

   return precision;
} /* end measure_get_precision() */

static void measure_errMessage(char *mName, char *mFunction, char *trigTarg, char *errMsg, int chk_only)
{

   if(!(chk_only)){
      printf("\tmeasure '%s'  failed\n", mName);
      printf("Error: measure  %s  %s(%s) :\n", mName, mFunction, trigTarg);
      printf("\t%s\n",errMsg);
   }
   return;
} /* end measure_errMessage() */


/* If you have a vector vm(out), extract 'm' to meas->m_vectype
   and v(out) to meas->m_vec (without 'm') */

static
void correct_vec(MEASUREPTR meas)
{
   char *vec, *vecfirst, newvec[BSIZE_SP];
   char *vec2, newvec2[BSIZE_SP];

   vec = meas->m_vec;
   /* return if not of type VM() etc */
   if ((*vec != 'v') || (!strstr(vec, "("))) return;
   
   if (*(++vec) != '(') {
      vecfirst = copy(meas->m_vec);
      vecfirst[1] = '\0';
      meas->m_vectype = *vec;
      sprintf(newvec, "%s%s", vecfirst, strstr(meas->m_vec, "("));
      tfree(meas->m_vec);
      tfree(vecfirst);
      meas->m_vec = copy(newvec);
   }

   vec2 = meas->m_vec2;
   if (vec2 && (*(++vec2) != '(')) {
      vecfirst = copy(meas->m_vec);
      vecfirst[1] = '\0';
      meas->m_vectype2 = *vec2;
      sprintf(newvec, "%s%s", vecfirst, strstr(meas->m_vec2, "("));
      tfree(meas->m_vec2);
      tfree(vecfirst);
      meas->m_vec2 = copy(newvec2);
   }
   return;
}

/* Returns a value from a complex vector *values, depending on meas->m_vectype */
static double get_value(
   MEASUREPTR meas,     /*in: pointer to mesurement structure */
   struct dvec *values, /*in: vector of complex values */
   int idx              /*in: index of vector value to be read out */
) {
   double ar, bi, tt;

   ar = values->v_compdata[idx].cx_real;
   bi = values->v_compdata[idx].cx_imag;

   if ((meas->m_vectype == 'm') || (meas->m_vectype == 'M'))
   {
      return sqrt(ar*ar + bi*bi); /* magnitude */
   }
   else if ((meas->m_vectype == 'r') || (meas->m_vectype == 'R'))
   {
      return ar;  /* real value */
   }
   else if ((meas->m_vectype == 'i') || (meas->m_vectype == 'I'))
   {
      return bi;  /* imaginary value */
   }
   else if ((meas->m_vectype == 'p') || (meas->m_vectype == 'P'))
   {
      return radtodeg(atan2(bi, ar)); /* phase (in degrees) */
   }
   else if ((meas->m_vectype == 'd') || (meas->m_vectype == 'D'))
   {
      tt = sqrt(ar*ar + bi*bi);  /* dB of magnitude */
      return 20.0 * log10(tt);
   }
   else
      return ar;  /* default: real value */
}

/* Returns interpolated value. If ac simulation, exploit vector type with complex data for y */
static double
measure_interpolate(
   struct dvec *xScale, /* in: vector of independent variables, if ac: complex vector,
                           but only real part used */
   struct dvec *values, /* in: vector of dependent variables, if ac: complex vector */
   int i,               /* in: index of first interpolation value */
   int j,               /* in: index of second interpolation value */
   double var_value,    /* in: variable, whose counterpart is sought by interpolation */
   char x_or_y ,        /* in: if 'x', then look for y, if 'y' then look for x */
   MEASUREPTR meas      /* pointer to measurement structure */
) {
   double slope;
   double yint;
   double result;

   if (cieq (meas->m_analysis,"ac")) {
      /* get values from complex y vector according to meas->m_vectype,
         x vector uses only real part of complex data (frequency).*/
      slope = (get_value(meas, values, j)  - get_value(meas, values, i)) /
         (xScale->v_compdata[j].cx_real - xScale->v_compdata[i].cx_real);
      yint  = get_value(meas, values, i) - slope*xScale->v_compdata[i].cx_real;
   }
   else {
      slope = (values->v_realdata[j] - values->v_realdata[i]) /
         (xScale->v_realdata[j] - xScale->v_realdata[i]);
      yint  = values->v_realdata[i] - slope*xScale->v_realdata[i];
   }

   if ( x_or_y == 'x' ) result = (var_value - yint)/slope;
   else                 result = slope*var_value + yint;

   return result;
} /* end measure_interpolate() */


/* -----------------------------------------------------------------
 * Function: Given an operation string returns back the measure type -
 * one of the enumerated type ANALSYS_TYPE_T.
 * ----------------------------------------------------------------- */
static ANALYSIS_TYPE_T measure_function_type( char *operation )
{
    char *mFunction ;           /* operation */
    ANALYSIS_TYPE_T mFunctionType ; /* type of requested function */

    mFunction = cp_unquote(operation);
    // Functions
    if (strcasecmp(mFunction,"DELAY")==0)
        mFunctionType = AT_DELAY;
    else if (strcasecmp(mFunction,"TRIG")==0)
        mFunctionType = AT_DELAY;
    else if (strcasecmp(mFunction,"TARG")==0)
        mFunctionType = AT_DELAY;
    else if (strcasecmp(mFunction,"FIND")==0)
        mFunctionType = AT_FIND;
    else if (strcasecmp(mFunction,"WHEN")==0)
        mFunctionType = AT_WHEN;
    else if (strcasecmp(mFunction,"AVG")==0)
        mFunctionType = AT_AVG;
    else if (strcasecmp(mFunction,"MIN")==0)
        mFunctionType = AT_MIN;
    else if (strcasecmp(mFunction,"MAX")==0)
        mFunctionType = AT_MAX;
    else if (strcasecmp(mFunction,"MIN_AT")==0)
        mFunctionType = AT_MIN_AT;
    else if (strcasecmp(mFunction,"MAX_AT")==0)
        mFunctionType = AT_MAX_AT;
    else if (strcasecmp(mFunction,"RMS")==0)
        mFunctionType = AT_RMS;
    else if (strcasecmp(mFunction,"PP")==0)
        mFunctionType = AT_PP;
    else if (strcasecmp(mFunction,"INTEG")==0)
        mFunctionType = AT_INTEG;
    else if (strcasecmp(mFunction,"DERIV")==0)
        mFunctionType = AT_DERIV;
    else if (strcasecmp(mFunction,"ERR")==0)
        mFunctionType = AT_ERR;
    else if (strcasecmp(mFunction,"ERR1")==0)
        mFunctionType = AT_ERR1;
    else if (strcasecmp(mFunction,"ERR2") == 0)
        mFunctionType = AT_ERR2;
    else if (strcasecmp(mFunction,"ERR3") == 0)
        mFunctionType = AT_ERR3;
    else
        mFunctionType = AT_UNKNOWN;

    return( mFunctionType) ;

} /* end measure_function_type() */

/* -----------------------------------------------------------------
 * Function: Parse the measurement line and extract any variables in
 * the statement and call com_save2 to instantiate the variable as a
 * measurement vector in the transient analysis.
 * ----------------------------------------------------------------- */
int measure_extract_variables( char *line )
{
  /* Various formats for measure statement:
   * .MEASURE {DC|AC|TRAN} result TRIG trig_variable VAL=val
   * + <TD=td> <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   * + <TRIG AT=time>
   * + TARG targ_variable VAL=val
   * + <TD=td> <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   * + <TRIG AT=time>
   *
   * .MEASURE {DC|AC|TRAN} result WHEN out_variable=val
   * + <TD=td> <FROM=val> <TO=val>
   * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   *
   * .MEASURE {DC|AC|TRAN} result WHEN out_variable=out_variable2
   * + <TD=td> <FROM=val> <TO=val>
   * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   *
   * .MEASURE {DC|AC|TRAN} result FIND out_variable WHEN out_variable2=val
   * + <TD=td> <FROM=val> <TO=val>
   * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   *
   * .MEASURE {DC|AC|TRAN} result FIND out_variable WHEN out_variable2=out_variable3
   * + <TD=td>
   * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   *
   * .MEASURE {DC|AC|TRAN} result FIND out_variable AT=val
   * + <FROM=val> <TO=val>
   *
   * .MEASURE {DC|AC|TRAN} result {AVG|MIN|MAX|MIN_AT|MAX_AT|PP|RMS} out_variable
   * + <TD=td> <FROM=val> <TO=val>
   *
   * .MEASURE {DC|AC|TRAN} result INTEG<RAL> out_variable
   * + <TD=td> <FROM=val> <TO=val>
   *
   * .MEASURE {DC|AC|TRAN} result DERIV<ATIVE> out_variable AT=val
   *
   * .MEASURE {DC|AC|TRAN} result DERIV<ATIVE> out_variable WHEN out_variable2=val
   * + <TD=td>
   * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   *
   * .MEASURE {DC|AC|TRAN} result DERIV<ATIVE> out_variable WHEN out_variable2=out_variable3
   * + <TD=td>
   * + <CROSS=# | CROSS=LAST> <RISE=#|RISE=LAST> <FALL=#|FALL=LAST>
   * ----------------------------------------------------------------- */

   int len ;                   /* length of string */
   int status ;                /* return status */
   char *item ;                /* parsing item */
   char *measure ;             /* measure keyword */
   char *analysis ;                /* analysis option */
   char *variable ;                /* variable to trace */
   wordlist *measure_var ;         /* wordlist of measurable */
   ANALYSIS_TYPE_T op ;            /* measure function type */

   status = TRUE;
   measure = gettok(&line);
   if(!(measure)){
      return(status) ;
   }
   analysis = gettok(&line);
   if(!(analysis)){
      return(status) ;
   }
   if( (strcasecmp(analysis,"DC")==0) ||
      (strcasecmp(analysis,"AC")==0) ||
      (strcasecmp(analysis,"TRAN")==0) ){
      analysis = copy(analysis) ;
   } else {
   /* sometimes operation is optional - for now just pick trans */
      analysis = copy("TRAN") ;
   }
   do {
      item = gettok(&line) ;
      if( item ){
         op = measure_function_type(item) ;
         if( op != AT_UNKNOWN ){
            /* We have a variable/complex variable coming next */
            variable = gettok(&line) ;
            if( variable ){
               len = strlen(item) ;
               if( item[len-1] == '=' ){
               } else {
                  measure_var = gettoks(variable) ;
                  com_save2(measure_var, analysis);
                  status = FALSE;
               }
            }
         }
      }
   } while(line && *line) ;

   return( status ) ;
} /* end measure_extract_variables() */

/* -----------------------------------------------------------------
 * Function: process a WHEN measurement statement which has been
 * parsed into a measurement structure.
 * ----------------------------------------------------------------- */
static void com_measure_when(
    MEASUREPTR meas     /* in : parsed measurement structure */
) {

   int i, first;
   int riseCnt = 0;
   int fallCnt = 0;
   int crossCnt = 0;
   int section = -1;
   int measurement_pending;
   int init_measured_value;
   double value, prevValue;
   double timeValue, prevTimeValue;

   enum ValSide { S_ABOVE_VAL, S_BELOW_VAL };
   enum ValEdge { E_RISING, E_FALLING };

   struct dvec *d, *dTime;

   d = vec_get(meas->m_vec);
   dTime = plot_cur->pl_scale;

   if (d == NULL) {
      fprintf(cp_err, "Error: no such vector as %s.\n", meas->m_vec);
      return;
   }

   if (dTime == NULL) {
       fprintf(cp_err, "Error: no such vector as time.\n");
       return;
   }

   prevValue =0;
   prevTimeValue =0;
   first =0;
   measurement_pending=0;
   init_measured_value=1;

   for (i=0; i < d->v_length; i++) {

//      value = d->v_realdata[i];
//      timeValue = dTime->v_realdata[i];

      if (cieq (meas->m_analysis,"ac")) {
         value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         timeValue = dTime->v_compdata[i].cx_real;
      }
      else if (cieq (meas->m_analysis,"sp")) {
         if (d->v_compdata)
            value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         else
            value = d->v_realdata[i];
         timeValue = dTime->v_realdata[i];
      }
      else {
         value = d->v_realdata[i];
         timeValue = dTime->v_realdata[i];
      }

      if (timeValue < meas->m_td)
         continue;

      if (first == 1) {
         // initialise
         crossCnt =0;
         if (value < meas->m_val) {
            section = S_BELOW_VAL;
            if ( (prevValue <= meas->m_val) && (value >= meas->m_val) ) {
               fallCnt =1;
               crossCnt =1;
            }

         } else {
            section = S_ABOVE_VAL;
            if ( (prevValue <= meas->m_val) && (value >= meas->m_val) ) {
               riseCnt =1;
               crossCnt =1;
            }
         }
         fflush( stdout ) ;
      }

      if (first > 1) {
         if ( (section == S_BELOW_VAL) && (value >= meas->m_val) ) {
            section = S_ABOVE_VAL;
            crossCnt++;
            riseCnt++;
            if( meas->m_fall != MEASURE_LAST_TRANSITION ){
               /* we can measure rise/cross transition if the user
                * has not requested a last fall transition */
               measurement_pending=1;
            }

         } else if ( (section == S_ABOVE_VAL) && (value <= meas->m_val) ) {
            section = S_BELOW_VAL;
            crossCnt++;
            fallCnt++;
            if( meas->m_rise != MEASURE_LAST_TRANSITION ){
               /* we can measure fall/cross transition if the user
                * has not requested a last rise transition */
               measurement_pending=1;
            }
         }

         if  ((crossCnt == meas->m_cross) || (riseCnt == meas->m_rise) || (fallCnt == meas->m_fall)) {
             /* user requested an exact match of cross, rise, or fall
              * exit when we meet condition */
            meas->m_measured = prevTimeValue + (meas->m_val - prevValue) * (timeValue - prevTimeValue) / (value - prevValue);
            return;
         }
         if  ( measurement_pending ){
            if( (meas->m_cross == MEASURE_DEFAULT) && (meas->m_rise == MEASURE_DEFAULT) && (meas->m_fall == MEASURE_DEFAULT) ){
               /* user didn't request any option, return the first possible case */
               meas->m_measured = prevTimeValue + (meas->m_val - prevValue) * (timeValue - prevTimeValue) / (value - prevValue);
               return;
            } else if( (meas->m_cross == MEASURE_LAST_TRANSITION) || (meas->m_rise == MEASURE_LAST_TRANSITION) || (meas->m_fall == MEASURE_LAST_TRANSITION) ){
               meas->m_measured = prevTimeValue + (meas->m_val - prevValue) * (timeValue - prevTimeValue) / (value - prevValue);
               /* no return - look for last */
               init_measured_value=0;
            }
            measurement_pending=0;
         }
      }
      first ++;

      prevValue = value;
      prevTimeValue = timeValue;
   }

   if ( init_measured_value ){
      meas->m_measured = 0.0e0;
   }
   return;
}

/* -----------------------------------------------------------------
 * Function: process an AT measurement statement which has been
 * parsed into a measurement structure.  We make sure to interpolate
 * the value when appropriate.
 * ----------------------------------------------------------------- */
static void measure_at(
    MEASUREPTR meas,        /* in : parsed "at" data */
    double at           /* in: time to perform measurement */
) {

   int i;
   double value, pvalue, svalue, psvalue;
   struct dvec *d, *dScale;

   psvalue = pvalue = 0;
   d = vec_get(meas->m_vec);
   dScale = plot_cur->pl_scale;

   if (d == NULL) {
      fprintf(cp_err, "Error: no such vector as %s.\n", meas->m_vec);
      return;
   }

   if (dScale == NULL) {
      fprintf(cp_err, "Error: no such vector time, frequency or dc.\n");
      return;
   }

   for (i=0; i < d->v_length; i++) {
      if (cieq (meas->m_analysis,"ac")) {
         value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         svalue = dScale->v_compdata[i].cx_real;
      }
      else if (cieq (meas->m_analysis,"sp")) {
         if (d->v_compdata)
            value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         else
            value = d->v_realdata[i];      	
         svalue = dScale->v_realdata[i];
      }
      else {
         value = d->v_realdata[i];
         svalue = dScale->v_realdata[i];
      }

      if ( (i > 0) && (psvalue <= at) && (svalue >= at) ) {
         meas->m_measured = pvalue + (at - psvalue) * (value - pvalue) / (svalue - psvalue);
        //  meas->m_measured = value;
         return;
      }

      psvalue = svalue;
      pvalue = value;
   }

   meas->m_measured = 0.0e0;
   return;
}


/* -----------------------------------------------------------------
 * Function: process an MIN, MAX, or AVG statement which has been
 * parsed into a measurement structure.  We should make sure to interpolate
 * the value here when we have m_from and m_to constraints * so this
 * function is slightly wrong.   Need to fix in future rev.
 * ----------------------------------------------------------------- */
static void measure_minMaxAvg(
   MEASUREPTR meas,                /* in : parsed measurement data request */
   ANALYSIS_TYPE_T mFunctionType   /* in: one of AT_AVG, AT_MIN, AT_MAX, AT_MIN_AT, AT_MAX_AT */
) {

   int i, avgCnt;
   struct dvec *d, *dScale;
   double value, svalue, mValue, mValueAt;
   int first;

   mValue =0;
   mValueAt = svalue =0;
   meas->m_measured = 0.0e0;
   meas->m_measured_at = 0.0e0;
   first =0;
   avgCnt =0;

   d = vec_get(meas->m_vec);
   if (d == NULL) {
      fprintf(cp_err, "Error: no such vector as %s.\n", meas->m_vec);
      return;
   }

   if (cieq (meas->m_analysis,"ac") || cieq (meas->m_analysis,"sp"))
      dScale = vec_get("frequency");
   else if (cieq (meas->m_analysis,"tran"))
      dScale = vec_get("time");
   else if (cieq (meas->m_analysis,"dc"))
      dScale = vec_get("v-sweep");
   else {/* error */
      fprintf(cp_err, "Error: no such analysis type as %s.\n", meas->m_analysis);
      return;
   }
   if (dScale == NULL) {
      fprintf(cp_err, "Error: no such vector as time, frquency or dc.\n");
      return;
   }

   for (i=0; i < d->v_length; i++) {
      if (cieq (meas->m_analysis,"ac")) {
         value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         svalue = dScale->v_compdata[i].cx_real;
      }
      else if (cieq (meas->m_analysis,"sp")) {
         if (d->v_compdata)
            value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         else
            value = d->v_realdata[i];      	
         svalue = dScale->v_realdata[i];
      }
      else {
         value = d->v_realdata[i];
         svalue = dScale->v_realdata[i];
      }

      if (svalue < meas->m_from)
         continue;

      if ((meas->m_to != 0.0e0) && (svalue > meas->m_to) )
         break;

      if (first ==0) {
         mValue = value;
         mValueAt = svalue;
         first =1;

      } else  {
         switch (mFunctionType) {
         case AT_MIN:
         case AT_MIN_AT: {
            if (value <= mValue) {
               mValue = value;
               mValueAt = svalue;
            }
            break;
         }
         case AT_MAX_AT:
         case AT_MAX: {
            if (value >= mValue) {
               mValue = value;
               mValueAt = svalue;
            }
            break;
         }
         case AT_AVG: {
            mValue = mValue + value;
            avgCnt ++;
            break;
         }
         default :
            fprintf(cp_err, "Error: improper min/max/avg call.\n");
         }

      }
   }

   switch (mFunctionType)
   {
      case AT_AVG: {
         meas->m_measured = (mValue / avgCnt);
         meas->m_measured_at = svalue;
         break;
      }
      case AT_MIN:
      case AT_MAX:
      case AT_MIN_AT:
      case AT_MAX_AT: {
         meas->m_measured = mValue;
         meas->m_measured_at = mValueAt;
         break;
      }
      default :
         fprintf(cp_err, "Error: improper min/max/avg call.\n");
   }
   return;
}

/* -----------------------------------------------------------------
 * Function: process an RMS or INTEG statement which has been
 * parsed into a measurement structure.  Here we do interpolate
 * the starting and stopping time window so the answer is correct.
 * ----------------------------------------------------------------- */
static void measure_rms_integral(
   MEASUREPTR meas,                    /* in : parsed measurement data request */
   ANALYSIS_TYPE_T mFunctionType   /* in: one of AT_RMS, or AT_INTEG */
) {
        int i;              /* counter */
   int xy_size ;           /* # of temp array elements */
   struct dvec *d, *xScale;        /* value and indpendent (x-axis) vectors */
   double value, xvalue;       /* current value and independent value */
   double *x ;         /* temp x array */
   double *y ;         /* temp y array */
   double toVal ;          /* to time value */
   double *width ;         /* temp width array */
   double sum1 ;           /* first sum */
   double sum2 ;           /* second sum */
   double sum3 ;           /* third sum */
   int first;

   xvalue =0;
   meas->m_measured = 0.0e0;
   meas->m_measured_at = 0.0e0;
   first =0;

   d = vec_get(meas->m_vec);
   if (d == NULL) {
      fprintf(cp_err, "Error: no such vector as %s.\n", meas->m_vec);
      return;
   }

   if (cieq (meas->m_analysis,"ac") || cieq (meas->m_analysis,"sp"))
      xScale = vec_get("frequency");
   else if (cieq (meas->m_analysis,"tran"))
      xScale = vec_get("time");
   else if (cieq (meas->m_analysis,"dc"))
      xScale = vec_get("v-sweep");
   else {/* error */
      fprintf(cp_err, "Error: no such analysis type as %s.\n", meas->m_analysis);
      return;
   }

   if (xScale == NULL) {
      fprintf(cp_err, "Error: no such vector as time.\n");
      return;
   }

   /* Allocate buffers for calculation. */
   x     = (double *) tmalloc(xScale->v_length * sizeof(double));
   y     = (double *) tmalloc(xScale->v_length * sizeof(double));
   width = (double *) tmalloc((xScale->v_length + 1) * sizeof(double));

   xy_size = 0 ;
   toVal = -1 ;
   /* create new set of values over interval [from, to] -- interpolate if necessary */
   for (i=0; i < d->v_length; i++) {
      if (cieq (meas->m_analysis,"ac")) {
         value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         xvalue = xScale->v_compdata[i].cx_real;
      }
      else {
         value = d->v_realdata[i];
         xvalue = xScale->v_realdata[i];
      }

      if (xvalue < meas->m_from)
         continue;

      if ((meas->m_to != 0.0e0) && (xvalue > meas->m_to) ){
         // interpolate ending value if necessary.
         if (!(AlmostEqualUlps( xvalue, meas->m_to, 100))){
            value = measure_interpolate( xScale, d, i-1, i, meas->m_to, 'y', meas );
            xvalue = meas->m_to ;
         }
         x[xy_size] = xvalue ;
         if (mFunctionType == AT_RMS)
            y[xy_size++] = value * value ;
         else
            y[xy_size++] = value ;
         toVal = xvalue ;
         break;
      }

      if (first == 0) {
         if( meas->m_from != 0.0e0 && (i > 0) ){
            // interpolate starting value.
            if (!(AlmostEqualUlps( xvalue, meas->m_from, 100))){
               value = measure_interpolate( xScale, d, i-1, i, meas->m_from, 'y' , meas);
               xvalue = meas->m_from ;
            }
         }
         meas->m_measured_at = xvalue ;
         first = 1;
      }
      x[xy_size] = xvalue ;
      if (mFunctionType == AT_RMS)
         y[xy_size++] = value * value ;
      else
         y[xy_size++] = value ;
   }

   // evaluate segment width
   for ( i = 0; i < xy_size-1; i++ ) width[i] = x[i+1] - x[i] ;
   width[i++] = 0;
   width[i++] = 0;

   // Compute Integral (area under curve)
   i = 0;
   sum1 = sum2 = sum3 = 0.0 ;
   while ( i < xy_size-1 ) {
      // Simpson's 3/8 Rule
      if ( AlmostEqualUlps( width[i], width[i+1], 100 ) &&
           AlmostEqualUlps( width[i], width[i+2], 100 ) ) {
         sum1 += 3*width[i] * (y[i] + 3*(y[i+1] + y[i+2]) + y[i+3]) / 8.0;
         i += 3;
      }
      // Simpson's 1/3 Rule
      else if ( AlmostEqualUlps( width[i], width[i+1], 100 ) ) {
         sum2 += width[i] * (y[i] + 4*y[i+1] + y[i+2]) / 3.0 ;
         i += 2;
      }
      // Trapezoidal Rule
      else if ( !AlmostEqualUlps( width[i], width[i+1], 100 ) ) {
         sum3 += width[i] * (y[i] + y[i+1]) / 2;
         i++;
      }
   }

   /* Now set the measurement values if not set */
   if( toVal < 0.0 ){
      if (cieq (meas->m_analysis,"ac")) {
         value = get_value(meas, d, i); //d->v_compdata[i].cx_real;
         xvalue = xScale->v_compdata[i].cx_real;
         toVal = xScale->v_compdata[d->v_length-1].cx_real;
      }
      else {
         toVal = xScale->v_realdata[d->v_length-1];
      }


   }
   meas->m_from = meas->m_measured_at ;
   meas->m_to = toVal ;

   if (mFunctionType == AT_RMS) {
      meas->m_measured = (sum1 + sum2 + sum3)/ (toVal - meas->m_measured_at) ;
      meas->m_measured = sqrt(meas->m_measured);

   } else {
      meas->m_measured = ( sum1 + sum2 + sum3 );
   }
   txfree(x); txfree(y); txfree(width);

} /* end measure_rms_integral() */

/* -----------------------------------------------------------------
 * Function: Wrapper function to process a RMS measurement.
 * ----------------------------------------------------------------- */
static void measure_rms(
    MEASUREPTR meas                 /* in : parsed measurement data request */
) {
   // RMS (root mean squared):
   // Calculates the square root of the area under the 'out_var2' curve
   //  divided be the period of interest
   measure_rms_integral(meas,AT_RMS) ;
   return;
}

/* -----------------------------------------------------------------
 * Function: Wrapper function to process a integration measurement.
 * ----------------------------------------------------------------- */
static void measure_integ(
   MEASUREPTR meas                 /* in : parsed measurement data request */
)  {
   // INTEGRAL INTEG
   measure_rms_integral(meas,AT_INTEG) ;
   return;
}

/* still some more work to do.... */
void measure_deriv( ) {
    // DERIVATIVE DERIV
   return;
}

// ERR Equations
void measure_ERR( ) {
   return;
}

void measure_ERR1( ) {
   return;
}

void measure_ERR2( ) {
   return;
}

void measure_ERR3( ) {
   return;
}

void com_dotmeasure( ) {
/* simulation info */
//      printf("*%s\n", plot_cur->pl_title);
//      printf("\t %s, %s\n", plot_cur->pl_name, plot_cur->pl_date); // missing temp

   return;
}

/* -----------------------------------------------------------------
 * Function: Given a measurement variable name, see if the analysis
 * has generated a measure vector for it.  Returns TRUE if it exists
 * or varname is NULL,  Return FALSE otherwise
 * ----------------------------------------------------------------- */
static int measure_valid_vector(
   char *varname         /* in: requested variable name */
) {

   struct dvec *d;         /* measurement vector */

   if(varname == NULL)
      return TRUE;
   d = vec_get(varname);
   if (d == NULL)
      return FALSE;

   return TRUE;
}

/* -----------------------------------------------------------------
 * Function: Given a wordlist and measurement structure, parse the
 * standard parameters such as RISE, FALL, VAL, TD, FROM, TO, etc.
 * in a measurement statement.   We also check the appropriate
 * variables found in the measurement statement.
 * ----------------------------------------------------------------- */
static int measure_parse_stdParams (
    MEASUREPTR meas,        /* in : measurement structure */
    wordlist *wl,       /* in : word list to parse */
    wordlist *wlBreak,      /* out: where we stopped parsing */
    char *errbuf        /* in/out: buffer where we write error messages */
) {

   int pCnt;
   char *p, *pName, *pValue;
   double *engVal, engVal1;

   pCnt = 0;
   while (wl != wlBreak) {
      p = wl->wl_word;
      pName = strtok(p, "=");
      pValue = strtok(NULL, "=");

      if (pValue == NULL) {
         if( strcasecmp(pName,"LAST")==0) {
            meas->m_cross = MEASURE_LAST_TRANSITION;
            meas->m_rise = -1;
            meas->m_fall = -1;
            pCnt ++;
            wl = wl->wl_next;
            continue ;
         } else {
            sprintf(errbuf,"bad syntax of ??\n");
            return 0;
         }
      }

      if( strcasecmp(pValue,"LAST")==0) {
         engVal1 = MEASURE_LAST_TRANSITION;
      } else {
         if (!(engVal = ft_numparse(&pValue, FALSE))) {
            sprintf(errbuf,"bad syntax of ??\n");
            return 0;
         }
         engVal1 = *engVal;  // What is this ??
      }

      if(strcasecmp(pName,"RISE")==0) {
         meas->m_rise = (int)engVal1;
         meas->m_fall = -1;
         meas->m_cross = -1;
      } else if(strcasecmp(pName,"FALL")==0) {
         meas->m_fall = (int)engVal1;
         meas->m_rise = -1;
         meas->m_cross = -1;
      } else if(strcasecmp(pName,"CROSS")==0) {
         meas->m_cross = (int)engVal1;
         meas->m_rise = -1;
         meas->m_fall = -1;
      } else if(strcasecmp(pName,"VAL")==0) {
         meas->m_val = engVal1;
      } else if(strcasecmp(pName,"TD")==0) {
         meas->m_td = engVal1;
      } else if(strcasecmp(pName,"FROM")==0) {
         meas->m_from = engVal1;
      } else if(strcasecmp(pName,"TO")==0) {
         meas->m_to = engVal1;
      } else if(strcasecmp(pName,"AT")==0) {
         meas->m_at = engVal1;
      } else {
         sprintf(errbuf,"no such parameter as '%s'\n",pName);
         return 0;
      }

      pCnt ++;
      wl = wl->wl_next;
   }

   if (pCnt == 0) {
      sprintf(errbuf,"bad syntax of ??\n");
      return 0;
   }

   // valid vector
   if (measure_valid_vector(meas->m_vec)==0) {
      sprintf(errbuf,"no such vector as '%s'\n", meas->m_vec);
      return 0;
   }

   // valid vector2
   if (meas->m_vec2 != NULL) {
      if (measure_valid_vector(meas->m_vec2)==0) {
         sprintf(errbuf,"no such vector as '%s'\n", meas->m_vec2);
         return 0;
      }
   }

   return 1;
}

/* -----------------------------------------------------------------
 * Function: Given a wordlist and measurement structure, parse a
 * FIND measurement statement.   Most of the work is done by calling
 * measure_parse_stdParams.
 * ----------------------------------------------------------------- */
static int measure_parse_find (
    MEASUREPTR meas,        /* in : measurement structure */
    wordlist *wl,       /* in : word list to parse */
    wordlist *wlBreak,      /* out: where we stopped parsing */
    char *errbuf        /* in/out: buffer where we write error messages */
) {

   int pCnt;
   char *p, *pName, *pVal;
   double *engVal, engVal1;

   meas->m_vec = NULL;
   meas->m_vec2 = NULL;
   meas->m_val = -1;
   meas->m_cross = -1;
   meas->m_fall = -1;
   meas->m_rise = -1;
   meas->m_td = 0;
   meas->m_from = 0.0e0;
   meas->m_to = 0.0e0;
   meas->m_at = -1;

   pCnt =0;
   while(wl != wlBreak) {
      p = wl->wl_word;

      if (pCnt == 0 ) {
         meas->m_vec= cp_unquote(wl->wl_word);
         /* correct for vectors like vm, vp etc. */
         if (cieq("ac", meas->m_analysis) || cieq("sp", meas->m_analysis))
            correct_vec(meas);
      } else if (pCnt == 1) {
         pName = strtok(p, "=");
         pVal = strtok(NULL, "=");

         if (pVal == NULL) {
            sprintf(errbuf,"bad syntax of WHEN\n");
            return 0;
         }

         if (strcasecmp(pName,"AT")==0) {
            if (!(engVal = ft_numparse(&pVal, FALSE))) {
               sprintf(errbuf,"bad syntax of WHEN\n");
               return 0;
            }

            engVal1 = *engVal;

            meas->m_at = engVal1;

         } else {
             sprintf(errbuf,"bad syntax of WHEN\n");
             return 0;
         }
      } else {
         if (measure_parse_stdParams(meas, wl, NULL, errbuf) == 0)
            return 0;
      }

      wl = wl->wl_next;
      pCnt ++;
   }

   return 1;
}

/* -----------------------------------------------------------------
 * Function: Given a wordlist and measurement structure, parse a
 * WHEN measurement statement.   Most of the work is done by calling
 * measure_parse_stdParams.
 * ----------------------------------------------------------------- */
static int measure_parse_when (
    MEASUREPTR meas,        /* in : measurement structure */
    wordlist *wl,       /* in : word list to parse */
    char *errBuf        /* in/out: buffer where we write error messages */
) {
   int pCnt, err = 0;
   char *p, *pVar1, *pVar2;
   meas->m_vec = NULL;
   meas->m_vec2 = NULL;
   meas->m_val = -1;
   meas->m_cross = -1;
   meas->m_fall = -1;
   meas->m_rise = -1;
   meas->m_td = 0;
   meas->m_from = 0.0e0;
   meas->m_to = 0.0e0;
   meas->m_at = -1;

   pCnt =0;
   while (wl) {
      p= wl->wl_word;

      if (pCnt == 0) {
         pVar1 = strtok(p, "=");
         pVar2 = strtok(NULL, "=");

         if (pVar2 == NULL) {
            sprintf(errBuf,"bad syntax\n");
            return 0;
         }

         meas->m_vec = copy(pVar1);
         /* correct for vectors like vm, vp etc. */
         if (cieq("ac", meas->m_analysis) || cieq("sp", meas->m_analysis)) 
            correct_vec(meas);
         if (measure_valid_vector(pVar2)==1) {
            meas->m_vec2 = copy(pVar2);
            /* correct for vectors like vm, vp etc. */
            if (cieq("ac", meas->m_analysis) || cieq("sp", meas->m_analysis)) 
               correct_vec(meas);
         }
         else
            meas->m_val = INPevaluate( &pVar2, &err, 1 );
      } else {
         if (measure_parse_stdParams(meas, wl, NULL, errBuf) == 0)
            return 0;
         break;
      }

      wl = wl->wl_next;
      pCnt ++;
   }
   return 1;
}


/* -----------------------------------------------------------------
 * Function: Given a wordlist and measurement structure, parse a
 * TRIGGER or TARGET clause of a measurement statement.   Most of the
 * work is done by calling measure_parse_stdParams.
 * ----------------------------------------------------------------- */
static int measure_parse_trigtarg (
    MEASUREPTR meas,        /* in : measurement structure */
    wordlist *words,        /* in : word list to parse */
    wordlist *wlTarg,       /* out : where we stopped parsing target clause */
    char *trigTarg,         /* in : type of clause */
    char *errbuf        /* in/out: buffer where we write error messages */
) {

   int pcnt;
   char *p;

   meas->m_vec = NULL;
   meas->m_vec2 = NULL;
   meas->m_cross = -1;
   meas->m_fall = -1;
   meas->m_rise = -1;
   meas->m_td = 0;
   meas->m_from = 0.0e0;
   meas->m_to = 0.0e0;
   meas->m_at = -1;

   pcnt =0;
      while (words != wlTarg) {
         p = words->wl_word;

         if ((pcnt == 0) && !ciprefix("at", p)) {
            meas->m_vec= cp_unquote(words->wl_word);
            /* correct for vectors like vm, vp etc. */
            if (cieq("ac", meas->m_analysis) || cieq("sp", meas->m_analysis))
               correct_vec(meas);
         } else if (ciprefix("at", p)) {
            if (measure_parse_stdParams(meas, words, wlTarg, errbuf) == 0)
               return 0;
         } else {
            if (measure_parse_stdParams(meas, words, wlTarg, errbuf) == 0)
               return 0;
         break;

      }

      words = words->wl_next;
      pcnt ++;
   }

   if (pcnt == 0) {
      sprintf(errbuf,"bad syntax of '%s'\n", trigTarg);
      return 0;
   }

   // valid vector
   if (measure_valid_vector(meas->m_vec)==0) {
      sprintf(errbuf,"no such vector as '%s'\n", meas->m_vec);
      return 0;
   }

   return 1;
}

/* -----------------------------------------------------------------
 * Function: Given a wordlist, extract the measurement statement,
 * process it, and return a result.  If out_line is furnished, we
 * format and copy the result it this string buffer.  The autocheck
 * variable allows us to check for "autostop".  This function is
 * called from measure.c.    We use the functions in this file because
 * the parsing is much more complete and thorough.
* ----------------------------------------------------------------- */
int
get_measure2(
   wordlist *wl,       /* in: a word list for us to process */
   double *result,     /* out : the result of the measurement */
   char *out_line,         /* out: formatted result - may be NULL */
   bool autocheck      /* in: TRUE if checking for "autostop"; FALSE otherwise */
) {
   wordlist *words, *wlTarg, *wlWhen;
   char errbuf[100];
   char *mAnalysis = NULL;             // analysis type
   char *mName = NULL;             // name given to the measured output
   char *mFunction = NULL;
   int precision;          // measurement precision
   int mFunctionType, wl_cnt;
   char *p;

   mFunctionType = -1;
   *result = 0.0e0;        /* default result */

   if (!wl) {
      printf("usage: measure .....\n");
      return MEASUREMENT_FAILURE;
   }

   if (!plot_cur || !plot_cur->pl_dvecs || !plot_cur->pl_scale) {
      fprintf(cp_err, "Error: no vectors available\n");
      return MEASUREMENT_FAILURE;
   }

   if (!ciprefix("tran", plot_cur->pl_typename) && !ciprefix("ac", plot_cur->pl_typename) &&
         !ciprefix("dc", plot_cur->pl_typename) && !ciprefix("sp", plot_cur->pl_typename)) {
      fprintf(cp_err, "Error: measure limited to tran, dc or ac analysis\n");
      return MEASUREMENT_FAILURE;
   }

   words =wl;
   wlTarg = NULL;
   wlWhen = NULL;

   if (!words) {
      fprintf(cp_err, "Error: no assignment found.\n");
      return MEASUREMENT_FAILURE;
   }

   precision = measure_get_precision() ;
   wl_cnt = 0;
   while (words) {

      switch(wl_cnt)
      {
         case 0:
            mAnalysis = cp_unquote(words->wl_word);
            break;
         case 1:
            mName = cp_unquote(words->wl_word);
            break;
         case 2:
         {
            mFunctionType = measure_function_type(words->wl_word);
            if ( mFunctionType == AT_UNKNOWN ){
               if(!(autocheck)){
                  printf("\tmeasure '%s'  failed\n", mName);
                  printf("Error: measure  %s  :\n", mName);
                  printf("\tno such function as '%s'\n", words->wl_word);
               }
               return MEASUREMENT_FAILURE;
            }
            break;
         }
         default:
         {
            p = words->wl_word;

            if (strcasecmp(p,"targ")==0)
               wlTarg = words;

               if (strcasecmp(p,"when")==0)
                  wlWhen = words;

               break;
         }
      }
      wl_cnt ++;
      words = words->wl_next;
   }

   if (wl_cnt < 3) {
      printf("\tmeasure '%s'  failed\n", mName);
      printf("Error: measure  %s  :\n", mName);
      printf("\tinvalid num params\n");
      return MEASUREMENT_FAILURE;
   }

    //------------------------


   words =wl;

   if (words)
      words = words->wl_next; // skip
   if (words)
      words = words->wl_next; // skip results name
   if (words)
      words = words->wl_next; // Function

    // switch here
   switch(mFunctionType)
   {
      case AT_DELAY:
      case AT_TRIG:
      {
         // trig parameters
         MEASUREPTR measTrig, measTarg;
         measTrig = (struct measure*)tmalloc(sizeof(struct measure));
         measTarg = (struct measure*)tmalloc(sizeof(struct measure));

         measTrig->m_analysis = measTarg->m_analysis = mAnalysis;

         if (measure_parse_trigtarg(measTrig, words , wlTarg, "trig", errbuf)==0) {
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         if ((measTrig->m_rise == -1) && (measTrig->m_fall == -1) && 
               (measTrig->m_cross == -1) && (measTrig->m_at == -1)) {
            sprintf(errbuf,"at, rise, fall or cross must be given\n");
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         while (words != wlTarg)
            words = words->wl_next; // hack

         if (words)
            words = words->wl_next; // skip targ

         if (measure_parse_trigtarg(measTarg, words , NULL, "targ", errbuf)==0) {
            measure_errMessage(mName, mFunction, "TARG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         if ((measTarg->m_rise == -1) && (measTarg->m_fall == -1) && 
               (measTarg->m_cross == -1)&& (measTarg->m_at == -1)) {
            sprintf(errbuf,"at, rise, fall or cross must be given\n");
            measure_errMessage(mName, mFunction, "TARG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // measure trig
         if (measTrig->m_at == -1)
            com_measure_when(measTrig);
         else
            measTrig->m_measured = measTrig->m_at;


         if (measTrig->m_measured == 0.0e0) {
             sprintf(errbuf,"out of interval\n");
             measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
             return MEASUREMENT_FAILURE;
         }
         // measure targ
         com_measure_when(measTarg);

         if (measTarg->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "TARG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // print results
         if( out_line ){
            sprintf(out_line,"%-20s=  %e targ=  %e trig=  %e\n", mName, (measTarg->m_measured - measTrig->m_measured), measTarg->m_measured, measTrig->m_measured);
         } else {
            printf("%-20s=  %e targ=  %e trig=  %e\n", mName, (measTarg->m_measured - measTrig->m_measured), measTarg->m_measured, measTrig->m_measured);
         }

         *result = (measTarg->m_measured - measTrig->m_measured);
         return MEASUREMENT_OK;
      }
      case AT_FIND:
      {
         MEASUREPTR meas, measFind;
         meas = (struct measure*)tmalloc(sizeof(struct measure));
         measFind = (struct measure*)tmalloc(sizeof(struct measure));

         meas->m_analysis = measFind->m_analysis = mAnalysis;

         if (measure_parse_find(meas, words, wlWhen, errbuf) == 0) {
            measure_errMessage(mName, mFunction, "FIND", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         if (meas->m_at == -1 ) {
            // find .. when statment
            while (words != wlWhen)
               words = words->wl_next; // hack

            if (words)
               words = words->wl_next; // skip targ

            if (measure_parse_when(measFind, words, errbuf) ==0) {
               measure_errMessage(mName, mFunction, "WHEN", errbuf, autocheck);
               return MEASUREMENT_FAILURE;
            }

            com_measure_when(measFind);

            if (measFind->m_measured == 0.0e0) {
               sprintf(errbuf,"out of interval\n");
               measure_errMessage(mName, mFunction, "WHEN", errbuf, autocheck);
               return MEASUREMENT_FAILURE;
            }

            measure_at(measFind, measFind->m_measured);
            meas->m_measured = measFind->m_measured;

         } else {
            measure_at(meas, meas->m_at);
         }

         if (meas->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "WHEN", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // print results
         if( out_line ){
            sprintf(out_line,"%-20s=  %e\n", mName, meas->m_measured);
         } else {
            printf("%-20s=  %e\n", mName, meas->m_measured);
         }
         *result = meas->m_measured;
         return MEASUREMENT_OK;
      }
      case AT_WHEN:
      {
         MEASUREPTR meas;
         meas = (struct measure*)tmalloc(sizeof(struct measure));
         meas->m_analysis = mAnalysis;
         if (measure_parse_when(meas, words, errbuf) ==0) {
            measure_errMessage(mName, mFunction, "WHEN", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         com_measure_when(meas);

         if (meas->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "WHEN", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // print results
         if( out_line ){
            sprintf(out_line,"%-20s=   %.*e\n", mName, precision, meas->m_measured);
         } else {
            printf("%-20s=  %e\n", mName, meas->m_measured);
         }

         *result = meas->m_measured;
         return MEASUREMENT_OK;
      }
      case AT_RMS:
      case AT_INTEG:
      {
         // trig parameters
         MEASUREPTR meas;
         meas = (struct measure*)tmalloc(sizeof(struct measure));
         meas->m_analysis = mAnalysis;
         if (measure_parse_trigtarg(meas, words , NULL, "trig", errbuf)==0) {
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // measure
         measure_rms_integral(meas,mFunctionType);

         if (meas->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck); // ??
            return MEASUREMENT_FAILURE;
         }

         if (meas->m_at == -1)
             meas->m_at = 0.0e0;

         // print results
         if( out_line ){
            sprintf(out_line,"%-20s=   %.*e from=  %.*e to=  %.*e\n", mName, precision, meas->m_measured, precision, meas->m_from, precision, meas->m_to);
         } else {
            printf("%-20s=  %.*e from=  %.*e to=  %.*e\n", mName, precision, meas->m_measured, precision, meas->m_from, precision, meas->m_to);
         }
         *result=meas->m_measured;
         return MEASUREMENT_OK;

      }
      case AT_AVG:
      {
        // trig parameters
         MEASUREPTR meas;
         meas = (struct measure*)tmalloc(sizeof(struct measure));

         meas->m_analysis = mAnalysis;

         if (measure_parse_trigtarg(meas, words , NULL, "trig", errbuf)==0) {
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // measure
         measure_minMaxAvg(meas, mFunctionType);
         if (meas->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck); // ??
            return MEASUREMENT_FAILURE;
         }

         if (meas->m_at == -1)
            meas->m_at = meas->m_from;

         // print results
         if( out_line ){
             sprintf(out_line,"%-20s=  %e from=  %e to=  %e\n", mName, meas->m_measured, meas->m_at, meas->m_measured_at);
         } else {
             printf("%-20s=  %e from=  %e to=  %e\n", mName, meas->m_measured, meas->m_at, meas->m_measured_at);
         }
         *result=meas->m_measured;
         return MEASUREMENT_OK;

      }
      case AT_MIN:
      case AT_MAX:
      case AT_MIN_AT:
      case AT_MAX_AT:
      {
         // trig parameters
         MEASUREPTR measTrig;
         measTrig = (struct measure*)tmalloc(sizeof(struct measure));
         measTrig->m_analysis = mAnalysis;
         if (measure_parse_trigtarg(measTrig, words , NULL, "trig", errbuf)==0) {
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // measure
         if ((mFunctionType == AT_MIN) || (mFunctionType == AT_MIN_AT))
            measure_minMaxAvg(measTrig, AT_MIN);
         else
            measure_minMaxAvg(measTrig, AT_MAX);

         if (measTrig->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck); // ??
            return MEASUREMENT_FAILURE;
         }

         if ((mFunctionType == AT_MIN) || (mFunctionType == AT_MAX)) {
            // print results
            if( out_line ){
               sprintf(out_line,"%-20s=  %e at=  %e\n", mName, measTrig->m_measured, measTrig->m_measured_at);
            } else {
               printf("%-20s=  %e at=  %e\n", mName, measTrig->m_measured, measTrig->m_measured_at);
            }
            *result=measTrig->m_measured;
         } else {
            // print results
            if( out_line ){
               sprintf(out_line,"%-20s=  %e with=  %e\n", mName, measTrig->m_measured_at, measTrig->m_measured);
            } else {
               printf("%-20s=  %e with=  %e\n", mName, measTrig->m_measured_at, measTrig->m_measured);
            }
            *result=measTrig->m_measured_at;
         }
         return MEASUREMENT_OK;
      }
      case AT_PP:
      {
         double minValue, maxValue;
         MEASUREPTR measTrig;
         measTrig = (struct measure*)tmalloc(sizeof(struct measure));
         measTrig->m_analysis = mAnalysis;
         if (measure_parse_trigtarg(measTrig, words , NULL, "trig", errbuf)==0) {
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck);
            return MEASUREMENT_FAILURE;
         }

         // measure min
         measure_minMaxAvg(measTrig, AT_MIN);
         if (measTrig->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck); // ??
            return MEASUREMENT_FAILURE;
         }
         minValue = measTrig->m_measured;

         // measure max
         measure_minMaxAvg(measTrig, AT_MAX);
         if (measTrig->m_measured == 0.0e0) {
            sprintf(errbuf,"out of interval\n");
            measure_errMessage(mName, mFunction, "TRIG", errbuf, autocheck); // ??
            return MEASUREMENT_FAILURE;
         }
         maxValue = measTrig->m_measured;

         // print results
         if( out_line ){
            sprintf(out_line,"%-20s=  %e from=  %e to=  %e\n", mName, (maxValue - minValue), measTrig->m_from, measTrig->m_to);
         } else {
            printf("%-20s=  %e from=  %e to=  %e\n", mName, (maxValue - minValue), measTrig->m_from, measTrig->m_to);
         }
         *result = (maxValue - minValue);
         return MEASUREMENT_OK;
      }
      
      case AT_DERIV:
      case AT_ERR:
      case AT_ERR1:
      case AT_ERR2:
      case AT_ERR3:
      {
         printf("\tmeasure '%s'  failed\n", mName);
         printf("Error: measure  %s  :\n", mName);
         printf("\tfunction '%s' currently not supported\n", mFunction);
         break;
      }
   }
   return MEASUREMENT_FAILURE;
}
