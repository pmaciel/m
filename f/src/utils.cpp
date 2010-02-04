
#include <cstdlib>
#include "utils.h"


// -- miscelaneous ----

// copy one array into another
void copy(double *U1, double *U2, int N)
{
  for (int i=0; i<N; ++i)
    U2[i] = U1[i];
}


// error handler
void nrerror(std::string msg)
{
  std::cerr << "ERROR: " << msg << std::endl;
  throw 42;
}


// 3d vector product
void vecprd3(double *U1, double *U2, double *VP)
{
  VP[0] = U1[1]*U2[2] - U1[2]*U2[1];
  VP[1] = U1[2]*U2[0] - U1[0]*U2[2];
  VP[2] = U1[0]*U2[1] - U1[1]*U2[0];
}


// -- memory allocation ----

/* allocate a double 3tensor with range t[nrow][ncol][ndep] */
double ***d3tensor(long nrow, long ncol, long ndep)
{
  double ***t = new double**[nrow+1];
  if (!t)
    nrerror("allocation failure 1 in d3tensor()");
  t += 1;

  t[0] = new double*[nrow*ncol+1];
  if (!t[0])
    nrerror("allocation failure 2 in d3tensor()");
  t[0] += 1;

  t[0][0] = new double[nrow*ncol*ndep+1];
  if (!t[0][0])
    nrerror("allocation failure 3 in d3tensor()");
  t[0][0] += 1;

  for (long j=1; j<ncol; j++)
    t[0][j] = t[0][j-1] + ndep;
  for (long i=1; i<nrow; i++) {
    t[i]    = t[i-1]    + ncol;
    t[i][0] = t[i-1][0] + ncol*ndep;
    for (long j=1; j<ncol; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  for (long i=0; i<nrow; ++i)
    for (long j=0; j<ncol; ++j)
      for (long k=0; k<ndep; ++k)
        t[i][j][k] = 0.;
  return t;
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m = (double **) malloc((size_t)((nrow+1)*sizeof(double*)));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;

  m[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
  if (!m[nrl])
    nrerror("allocation failure 2 in matrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)
    m[i]=m[i-1]+ncol;

  return m;
}


/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh)
{
  double *v = (double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v-nl+1;
}


// -- memory deallocation ----

/* free a double d3tensor allocated by d3tensor() */
void free_d3tensor(double ***t)
{
  delete[] (t[0][0]-1);
  delete[] (t[0]-1);
  delete[] (t-1);
}


/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((char*) (m[nrl]+ncl-1));
  free((char*) (m+nrl-1));
}


/* free a double vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh)
{
  free((char*) (v+nl-1));
}

