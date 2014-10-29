#ifndef UTIL_H__
#define UTIL_H__



#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <fftw3.h>



#define PI 3.14159265358979323846

#define KB    1.3806488e-23       /* Boltzmann constant: m^2*kg/s^2/K = J/K */
#define NA    6.02214129e23       /* Avogadro constant: mol^{-1} */
#define EC    1.602176565e-19     /* elementary charge: A*s */
#define EPS0  8.854187817620e-12  /* vacuum permittivity A^2*s^4/kg/m^3 */
#define CALPJ 4.184               /* calories per joule */

#define KC    (1/(4*PI*EPS0))
#define KE2   (NA*KC*EC*EC*1e7)   /* NA e^2 / (4 PI EPS0) in angstrom*kJ/mol, 1389.354578 */
#define KE2C  (KE2/CALPJ)         /* NA e^2 / (4 PI EPS0) in angstrom*kcal/mol, 322.0637137 */

#define KBNA  (KB*NA*0.001) /* Boltzmann constant in kJ/mol/K,   0.00831446214547 */
#define KBNAC (KBNA/CALPJ)  /* Boltzmann constant in kcal/mol/K, 0.00198720414667 */

#define KE2PK (KE2/KBNA)    /* e^2/(4 Pi EPS0)/kB in angstrom*K, 167100.9579 */



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %d\n", #x, (int) (n)); \
    exit(1); } }
#endif

#define newarr(x, n) x = (double *) fftw_malloc(sizeof(x[0]) * n)
#define delarr(x) fftw_free(x)

/* allocate a two-dimensional array */
#define newarr2d(x, m, n) { int i_; \
  x = malloc(sizeof(x[0]) * m); \
  for ( i_ = 0; i_ < m; i_++ ) newarr( x[i_], n ); }

/* free a two-dimensional array */
#define delarr2d(x, m) { int i_; \
  for ( i_ = 0; i_ < m; i_++ ) delarr( x[i_] ); \
  free(x); x = NULL; }

/* copy array */
#define cparr2d(x, y, m, n) { int i_, k_; \
  for ( i_ = 0; i_ < m; i_++ ) \
    for ( k_ = 0; k_ < n; k_++ ) \
      x[i_][k_] = y[i_][k_]; }



fftw_plan fftwplan;
double *fft_arr;
double *fft_ri, *fft_ki;
double fft_dr, fft_dk;
int fft_npt;



/* initialize fftw */
static void initfftw(double rmax, int npt)
{
  int i;

  fft_npt = npt;
  fft_dr = rmax / npt;
  fft_dk = (2*PI) / (2*npt*fft_dr);
  newarr(fft_arr, npt);
  fftwplan = fftw_plan_r2r_1d(npt, fft_arr, fft_arr,
      FFTW_RODFT11, FFTW_ESTIMATE);
  if ( fftwplan == NULL ) exit(1);
  newarr(fft_ri, npt);
  newarr(fft_ki, npt);
  for ( i = 0; i < npt; i++ ) {
    fft_ri[i] = (i + .5) * fft_dr;
    fft_ki[i] = (i + .5) * fft_dk;
  }
}



/* clean up fftw */
static void donefftw(void)
{
  delarr(fft_arr);
  fftw_destroy_plan(fftwplan);
  fftw_cleanup();
  free(fft_ri);
  free(fft_ki);
}



/* c(r) --> c(k)
 * c(k) = 2 Pi/k Int 2 r c(r) sin(kr) dr */
#define sphr_r2k(cr, ck, ns, prmask) \
  sphrt(cr, ck, fft_ri, fft_ki, fft_dr*(2*PI), ns, prmask)

/* t(k) --> t(r)
 * t(r) = 2 Pi/r/(2 Pi)^3 Int 2 k t(k) sin(kr) dk */
#define sphr_k2r(tk, tr, ns, prmask) \
  sphrt(tk, tr, fft_ki, fft_ri, fft_dk/(4*PI*PI), ns, prmask)

/* f(k) = fac/k Int 2 r f(r) sin(kr) dr */
static void sphrt(double **fr, double **fk,
    const double *ri, const double *ki, double fac,
    int ns, const int *prmask)
{
  int i, j, ij, l;

  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i * ns + j;
      if ( prmask && !prmask[ij] ) continue;
      for ( l = 0; l < fft_npt; l++ )
        fft_arr[l] = fr[ij][l] * ri[l];
      fftw_execute(fftwplan);
      for ( l = 0; l < fft_npt; l++ ) {
        fk[ij][l] = fft_arr[l] * fac / ki[l];
        if ( j > i ) fk[j*ns + i][l] = fk[ij][l];
      }
    }
  }
}



/* compute the inverse matrix b = a^(-1), by Gaussian elimination */
static int invmat(int n, double *a, double *b)
{
  int i, j, k, ip;
  double x;

  /* initialize b as the identity matrix */
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ )
      b[i*n + j] = (i == j);

  /* Gaussian elimination */
  for ( i = 0; i < n; i++ ) {
    /* choose the pivot as the largest element of column i */
    x = fabs( a[(ip = i)*n + i] );
    for ( k = ip + 1; k < n; k++ )
      if ( fabs( a[k*n + i] ) > x )
        x = fabs( a[(ip = k)*n + i] );

    /* swap the pivot (ip'th) row with the present row i */
    for ( k = i; k < n; k++ )
      x = a[i*n + k], a[i*n + k] = a[ip*n + k], a[ip*n + k] = x;
    for ( k = 0; k < n; k++ )
      x = b[i*n + k], b[i*n + k] = b[ip*n + k], b[ip*n + k] = x;

    /* normalize this row */
    x = a[i*n + i];
    if ( fabs(x) < DBL_MIN ) {
      fprintf(stderr, "Error: singular matrix of %dx%d\n", n, n);
      return -1;
    }
    for ( k = i; k < n; k++ ) a[i*n + k] /= x;
    for ( k = 0; k < n; k++ ) b[i*n + k] /= x;

    /* use the pivot row to zero the rest rows */
    for ( j = i + 1; j < n; j++ ) {
      x = a[j*n + i];
      for ( k = i; k < n; k++ )
        a[j*n + k] -= x * a[i*n + k];
      for ( k = 0; k < n; k++ )
        b[j*n + k] -= x * b[i*n + k];
    }
  }

  /* now that the matrix is upper triangular
   * make it diagonal */
  for ( i = n - 1; i >= 0; i-- ) {
    /* note a[i*n + i] should be 1 now */
    for ( j = 0; j < i; j++ ) {
      x = a[j*n + i];
      for ( k = 0; k < n; k++ )
        b[j*n + k] -= b[i*n + k] * x;
    }
  }
  return 0;
}



/* solve the linear equation: a x = b */
static int linsolve(int n, double *a, double *x, double *b)
{
  int i, j, k, ip;
  double y;

  for ( i = 0; i < n; i++ ) x[i] = 0;

  for ( i = 0; i < n; i++ ) {
    /* 1. select the pivot of the ith column
     * pivot: the maximal element a(j = i..n-1, i)  */
    y = fabs( a[(ip = i)*n + i] );
    for ( j = i + 1; j < n; j++ )
      if ( fabs( a[j*n + i] ) > y )
        y = fabs( a[(ip = j)*n + i] );

    /* 2. swap the pivot (ip'th) row with the ith row */
    if ( ip != i ) {
      y = b[ip], b[ip] = b[i], b[i] = y;
      for ( j = i; j < n; j++ )
        y = a[i*n + j], a[i*n + j] = a[ip*n + j], a[ip*n + j] = y;
    }
    y = a[i*n + i];
    if ( fabs(y) < DBL_MIN ) {
      fprintf(stderr, "singular matrix on %dth row\n", i);
      return -1;
    }

    /* 3. normalize the ith row */
    b[i] /= y;
    for ( k = i; k < n; k++ )
      a[i*n + k] /= y;

    /* 4. use the pivot row to eliminate the following rows */
    for ( j = i + 1; j < n; j++ ) { /* for rows */
      y = a[j*n + i];
      b[j] -= y * b[i];
      for ( k = i; k < n; k++ )
        a[j*n + k] -= y * a[i*n + k];
    }
  }

  /* 5. now that the matrix is upper-triangular
   *    solve for x */
  for ( i = n - 1; i >= 0; i-- ) {
    x[i] = b[i] / a[i*n + i];
    for ( j = 0; j < i; j++ ) {
      b[j] -= a[j*n + i] * x[i];
    }
  }
  return 0;
}



#endif /* UTIL_H__ */

