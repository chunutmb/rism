/* utility routines */



"use strict";



/* Numerical constants
 * default units:
 * length: angstrom
 * energy: kJ
 * */



var PI          = Math.PI;
var INFTY       = 1e30;
var errinf      = 1e20;

var CAL_TO_J    = 4.184;              // calories to joules
var J_TO_CAL    = 1/CAL_TO_J;         // joules to calories
var CAL_TO_KJ   = CAL_TO_J*1e-3;      // calories to kilo joules
var J_TO_KCAL   = J_TO_CAL*1e-3;      // joules to kilo calories
var KCAL_TO_J   = CAL_TO_J*1e3;       // kilo calories to joules
var KJ_TO_CAL   = J_TO_CAL*1e3;       // kilo joules to calories
var NA          = 6.02214129e23;      // Avogadro constant in mol^{-1}
var EC          = 1.602176565e-19;    // elementary charge in A*s
var EPS0_SI     = 8.854187817620e-12; // vacuum permittivity in A^2*s^4/kg/m^3

var KB_SI       = 1.3806488e-23;      // Boltzmann constant in m^2*kg/s^2/K = J/K
var KB_KJ       = KB_SI*0.001;        // Boltzmann constant in kJ/K, 1.3806488e-26
var KB_J        = KB_SI;              // Boltzmann constant in J/K, 1.3806488e-23
var KB_ERG      = KB_SI*1e7;          // Boltzmann constant in erg/K, 1.3806488e-16
var KB_KCAL     = KB_KJ/CAL_TO_J;     // Boltzmann constant in kcal/K, 3.29982982791587e-27
var KB_CAL      = KB_J/CAL_TO_J;      // Boltzmann constant in cal/K, 3.29982982791587e-24
var KB          = KB_KJ;              // Boltzmann constant in the default unit

var KBNA_SI     = KB_SI*NA;           // Boltzmann constant in J/mol/K, 8.31446214547
var KBNA_KJ     = KB_KJ*NA;           // Boltzmann constant in kJ/mol/K, 0.00831446214547
var KBNA_J      = KBNA_SI;            // Boltzmann constant in J/mol/K, 8.31446214547
var KBNA_ERG    = KBNA_SI*1e7;        // Boltzmann constant in erg/mol/K, 8.31446214547e7
var KBNA_KCAL   = KBNA_KJ/CAL_TO_J;   // Boltzmann constant in kcal/mol/K, 0.00198720414667
var KBNA_CAL    = KBNA_J/CAL_TO_J;    // Boltzmann constant in cal/mol/K, 1.98720414667
var KBNA        = KBNA_KJ;            // Boltzmann constant in the default unit
var KBNAC       = KBNA_KCAL;

var KC_SI       = 1/(4*PI*EPS0_SI);   // 1 / (4 PI EPS0) in kg*m^3/A^2/s^4 = J*m/A^2/s^2, 8.9875517873686e9

var KE2_SI      = KC_SI*EC*EC;       // e^2 / (4 PI EPS0) in m*J, 2.3070773523707e-28
var KE2_AKJ     = KE2_SI*1e7;        // e^2 / (4 PI EPS0) in angstrom*kJ, 2.3070773523707e-21
var KE2_AJ      = KE2_SI*1e10;       // e^2 / (4 PI EPS0) in angstrom*J, 2.3070773523707e-18
var KE2_AERG    = KE2_SI*1e17;       // e^2 / (4 PI EPS0) in angstrom*erg, 2.3070773523707e-11
var KE2_AKCAL   = KE2_AKJ/CAL_TO_J;  // e^2 / (4 PI EPS0) in angstrom*kcal, 5.514047209298996e-22
var KE2_ACAL    = KE2_AJ/CAL_TO_J;   // e^2 / (4 PI EPS0) in angstrom*cal, 5.514047209298996e-19
var KE2         = KE2_AKJ;

var KE2NA_SI    = KE2_SI*NA;         // e^2 / (4 PI EPS0) * NA in m*J/mol, 1.389354578e-7
var KE2NA_AKJ   = KE2_AKJ*NA;        // e^2 / (4 PI EPS0) * NA in angstrom*kJ/mol, 1389.354578
var KE2NA_AJ    = KE2_AJ*NA;         // e^2 / (4 PI EPS0) * NA in angstrom*J/mol, 1.389354578e6
var KE2NA_AERG  = KE2_AERG*NA;       // e^2 / (4 PI EPS0) * NA in angstrom*J/mol, 1.389354578e13
var KE2NA_AKCAL = KE2_AKCAL*NA;      // e^2 / (4 PI EPS0) * NA in angstrom*kcal/mol, 322.0637137
var KE2NA_ACAL  = KE2_ACAL*NA;       // e^2 / (4 PI EPS0) * NA in angstrom*cal/mol, 3.220637137e5
var KE2NA       = KE2NA_AKJ;         // e^2 / (4 PI EPS0) * NA in the default unit
var KE2NAC      = KE2NA_AKCAL;

var KE2PK_SI    = KE2_SI/KB_SI;      // e^2/(4 Pi EPS0)/kB in m*K, 1.671009566e-5
var KE2PK_A     = KE2PK_SI*1e10;     // e^2/(4 Pi EPS0)/kB in angstrom*K, 167100.9566
var KE2PK       = KE2PK_A;           // e^2/(4 Pi EPS0)/kB in the default unit

var ERG_TO_J           = (1e-7);                // 1e-7
var ERG_TO_CAL         = (1e-7*J_TO_CAL);       // 2.390057361376673e-8
var ERG_TO_KJ          = (1e-10);               // 1e-10
var ERG_TO_KCAL        = (1e-10*J_TO_CAL);      // 2.390057361376673e-11
var ERG_TO_JPMOL       = (NA*1e-7);             // 6.02214129e16
var ERG_TO_CALPMOL     = (NA*1e-7*J_TO_CAL);    // 1.43932631214914e16
var ERG_TO_KJPMOL      = (NA*1e-10);            // 6.02214129e13
var ERG_TO_KCALPMOL    = (NA*1e-10*J_TO_CAL);   // 1.43932631214914e13

var J_TO_ERG           = (1/ERG_TO_J);          // 1e7
var CAL_TO_ERG         = (1/ERG_TO_CAL);        // 4.184e7
var KJ_TO_ERG          = (1/ERG_TO_KJ);         // 1e10
var KCAL_TO_ERG        = (1/ERG_TO_KCAL);       // 4.184e10
var JPMOL_TO_ERG       = (1/ERG_TO_JPMOL);      // 1.66053892103219e-17
var CALPMOL_TO_ERG     = (1/ERG_TO_CALPMOL);    // 6.947694845598683e-17
var KJPMOL_TO_ERG      = (1/ERG_TO_KJPMOL);     // 1.66053892103219e-14
var KCALPMOL_TO_ERG    = (1/ERG_TO_KCALPMOL);   // 6.947694845598683e-14


/* Memory allocation routines */



var fft_arr;
var fft_ri, fft_ki;
var fft_dr, fft_dk;
var fft_npt;



/* initialize fftw */
function initfftw(rmax, npt)
{
  var i;

  fft_npt = npt;
  fft_dr = rmax / npt;
  fft_dk = (2*PI) / (2*npt*fft_dr);
  fft_arr = newarr(npt);
  fft_ri = newarr(npt);
  fft_ki = newarr(npt);
  for ( i = 0; i < npt; i++ ) {
    fft_ri[i] = (i + .5) * fft_dr;
    fft_ki[i] = (i + .5) * fft_dk;
  }
}



/* f(k) = fac/k Int 2 r f(r) sin(kr) dr */
function sphrt(fr, fk, ri, ki, fac, ns, prmask)
{
  var i, j, ij, l;

  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i * ns + j;
      if ( prmask !== null && !prmask[ij] ) continue;
      for ( l = 0; l < fft_npt; l++ )
        fft_arr[l] = fr[ij][l] * ri[l];
      sint11(fft_arr, fft_npt);
      for ( l = 0; l < fft_npt; l++ ) {
        fk[ij][l] = fft_arr[l] * fac / ki[l];
        if ( j > i ) fk[j*ns + i][l] = fk[ij][l];
      }
    }
  }
}



/* c(r) --> c(k)
 * c(k) = 2 Pi/k Int 2 r c(r) sin(kr) dr */
function sphr_r2k(cr, ck, ns, prmask)
{
  sphrt(cr, ck, fft_ri, fft_ki, fft_dr*(2*PI), ns, prmask);
}

/* t(k) --> t(r)
 * t(r) = 2 Pi/r/(2 Pi)^3 Int 2 k t(k) sin(kr) dk */
function sphr_k2r(tk, tr, ns, prmask)
{
  sphrt(tk, tr, fft_ki, fft_ri, fft_dk/(4*PI*PI), ns, prmask);
}





/* Linear algebra routines */



/* multiply two matrices c = a * b */
function matmul(c, a, b, n)
{
  var i, j, k, x;

  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ ) {
      x = 0;
      for ( k = 0; k < n; k++ )
        x += a[i*n + k] * b[k*n + j];
      c[i*n + j] = x;
    }
  return c;
}



/* compute the inverse matrix b = a^(-1), by Gaussian elimination */
function invmat(a, b, n)
{
  var i, j, k, ip, x;

  // initialize b as the identity matrix
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ )
      b[i*n + j] = (i == j);

  // Gaussian elimination
  for ( i = 0; i < n; i++ ) {
    // choose the pivot as the largest element of column i
    x = Math.abs( a[(ip = i)*n + i] );
    for ( k = ip + 1; k < n; k++ )
      if ( Math.abs( a[k*n + i] ) > x )
        x = Math.abs( a[(ip = k)*n + i] );

    // swap the pivot (ip'th) row with the present row i
    for ( k = i; k < n; k++ )
      x = a[i*n + k], a[i*n + k] = a[ip*n + k], a[ip*n + k] = x;
    for ( k = 0; k < n; k++ )
      x = b[i*n + k], b[i*n + k] = b[ip*n + k], b[ip*n + k] = x;

    // normalize this row
    x = a[i*n + i];
    if ( Math.abs(x) <= 0 ) {
      console.log("Error: singular matrix of " + n + "x" + n);
      return -1;
    }
    for ( k = i; k < n; k++ ) a[i*n + k] /= x;
    for ( k = 0; k < n; k++ ) b[i*n + k] /= x;

    // use the pivot row to zero the rest rows
    for ( j = i + 1; j < n; j++ ) {
      x = a[j*n + i];
      for ( k = i; k < n; k++ )
        a[j*n + k] -= x * a[i*n + k];
      for ( k = 0; k < n; k++ )
        b[j*n + k] -= x * b[i*n + k];
    }
  }

  // now that the matrix is upper triangular
  // make it diagonal
  for ( i = n - 1; i >= 0; i-- ) {
    // note a[i*n + i] should be 1 now
    for ( j = 0; j < i; j++ ) {
      x = a[j*n + i];
      for ( k = 0; k < n; k++ )
        b[j*n + k] -= b[i*n + k] * x;
    }
  }
  return 0;
}



// solve A x = b by L U decomposition
// the matrix `a' will be destroyed
// the vector `b' will be replaced by `x' on return
function lusolve(a, b, n, tiny)
{
  var i, j, k, ip = 0, x, max;

  for (i = 0; i < n; i++) {  // normalize each equation
    for (max = 0.0, j = 0; j < n; j++)
      if ((x = Math.abs(a[i*n + j])) > max)
        max = x;
    if (max <= 0) return -1;
    for (x = 1./max, j = 0; j < n; j++)
      a[i*n + j] *= x;
    b[i] *= x;
  }

  // step 1: A = L U, column by column
  for (j = 0; j < n; j++) {
    // matrix U
    for (i = 0; i < j; i++) {
      for (x = a[i*n + j], k = 0; k < i; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
    }

    // matrix L, diagonal of L are 1
    max = 0.0;
    ip = j;
    for (i = j; i < n; i++) {
      for (x = a[i*n + j], k = 0; k < j; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
      if (Math.abs(x) > max) {
        max = Math.abs(x);
        ip = i;
      }
    }

    if (j != ip) { // swap the pivot row with the jth row
      for (k = 0; k < n; k++)
        x = a[ip*n + k], a[ip*n + k] = a[j*n + k], a[j*n + k] = x;
      x = b[ip], b[ip] = b[j], b[j] = x;
    }
    if (Math.abs(a[j*n + j]) < tiny)
      a[j*n + j] = tiny;
    // divide by the pivot element, for the L matrix
    if (j != n - 1)
      for (x = 1./a[j*n + j], i = j + 1; i < n; i++)
        a[i*n + j] *= x;
  }

  // step 2: solve the equation L U x = b
  for (i = 0; i < n; i++) { // L y = b
    x = b[i];
    for (j = 0; j < i; j++) x -= a[i*n + j] * b[j];
    b[i] = x;
  }
  for (i = n - 1; i >= 0; i--) { // U x = y
    x = b[i];
    for (j = i + 1; j < n; j++) x -= a[i*n + j] * b[j];
    b[i] = x / a[i*n + i];
  }
  return 0;
}





/* String routines */



/* check if a string is a number */
function striscnum(s)
{
  return !isNaN(s);
}



/* check if string `s' starts with a substring `t' */
function strstartswith(s, t)
{
  return s.substr(0, t.length) == t;
}



/* compare two strings, insensitive to the cases */
function strieq(s, t)
{
  return s.toLowerCase() == t.toLowerCase();
}



/* convert to lower cases, ignore punctuations */
function striiconv(s)
{
  return s.toLowerCase().replace(/[^a-z0-9]+/ig, "");
}



/* compare two strings, insensitive to the cases, ignore punctuations */
function striieq(s, t)
{
  return striiconv(s) == striiconv(t);
}



/* Special functions */



function erf(x) {
    // constants
    var a1 =  0.254829592;
    var a2 = -0.284496736;
    var a3 =  1.421413741;
    var a4 = -1.453152027;
    var a5 =  1.061405429;
    var p  =  0.3275911;

    // Save the sign of x
    var sign = 1;
    if (x < 0) {
      sign = -1;
    }
    x = Math.abs(x);

    // A&S formula 7.1.26
    var t = 1.0/(1.0 + p*x);
    var y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);

    return sign*y;
}

