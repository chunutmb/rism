/* solve A x = b by L U decomposition
 * on return, matrix `a' is destroyed, and vector `b' becomes `x' */
function lusolve(a, b, n, tiny)
{
  var i, j, k, ip = 0;
  var x, max;

  for (i = 0; i < n; i++) {  // normalize each equation
    for (max = 0.0, j = 0; j < n; j++)
      if ((x = Math.abs(a[i*n + j])) > max)
        max = x;
    if (max < 1e-14) return -1;
    for (x = 1.0/max, j = 0; j < n; j++)
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
    for (ip = i = j; i < n; i++) {
      for (x = a[i*n + j], k = 0; k < j; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
      if (Math.abs(x) >= max) {
        max = Math.abs(x);
        ip = i;
      }
    }

    if (j != ip) { // swap the pivot row with the jth row
      for (k = 0; k < n; k++)
        x = a[ip*n + k], a[ip*n + k] = a[j*n + k], a[j*n + k] = x;
      x = b[ip], b[ip] = b[j], b[j] = x;
    }
    if (Math.abs(a[j*n + j]) < tiny) a[j*n + j] = tiny;
    // divide by the pivot element, for the L matrix
    if (j != n - 1)
      for (x = 1.0/a[j*n + j], i = j + 1; i < n; i++)
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



/* invert matrix `a' as b = a^(-1)
 * on return, matrix `a' is destroyed */
function luinv(a, b, n, tiny)
{
  var i, j, k, ip = 0;
  var x, max;

  // initialize the matrix as the identity matrix
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      b[i*n + j] = (i == j);

  for (i = 0; i < n; i++) {  // normalize each equation
    for (max = 0.0, j = 0; j < n; j++)
      if ((x = Math.abs(a[i*n + j])) > max)
        max = x;
    if (max < tiny) return 1;
    for (x = 1.0/max, j = 0; j < n; j++)
      a[i*n + j] *= x;
    b[i*n + i] *= x;
  }

  // solve A = L U, column by column
  for (j = 0; j < n; j++) {
    // matrix U
    for (i = 0; i < j; i++) {
      for (x = a[i*n + j], k = 0; k < i; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
    }

    // matrix L, diagonal of L are 1
    ip = j;
    max = 0;
    for (i = j; i < n; i++) {
      for (x = a[i*n + j], k = 0; k < j; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
      if (Math.abs(x) >= max) {
        max = Math.abs(x);
        ip = i;
      }
    }

    if (j != ip) { // swap the pivot row with the jth row
      for (k = 0; k < n; k++) {
        x = a[ip*n + k], a[ip*n + k] = a[j*n + k], a[j*n + k] = x;
        x = b[ip*n + k], b[ip*n + k] = b[j*n + k], b[j*n + k] = x;
      }
    }
    if (Math.abs(a[j*n + j]) < tiny)
      a[j*n + j] = tiny;
    // divide by the pivot element, for the L matrix
    if (j != n - 1)
      for (x = 1.0/a[j*n + j], i = j + 1; i < n; i++)
        a[i*n + j] *= x;
  }

  for ( k = 0; k < n; k++ ) {
    // solve L U x = b_k
    for (i = 0; i < n; i++) { // L y = b
      x = b[i*n + k];
      for (j = 0; j < i; j++) x -= a[i*n + j] * b[j*n + k];
      b[i*n + k] = x;
    }
    for (i = n - 1; i >= 0; i--) { // U x = y
      x = b[i*n + k];
      for (j = i + 1; j < n; j++) x -= a[i*n + j] * b[j*n + k];
      b[i*n + k] = x / a[i*n + i];
    }
  }
  return 0;
}


