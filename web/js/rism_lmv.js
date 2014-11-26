/* Labik-Malijevsky-Vonka solver */



// compute Cjk
function getCjk(Cjk, npt, M, ns, der, costab, dp)
{
  var i, j, ij, ji, m, k, l;

  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;

      for ( m = 1; m < 3*M - 1; m++ ) {
        for ( dp[m] = 0, l = 0; l < npt; l++ )
          dp[m] += der[ij][l] * costab[m][l];
        dp[m] /= npt;
      }

      for ( m = 0; m < M; m++ )
        for ( k = 0; k < M; k++ )
          Cjk[ij][m*M+k] = dp[k-m+M] - dp[k+m+M];

      if ( j == i ) continue;

      ji = j*ns + i;
      for ( m = 0; m < M*M; m++ )
        Cjk[ji][m] = Cjk[ij][m];
    }
  }
}



// compute the Jacobian matrix for the Newton-Raphson method
function getjacob(mat, b, M, npr, ns,
    prmask, ntk, tk, Cjk, invwc1w)
{
  var i1, j1, ipr1, m1, id1, i2, j2, ipr2, m2, id2, Mp;
  var y;

  Mp = M * npr;
  for ( m1 = 0; m1 < M; m1++ )
    for ( ipr1 = 0, i1 = 0; i1 < ns; i1++ )
      for ( j1 = i1; j1 < ns; j1++, ipr1++ ) {
        id1 = m1*npr + ipr1;
        if ( prmask[i1*ns + j1] )
          b[id1] = fft_ki[m1] * (ntk[i1*ns+j1][m1] - tk[i1*ns+j1][m1]);
        else
          b[id1] = 0;

        for ( m2 = 0; m2 < M; m2++ )
          for ( ipr2 = 0, i2 = 0; i2 < ns; i2++ )
            for ( j2 = i2; j2 < ns; j2++, ipr2++ ) {
              id2 = m2*npr + ipr2;
              y = (ipr1 == ipr2 ? (m1 == m2) + Cjk[i1*ns+j1][m1*M+m2] : 0)
                - invwc1w[i2*ns+i1][m1] * Cjk[i2*ns+j2][m1*M+m2]
                * invwc1w[j2*ns+j1][m2];
              mat[id1*Mp + id2] = y;
            }
      }
}



// Reference:
// Stanislav Labik, Anatol Malijevsky, Petr Vonka
// A rapidly convergent method of solving the OZ equation
// Molecular Physics, 1985, Vol. 56, No. 3, 709-715
function iter_lmv(vrsr, wk, cr, der, ck, vklr,
    tr, tk, ntk, invwc1w, uv)
{
  var i, j, l, ij, it, M, npr, ipr, Mp;
  var Cjk, mat, b, costab, dp;
  var y, err1 = 0, err2 = 0, err = 0, errp = errinf, dmp;

  // initialize t(k) and t(r)
  sphr_r2k(cr, ck, ns, null);
  oz(ck, vklr, tk, wk, null); // c(k) --> t(k)
  sphr_k2r(tk, tr, ns, null);

  // set the optimal M
  M = get_int("lmv_M", 0);
  M = (M > 0) ? M : (2 * rmax/sigma[ns-1]);
  if ( M >= npt ) M = npt;
  if ( verbose ) console.log("select M", M);

  // set the damping factor
  dmp = get_float("lmv_damp", 1);

  npr = ns * (ns + 1) / 2;
  Mp = M * npr;

  /* initialize the cosine table */
  if ( M > 0 ) {
    Cjk = newarr2d(ns*ns, M*M);
    mat = newarr(Mp*Mp);
    b = newarr(Mp);
    dp = newarr(3*M);
    costab = newarr2d(3*M, npt);
    for ( j = 0; j < 3*M; j++ )
      for ( i = 0; i < npt; i++ )
        costab[j][i] = Math.cos(PI*(i+.5)*(j-M)/npt);
  }

  for ( it = 0; it < itmax; it++ ) {
    // compute c(r) and c(k) from the closure
    err = closure(null, der, vrsr, cr, tr, uv.prmask, true, 1.0);
    sphr_r2k(cr, ck, ns, null);

    // compute Cjk
    getCjk(Cjk, npt, M, ns, der, costab, dp);

    // compute the new t(k)
    oz(ck, vklr, ntk, wk, invwc1w);

    // compute the Jacobian matrix for the Newton-Raphson method
    getjacob(mat, b, M, npr, ns, uv.prmask, ntk, tk, Cjk, invwc1w);

    if ( lusolve(mat, b, Mp, 1e-14) != 0 )
      throw new Error("LU solve failed: stage " + uv.stage + ", it " + it);

    // compute the new t(k)
    err1 = err2 = 0;
    for ( ipr = 0, i = 0; i < ns; i++ )
      for ( j = i; j < ns; j++, ipr++ ) {
        ij = i*ns + j;
        if ( !uv.prmask[ij] ) continue;

        for ( l = 0; l < npt; l++ ) {
          if ( l < M ) {
            // use the Newton-Raphson method to solve for t(k) of small k
            y = b[l*npr+ipr] / fft_ki[l];
            if ( Math.abs(y) > err1 ) err1 = Math.abs(y);
          } else {
            // use the OZ relation to solve for t(k) of large k
            y = ntk[ij][l] - tk[ij][l];
            if ( Math.abs(y) > err2 ) err2 = Math.abs(y);
          }

          tk[ij][l] += dmp * y;
          if ( j > i ) tk[j*ns+i][l] = tk[ij][l];
        }
      }

    sphr_k2r(tk, tr, ns, null);

    if ( verbose )
      console.log("it", it, "M", M, "err", errp, "->", err,
        "tk_err", err1, "/", err2, "damp", dmp);

    if ( err < tol ) {
      // switch between stages
      if ( uv_switch(uv) != 0 ) break;
      if ( uv.stage == SOLUTE_SOLUTE
        && uv.infdil && uv.atomicsolvent ) {
        step_picard(null, null, vrsr, wk,
            cr, ck, vklr, tr, tk, uv.prmask, 1, 1.);
        break; // no need to iterate further
      }
      it = -1;
      err = errinf;
    }
    errp = err;
  }
  return [err, it];
}



