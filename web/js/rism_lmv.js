/* Labik-Malijevsky-Vonka (LMV) solver */



"use strict";



function LMV(ns, npt, M, ki)
{
  this.ns = ns;
  var ns2 = ns * ns;
  this.npr = ns * (ns + 1) / 2;
  this.npt = npt;
  if ( M >= npt ) M = npt;
  if ( verbose ) console.log("select M", M);
  this.M = M;
  this.ki = ki;
  this.tr1      = newarr2d(ns2, npt);
  this.tk1      = newarr2d(ns2, npt);
  this.invwc1w  = newarr2d(ns2, npt);
  this.der      = newarr2d(ns2, npt);
  this.crbest   = newarr2d(ns2, npt);
  this.errmin   = errinf;
  /* initialize the cosine table */
  if ( M > 0 ) {
    var Mp = M * this.npr;
    this.Cjk = newarr2d(ns2, M * M);
    this.mat = newarr(Mp * Mp);
    this.a = newarr(Mp);
    this.dp = newarr(3 * M);
    this.costab = newarr2d(3 * M, npt);
    for ( var j = 0; j < 3*M; j++ )
      for ( var i = 0; i < npt; i++ )
        this.costab[j][i] = Math.cos(PI*(i+.5)*(j-M)/npt);
  }
}



/* save a good cr */
LMV.prototype.savebest = function(cr, err)
{
  if ( err < this.errmin ) {
    cparr2d(this.crbest, cr, this.ns * this.ns, this.npt);
    this.errmin = err;
  }
}



/* compute Cjk = d ck / d tk */
LMV.prototype.getCjk = function()
{
  var ns = this.ns, npt = this.npt, M = this.M;
  var k, m, l;

  for ( var i = 0; i < ns; i++ ) {
    for ( var j = i; j < ns; j++ ) {
      var ij = i*ns + j;

      for ( m = 1; m < 3*M - 1; m++ ) {
        for ( this.dp[m] = 0, l = 0; l < npt; l++ )
          this.dp[m] += der[ij][l] * this.costab[m][l];
        this.dp[m] /= npt;
      }

      for ( m = 0; m < M; m++ )
        for ( k = 0; k < M; k++ )
          this.Cjk[ij][m*M+k] = this.dp[k-m+M] - this.dp[k+m+M];

      if ( j == i ) continue;

      var ji = j*ns + i;
      for ( m = 0; m < M*M; m++ )
        this.Cjk[ji][m] = this.Cjk[ij][m];
    }
  }
}



/* compute the Jacobian matrix for the Newton-Raphson method */
LMV.prototype.getjacob = function(tk, uv)
{
  var ns = this.ns, npr = this.npr, M = this.M;
  var Mp = M * npr;
  for ( var m1 = 0; m1 < M; m1++ )
    for ( var ipr1 = 0, i1 = 0; i1 < ns; i1++ )
      for ( var j1 = i1; j1 < ns; j1++, ipr1++ ) {
        var ij1 = i1*ns + j1;
        var id1 = m1*npr + ipr1;
        this.a[id1] = uv.prmask[ij1] ? this.ki[m1] * (this.tk1[ij1][m1] - tk[ij1][m1]) : 0;

        for ( var m2 = 0; m2 < M; m2++ )
          for ( var ipr2 = 0, i2 = 0; i2 < ns; i2++ )
            for ( var j2 = i2; j2 < ns; j2++, ipr2++ ) {
              var id2 = m2*npr + ipr2;
              this.mat[id1*Mp + id2] =
                (ipr1 == ipr2 ? (m1 == m2) + this.Cjk[ij1][m1*M+m2] : 0)
                - this.invwc1w[i2*ns+i1][m1] * this.Cjk[i2*ns+j2][m1*M+m2]
                * this.invwc1w[j2*ns+j1][m2];
            }
      }
}



/* update tk */
LMV.prototype.update = function(tk, dmp, uv)
{
  var ns = this.ns, npr = this.npr, npt = this.npt, M = this.M;
  var i, j, ij, l, ipr;
  var del;

  // compute d ck / d tk
  this.getCjk();

  // compute the Jacobian matrix for the Newton-Raphson method
  this.getjacob(tk, uv);

  if ( lusolve(this.mat, this.a, npr * M, 1e-10) !== 0 ) {
    console.log("LU solve failed: stage", uv.stage);
    return -1;
  }

  // compute the new t(k)
  this.err1 = this.err2 = 0;
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++, ipr++ ) {
      ij = i*ns + j;
      if ( !uv.prmask[ij] ) continue;

      for ( l = 0; l < npt; l++ ) {
        if ( l < M ) {
          // use the Newton-Raphson method to solve for t(k) of small k
          del = this.a[l*npr+ipr] / this.ki[l];
          if ( Math.abs(del) > this.err1 ) this.err1 = Math.abs(del);
        } else {
          // use the OZ relation to solve for t(k) of large k
          del = this.tk1[ij][l] - tk[ij][l];
          if ( Math.abs(del) > this.err2 ) this.err2 = Math.abs(del);
        }

        tk[ij][l] += dmp * del;
      }
      if ( j > i ) cparr(tk[j*ns + i], tk[ij], npt);
    }

  return 0;
}



/* LMV solver */
function iter_lmv(vrsr, wk, cr, der, ck, vklr,
    tr, tk, ntk, invwc1w, uv)
{
  var it;
  var err = 0, errp = errinf, dmp;

  // set the optimal M
  var M = get_int("lmv_M", 0);
  M = (M > 0) ? M : (2 * rmax/sigma[ns-1]);

  var lmv = new LMV(ns, npt, M, fft_ki);
  cparr2d(lmv.crbest, cr, ns2, npt);

  // initialize t(k) and t(r)
  step_picard(null, vrsr, wk, cr, ck, vklr, tr, tk, null, 0.);
  cparr2d(lmv.tk1, tk, ns2, npt);

  // set the damping factor
  dmp = get_float("lmv_damp", 1);

  for ( it = 0; it <= itmax; it++ ) {
    // compute the error of the current c(r) and c(k)
    sphr_k2r(lmv.tk1, lmv.tr1, ns, uv.prmask);
    err = closure(null, null, vrsr, cr, lmv.tr1, uv.prmask, 0.);
    lmv.savebest(cr, err);

    closure(null, lmv.der, vrsr, cr, tr, uv.prmask, 1.);
    sphr_r2k(cr, ck, ns, uv.prmask);
    oz(ck, vklr, lmv.tk1, wk, lmv.invwc1w);

    // compute the new tk
    if ( lmv.update(tk, dmp, uv) != 0 ) break;
    sphr_k2r(tk, tr, ns, uv.prmask);

    if ( verbose )
      console.log("it", it, "M", M, "err", errp, "->", err,
        "tk_err", lmv.err1, "/", lmv.err2, "damp", dmp);

    if ( err < tol || it === itmax ) {
      // use the best cr discovered so far
      cparr2d(cr, lmv.crbest, ns2, npt);
      // update the corresponding ck, tr, tk, and the error
      err = step_picard(null, vrsr, wk,
          cr, ck, vklr, tr, tk, uv.prmask, 0.);
      // switch between stages
      if ( uv.switchstage() !== 0 ) {
        if ( uv.uu1step )
          step_picard(null, vrsr, wk, cr, ck, vklr, tr, tk, uv.prmask, 1.);
        break; // no need to iterate further
      }
      err = step_picard(null, vrsr, wk,
          cr, ck, vklr, tr, tk, uv.prmask, 0.);
      cparr2d(lmv.tk1, tk, ns2, npt);
      lmv.errmin = err;
      it = -1;
    }
    errp = err;
  }
  return [err, it];
}



