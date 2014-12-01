


/* modified direct inversion of the iterative subspace (MDIIS) method */



function MDIIS(ns, npt, mnb)
{
  var ns2, mnb1;
  ns2 = ns * (ns + 1) / 2;
  this.ns    = ns;
  this.npt   = npt;
  this.mnb   = mnb; // maximal number of bases
  mnb1 = mnb + 1;
  this.cr    = newarr2d(mnb1, ns2 * npt); // basis
  this.res   = newarr2d(mnb1, ns2 * npt); // residues
  this.mat   = newarr(mnb1 * mnb1); // correlations of residues
  this.mat2  = newarr(mnb1 * mnb1); // temporary matrix for LU decomposition
  this.coef  = newarr(mnb1); // coefficients
}



/* solve the coefficients of combination */
MDIIS.prototype.solve = function()
{
  var nb = this.nb, nb1 = this.nb + 1, mnb1 = this.mnb + 1, i, j;

  for ( i = 0; i < nb; i++ ) this.coef[i] = 0;
  this.coef[nb] = -1;
  /* copy the matrix, for the content is to be destroyed */
  for ( i = 0; i < nb1; i++ )
    for ( j = 0; j < nb1; j++ )
      this.mat2[i*nb1 + j] = this.mat[i*mnb1 + j];
  for ( i = 0; i < nb1; i++ )
    this.mat2[i*nb1 + nb] = this.mat2[nb*nb1 + i] = -1;
  this.mat2[nb*nb1 + nb] = 0;
  i = lusolve(this.mat2, this.coef, nb1, 1e-20);
  return 0;
}



/* construct the new c(r) */
MDIIS.prototype.gencr = function(cr, damp, uv)
{
  var ns = this.ns, npt = this.npt, nb = this.nb, npr = uv.npr;
  var i, j, ij, ipr, k, l, il, y, coef;

  for ( il = 0; il < npr * npt; il++ )
    this.cr[nb][il] = 0;
  for ( k = 0; k < nb; k++ ) {
    coef = this.coef[k];
    for ( il = 0; il < npr * npt; il++ )
      this.cr[nb][il] += coef * (this.cr[k][il] + damp * this.res[k][il]);
  }
  /* save this.cr[nb] to c(r) */
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv.prmask[ij = i*ns + j] ) continue;
      for ( l = 0; l < npt; l++ ) {
        cr[ij][l] = y = this.cr[nb][ipr*npt + l];
        if ( j > i) cr[j*ns + i][l] = y;
      }
      ipr++;
    }
}



/* compute the dot product */
function getdot(a, b, n)
{
  var i, x = 0;

  for ( i = 0; i < n; i++ ) x += a[i] * b[i];
  return x / n;
}



/* build the residue correlation matrix */
MDIIS.prototype.build = function(cr, res, uv)
{
  var i, j, ipr, l, ib, id, mnb, mnb1, ns = this.ns, npt = this.npt;

  this.nb = 1;
  mnb = this.mnb;
  mnb1 = this.mnb + 1;

  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv.prmask[i*ns + j] ) continue;
      for ( l = 0; l < npt; l++ ) {
        id = ipr*npt + l;
        this.cr[0][id] = cr[i*ns + j][l];
        this.res[0][id] = res[id];
      }
      ipr++;
    }

  this.mat[0] = getdot(this.res[0], this.res[0], uv.npr * npt);
  for ( ib = 0; ib < mnb; ib++ )
    this.mat[ib*mnb1 + mnb] = this.mat[mnb*mnb1 + ib] = -1;
  this.mat[mnb*mnb1 + mnb] = 0;
  return 0;
}



/* replace base ib by cr */
MDIIS.prototype.update = function(cr, res, uv)
{
  var i, j, l, id, ib, nb, mnb1, ns = this.ns, npt = this.npt;
  var dot, max;

  nb = this.nb;
  mnb1 = this.mnb + 1;

  if ( nb < this.mnb ) {
    ib = nb;
    this.nb = ++nb;
  } else {
    /* choose the base with the largest residue */
    ib = 0;
    for ( i = 1; i < nb; i++ )
      if ( this.mat[i*mnb1+i] > this.mat[ib*mnb1 + ib] )
        ib = i;
    max = this.mat[ib*mnb1 + ib];

    dot = getdot(res, res, uv.npr * npt);
    if ( dot > max ) {
      var reset = Math.sqrt(dot) < 1;
      if ( verbose )
        console.log("MDIIS: bad basis,", dot, ">", max, "reset", reset);
      if ( reset ) {
        this.build(cr, res, uv);
        return 1;
      }
    }
  }

  /* replace base ib by cr */
  for ( id = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv.prmask[i*ns + j] ) continue;
      for ( l = 0; l < npt; l++, id++ ) {
        this.cr[ib][id] = cr[i*ns + j][l];
        this.res[ib][id] = res[id];
      }
    }

  /* update the residue correlation matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ )
    this.mat[i*mnb1 + ib] = this.mat[ib*mnb1 + i]
      = getdot(this.res[i], res, uv.npr * npt);
  return ib;
}



function iter_mdiis(vrsr, wk, cr, ck, vklr, tr, tk, uv)
{
  var it, ibp = 0, ib;
  var err, errp = errinf, damp, res;

  var damp = get_float("mdiis_damp", 0.5);
  var nbases = get_int("mdiis_nbases", 5);

  /* open the mdiis object if needed */
  mdiis = new MDIIS(ns, npt, nbases);
  /* use the space of the last array for `res' */
  res = mdiis.res[mdiis.mnb];

  /* construct the initial base set */
  step_picard(res, null, vrsr, wk,
      cr, ck, vklr, tr, tk, uv.prmask, false, 0.);
  mdiis.build(cr, res, uv);

  for ( it = 0; it < itmax; it++ ) {
    mdiis.solve();
    mdiis.gencr(cr, damp, uv);
    err = step_picard(res, null, vrsr, wk,
        cr, ck, vklr, tr, tk, uv.prmask, false, 0.);
    ib = mdiis.update(cr, res, uv);

    if ( verbose )
      console.log("it", it, "err", errp, "->", err, "ib", ibp, "->", ib);
    if ( err < tol ) {
      if ( uv.switchstage() != 0 ) break;
      /* if all solutes are of zero density, break the loop */
      if ( uv.stage == SOLUTE_SOLUTE
        && uv.infdil && uv.atomicsolute ) {
        step_picard(res, null, vrsr, wk,
            cr, ck, vklr, tr, tk, uv.prmask, true, 1.);
        break;
      }
      /* reset the bases */
      step_picard(res, null, vrsr, wk,
          cr, ck, vklr, tr, tk, uv.prmask, false, 0.);
      mdiis_build(mdiis, cr, res, uv);
      it = -1;
      err = errinf;
    }
    ibp = ib;
    errp = err;
  }
  return [err, it];
}



