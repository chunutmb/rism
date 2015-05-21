/* modified direct inversion of the iterative subspace (MDIIS) method */



"use strict";



function MDIIS(ns, npt, mnb)
{
  var ns2 = ns * (ns + 1) / 2;
  this.ns    = ns;
  this.npt   = npt;
  this.mnb   = mnb; // maximal number of bases
  var mnb1 = mnb + 1;
  this.cr    = newarr2d(mnb1, ns2 * npt); // basis
  this.res   = newarr2d(mnb1, ns2 * npt); // residues
  this.mat   = newarr(mnb * mnb); // correlations of residues
  this.mat2  = newarr(mnb1 * mnb1); // temporary matrix for LU decomposition
  this.coef  = newarr(mnb1); // coefficients
  this.crbest = newarr(ns2 * npt); // best function
  this.errmin = errinf;
}



/* solve the coefficients of combination */
MDIIS.prototype.solve = function()
{
  var nb = this.nb, nb1 = this.nb + 1, mnb = this.mnb, i, j;

  for ( i = 0; i < nb; i++ ) {
    this.coef[i] = 0;
  }
  this.coef[nb] = -1;
  // copy the matrix, for the content is to be destroyed
  for ( i = 0; i < nb1; i++ ) {
    for ( j = 0; j < nb1; j++ ) {
      this.mat2[i*nb1 + j] = this.mat[i*mnb + j];
    }
  }
  for ( i = 0; i < nb1; i++ ) {
    this.mat2[i*nb1 + nb] = this.mat2[nb*nb1 + i] = -1;
  }
  this.mat2[nb*nb1 + nb] = 0;
  i = lusolve(this.mat2, this.coef, nb1, 1e-20);
  return 0;
};



/* save src in cr */
MDIIS.prototype.copyout = function(cr, src, uv)
{
  var ns = this.ns, npt = this.npt;
  var i, j, ij, ipr, l;

  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv.prmask[ij = i*ns + j] ) continue;
      for ( l = 0; l < npt; l++ )
        cr[ij][l] = src[ipr*npt + l];
      ipr++;
      cparr(cr[j*ns + i], cr[ij], npt);
    }
};



/* construct the new c(r) */
MDIIS.prototype.gencr = function(cr, damp, uv)
{
  var ns = this.ns, npt = this.npt, nb = this.nb, npr = uv.npr;
  var i, j, ij, ipr, k, l, il, coef;

  for ( il = 0; il < npr * npt; il++ )
    this.cr[nb][il] = 0;
  for ( k = 0; k < nb; k++ ) {
    coef = this.coef[k];
    for ( il = 0; il < npr * npt; il++ )
      this.cr[nb][il] += coef * (this.cr[k][il] + damp * this.res[k][il]);
  }
  // cr = this.cr[nb]
  this.copyout(cr, this.cr[nb], uv);
};



/* compute the dot product */
MDIIS.prototype.getdot = function(a, b, n)
{
  var i, x = 0;

  for ( i = 0; i < n; i++ ) x += a[i] * b[i];
  return x / n;
};



/* build the residue correlation matrix */
MDIIS.prototype.build = function(cr, res, uv)
{
  var i, j, ipr, l, id, ns = this.ns, npt = this.npt;

  this.nb = 1;

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

  this.mat[0] = this.getdot(this.res[0], this.res[0], uv.npr * npt);
  return 0;
};




/* replace base ib by cr */
MDIIS.prototype.update = function(cr, res, err, uv)
{
  var i, j, l, id, ib, nb, mnb, ns = this.ns, npt = this.npt;
  var dot, max;

  nb = this.nb;
  mnb = this.mnb;

  // save this function if it achieves the minimal error so far
  if ( err < this.errmin ) {
    cparr(this.crbest, this.cr[nb], uv.npr * npt);
    this.errmin = err;
  }

  // choose the base with the largest residue
  ib = 0;
  for ( i = 1; i < nb; i++ )
    if ( this.mat[i*mnb+i] > this.mat[ib*mnb + ib] )
      ib = i;
  max = this.mat[ib*mnb + ib];

  dot = this.getdot(res, res, uv.npr * npt);
  if ( dot > max ) {
    if ( nb >= 2 ) {
      // remove the base with the largest residue
      var jb = nb - 1;
      if ( ib != jb ) {
        // move the last base to position ib
        for ( l = 0, i = 0; i < ns; i++ ) {
          for ( j = i; j < ns; j++ ) {
            if ( !uv.prmask[i*ns + j] ) continue;
            l += npt;
          }
        }

        for ( id = 0; id < l; id++ ) {
          this.cr[ib][id] = this.cr[jb][id];
          this.res[ib][id] = this.res[jb][id];
        }

        for ( i = 0; i < nb - 1; i++ ) {
          if ( i == ib ) continue;
          this.mat[i*mnb + ib] = this.mat[i*mnb + jb];
          this.mat[ib*mnb + i] = this.mat[jb*mnb + i];
        }
        this.mat[ib*mnb + ib] = this.mat[jb*mnb + jb];
      }
      this.nb--;
      return -1;
    } else {
      this.build(cr, res, uv);
      return 0;
    }
  }

  if ( nb < mnb ) {
    ib = nb;
    this.nb = ++nb;
  }

  // replace base ib by cr
  for ( id = 0, i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      if ( !uv.prmask[i*ns + j] ) continue;
      for ( l = 0; l < npt; l++, id++ ) {
        this.cr[ib][id] = cr[i*ns + j][l];
        this.res[ib][id] = res[id];
      }
    }
  }

  // update the residue correlation matrix
  // note: we do not need to update the last row & column
  for ( i = 0; i < nb; i++ ) {
    this.mat[i*mnb + ib] = this.mat[ib*mnb + i]
      = this.getdot(this.res[i], res, uv.npr * npt);
  }
  return ib;
};



function iter_mdiis(vrsr, wk, cr, ck, vklr, tr, tk, uv)
{
  var it, ibp = 0, ib;
  var err, errp = errinf, res;

  var damp = get_float("mdiis_damp", 0.5);
  var nbases = get_int("mdiis_nbases", 5);

  // open an mdiis object
  var mdiis = new MDIIS(ns, npt, nbases);
  // use the space of the last array for `res'
  res = mdiis.res[mdiis.mnb];

  // construct the initial base set
  step_picard(res, vrsr, wk, cr, ck, vklr, tr, tk, uv.prmask, 0.);
  mdiis.build(cr, res, uv);

  for ( it = 0; it <= itmax; it++ ) {
    mdiis.solve();
    mdiis.gencr(cr, damp, uv);
    err = step_picard(res, vrsr, wk, cr, ck, vklr, tr, tk, uv.prmask, 0.);
    ib = mdiis.update(cr, res, err, uv);

    if ( verbose ) {
      console.log("it", it, "err", errp, "->", err, "ib", ibp, "->", ib);
    }

    if ( err < tol || it == itmax ) {
      mdiis.copyout(cr, mdiis.crbest, uv);
      err = step_picard(res, vrsr, wk, cr, ck, vklr, tr, tk, uv.prmask, 0.);
      if ( uv.switchstage() !== 0 ) {
        if ( uv.uu1step )
          step_picard(res, vrsr, wk, cr, ck, vklr, tr, tk, uv.prmask, 1.);
        break;
      }
      // reset the basis
      err = step_picard(res, vrsr, wk, cr, ck, vklr, tr, tk, uv.prmask, 0.);
      mdiis.build(cr, res, uv);
      mdiis.errmin = err;
      it = -1;
    }

    if ( ib < 0 ) {
      continue;
    }

    ibp = ib;
    errp = err;
  }
  return [err, it];
};



