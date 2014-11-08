#ifndef MDIIS_H__
#define MDIIS_H__



/* modified direct inversion of the iterative subspace (MDIIS) method */



typedef struct {
  int ns;
  int npt;
  int mnb; /* maximal number of bases */
  int nb; /* number of bases */
  double **cr;  /* basic sets */
  double **res; /* residues */
  double *mat; /* correlations of residues */
  double *mat2; /* temporary matrix for LU decomposition */
  double *coef; /* coefficients */
} mdiis_t;



/* open an mdiis object */
static mdiis_t *mdiis_open(int ns, int npt, int mnb)
{
  mdiis_t *m;
  int ns2, mnb1;

  xnew(m, 1);
  m->ns = ns;
  ns2 = ns * (ns + 1) / 2;
  m->npt = npt;
  m->mnb = mnb;
  m->nb = 0;
  mnb1 = mnb + 1;
  newarr2d(m->cr,   mnb1, ns2 * npt);
  newarr2d(m->res,  mnb1, ns2 * npt);
  newarr(m->mat,    mnb1 * mnb1);
  newarr(m->mat2,   mnb1 * mnb1);
  newarr(m->coef,   mnb1);
  return m;
}



/* close the mdiis object */
static void mdiis_close(mdiis_t *m)
{
  if ( m == NULL ) return;
  delarr2d(m->cr,   m->mnb + 1);
  delarr2d(m->res,  m->mnb + 1);
  delarr(m->mat);
  delarr(m->mat2);
  delarr(m->coef);
  free(m);
}



/* solve the coefficients of combination */
static int mdiis_solve(mdiis_t *m)
{
  int nb = m->nb, nb1 = m->nb + 1, mnb1 = m->mnb + 1, i, j;

  for ( i = 0; i < nb; i++ ) m->coef[i] = 0;
  m->coef[nb] = -1;
  /* copy the matrix, for the content is to be destroyed */
  for ( i = 0; i < nb1; i++ )
    for ( j = 0; j < nb1; j++ )
      m->mat2[i*nb1 + j] = m->mat[i*mnb1 + j];
  for ( i = 0; i < nb1; i++ )
    m->mat2[i*nb1 + nb] = m->mat2[nb*nb1 + i] = -1;
  m->mat2[nb*nb1 + nb] = 0;
  i = lusolve(m->mat2, m->coef, nb1, 1e-20);
  return 0;
}



/* construct the new c(r) */
static void mdiis_gencr(mdiis_t *m, double **cr, double damp,
    const uv_t *uv)
{
  int ns = m->ns, npt = m->npt, nb = m->nb, npr = uv->npr;
  int i, j, ij, ipr, k, l, il;

  for ( il = 0; il < npr * npt; il++ )
    m->cr[nb][il] = 0;
  for ( k = 0; k < nb; k++ ) {
    double coef = m->coef[k];
    for ( il = 0; il < npr * npt; il++ )
      m->cr[nb][il] += coef * (m->cr[k][il] + damp * m->res[k][il]);
  }
  /* save m->cr[nb] to c(r) */
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv->prmask[ij = i*ns + j] ) continue;
      for ( l = 0; l < npt; l++ ) {
        double y;
        cr[ij][l] = y = m->cr[nb][ipr*npt + l];
        if ( j > i) cr[j*ns + i][l] = y;
      }
      ipr++;
    }
}



/* compute the dot product */
static double getdot(double *a, double *b, int n)
{
  int i;
  double x = 0;

  for ( i = 0; i < n; i++ ) x += a[i] * b[i];
  return x / n;
}



/* build the residue correlation matrix */
static int mdiis_build(mdiis_t *m, double **cr, double *res,
    const uv_t *uv)
{
  int i, j, ipr, l, ib, id, mnb, mnb1, ns = m->ns, npt = m->npt;

  m->nb = 1;
  mnb = m->mnb;
  mnb1 = m->mnb + 1;

  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv->prmask[i*ns + j] ) continue;
      for ( l = 0; l < npt; l++ ) {
        id = ipr*npt + l;
        m->cr[0][id] = cr[i*ns + j][l];
        m->res[0][id] = res[id];
      }
      ipr++;
    }

  m->mat[0] = getdot(m->res[0], m->res[0], uv->npr * npt);
  for ( ib = 0; ib < mnb; ib++ )
    m->mat[ib*mnb1 + mnb] = m->mat[mnb*mnb1 + ib] = -1;
  m->mat[mnb*mnb1 + mnb] = 0;
  return 0;
}



/* replace base ib by cr */
static int mdiis_update(mdiis_t *m, double **cr, double *res,
    const uv_t *uv)
{
  int i, j, l, id, ib, nb, mnb1, ns = m->ns, npt = m->npt;
  double dot, max;

  nb = m->nb;
  mnb1 = m->mnb + 1;

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  } else {
    /* choose the base with the largest residue */
    ib = 0;
    for ( i = 1; i < nb; i++ )
      if ( m->mat[i*mnb1+i] > m->mat[ib*mnb1 + ib] )
        ib = i;
    max = m->mat[ib*mnb1 + ib];

    dot = getdot(res, res, uv->npr * npt);
    if ( dot > max ) {
      int reset = sqrt(dot) < 1;
      if ( verbose )
        fprintf(stderr, "MDIIS: bad basis, %g is greater than %g%s\n",
          dot, max, reset ? ", reset" : "");
      if ( reset ) {
        mdiis_build(m, cr, res, uv);
        return 1;
      }
    }
  }

  /* replace base ib by cr */
  for ( id = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !uv->prmask[i*ns + j] ) continue;
      for ( l = 0; l < npt; l++, id++ ) {
        m->cr[ib][id] = cr[i*ns + j][l];
        m->res[ib][id] = res[id];
      }
    }

  /* update the residue correlation matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ )
    m->mat[i*mnb1 + ib] = m->mat[ib*mnb1 + i]
      = getdot(m->res[i], res, uv->npr * npt);
  return ib;
}



static double iter_mdiis(model_t *model,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, int *niter)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib, ns = model->ns, npt = model->npt;
  double err, errp = errinf, damp = model->mdiis.damp, *res;
  uv_t *uv;

  /* initialize the manager for solvent-solvent iteraction */
  uv = uv_open(model);

  /* open the mdiis object if needed */
  mdiis = mdiis_open(ns, npt, model->mdiis.nbases);
  /* use the space of the last array for `res' */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial base set */
  step_picard(model, res, NULL, vrsr, wk,
      cr, ck, vklr, tr, tk, uv->prmask, 0, 0.);
  mdiis_build(mdiis, cr, res, uv);

  for ( it = 0; it < model->itmax; it++ ) {
    mdiis_solve(mdiis);
    mdiis_gencr(mdiis, cr, damp, uv);
    err = step_picard(model, res, NULL, vrsr, wk,
        cr, ck, vklr, tr, tk, uv->prmask, 0, 0.);
    ib = mdiis_update(mdiis, cr, res, uv);

    if ( verbose )
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d\n",
          it, errp, err, ibp, ib);
    if ( err < model->tol ) {
      if ( uv_switch(uv) != 0 ) break;
      /* if all solutes are of zero density, break the loop */
      if ( uv->stage == SOLUTE_SOLUTE && uv->infdil ) {
        step_picard(model, res, NULL, vrsr, wk,
            cr, ck, vklr, tr, tk, uv->prmask, 1, 1.);
        break;
      }
      /* reset the bases */
      step_picard(model, res, NULL, vrsr, wk,
          cr, ck, vklr, tr, tk, uv->prmask, 0, 0.);
      mdiis_build(mdiis, cr, res, uv);
      it = -1;
      err = errinf;
    }
    ibp = ib;
    errp = err;
  }
  *niter = it;
  mdiis_close(mdiis);

  uv_close(uv);
  return err;
}



#endif /* MDIIS_H__ */

