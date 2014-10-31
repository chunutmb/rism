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
  ns2 = ns * ns;
  m->npt = npt;
  m->mnb = mnb;
  m->nb = 0;
  mnb1 = mnb + 1;
  newarr2d(m->cr, mnb, ns2 * npt);
  newarr2d(m->res, mnb1, ns2 * npt);
  newarr(m->mat, mnb1 * mnb1);
  newarr(m->mat2, mnb1 * mnb1);
  newarr(m->coef, mnb1);
  return m;
}



/* close the mdiis object */
static void mdiis_close(mdiis_t *m)
{
  if ( m == NULL ) return;
  delarr2d(m->cr, m->mnb);
  delarr2d(m->res, m->mnb + 1);
  delarr(m->mat);
  delarr(m->mat2);
  delarr(m->coef);
  free(m);
}



/* compute residue from direct iteration (Picard) */
static double getres(model_t *model,
    double *res, double **fr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk)
{
  int ns = model->ns, npt = model->npt, i, j, ij, l;
  double y, err = 0;

  sphr_r2k(cr, ck, ns, NULL);
  oz(model, ck, vklr, tk, wk, NULL);
  sphr_k2r(tk, tr, ns, NULL);
  for ( i = 0; i < ns; i++ )
    for ( j = 0; j < ns; j++ )
      for ( ij = i*ns + j, l = 0; l < npt; l++ ) {
        y = getcr(tr[ij][l], fr[ij][l], NULL, model->ietype) - cr[ij][l];
        res[ij*npt + l] = y;
        if ( fabs(y) > err ) err = fabs(y);
      }
  return err;
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
static void mdiis_gencr(mdiis_t *m, double **cr, double damp)
{
  int ns = m->ns, npt = m->npt, nb = m->nb;
  int i, k, l, il;

  for ( i = 0; i < ns * ns; i++ )
    for ( l = 0; l < npt; l++ )
      cr[i][l] = 0;
  for ( k = 0; k < nb; k++ ) {
    double coef = m->coef[k];
    for ( il = 0, i = 0; i < ns * ns; i++ )
      for ( l = 0; l < npt; l++, il++ )
        cr[i][l] += coef * (m->cr[k][il] + damp * m->res[k][il]);
  }
}



/* compute the dot product */
static double getdot(double *a, double *b, int n)
{
  int i;
  double x = 0;

  for ( i = 0; i < n; i++ ) x += a[i] * b[i];
  return x;
}



/* build the matrix */
static int mdiis_build(mdiis_t *m, double **cr, double *res)
{
  int i, j, ib, id, mnb, mnb1, ns = m->ns, npt = m->npt, nps2;

  m->nb = 1;
  mnb = m->mnb;
  mnb1 = m->mnb + 1;
  nps2 = ns * ns * npt;

  for ( id = 0, j = 0; j < ns * ns; j++ )
    for ( i = 0; i < npt; i++, id++ ) {
      m->cr[0][id] = cr[j][i];
      m->res[0][id] = res[id];
    }

  m->mat[0] = getdot(m->res[0], m->res[0], nps2);
  for ( ib = 0; ib < mnb; ib++ )
    m->mat[ib*mnb1 + mnb] = m->mat[mnb*mnb1 + ib] = -1;
  m->mat[mnb*mnb1 + mnb] = 0;
  return 0;
}



/* replace base ib by cr */
static int mdiis_update(mdiis_t *m, double **cr, double *res)
{
  int i, j, id, ib, nb, mnb1, ns = m->ns, npt = m->npt, nps2;
  double dot, max;

  nb = m->nb;
  mnb1 = m->mnb + 1;
  nps2 = ns * ns * npt;

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

    dot = getdot(res, res, nps2);
    if ( dot > max ) {
      fprintf(stderr, "mdiis warning: %g is greater than %g, resetting\n", dot, max);
      mdiis_build(m, cr, res);
      return 1;
    }
  }

  /* replace base ib by cr */
  for ( id = 0, j = 0; j < ns * ns; j++)
    for ( i = 0; i < npt; i++, id++ ) {
      m->cr[ib][id] = cr[j][i];
      m->res[ib][id] = res[id];
    }

  /* update the matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ )
    m->mat[i*mnb1 + ib] = m->mat[ib*mnb1 + i]
      = getdot(m->res[i], res, nps2);
  return ib;
}



static double iter_mdiis(model_t *model,
    double **fr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, int *niter)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib;
  double err, errp = errinf, damp = model->mdiisdamp, *res;

  /* open the mdiis object if needed */
  mdiis = mdiis_open(model->ns, model->npt, model->nbases);
  /* use the space of the last array for `res' */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial base set */
  getres(model, res, fr, wk, cr, ck, vklr, tr, tk);
  mdiis_build(mdiis, cr, res);

  for ( it = 0; it < model->itmax; it++ ) {
    mdiis_solve(mdiis);
    mdiis_gencr(mdiis, cr, damp);
    err = getres(model, res, fr, wk, cr, ck, vklr, tr, tk);
    ib = mdiis_update(mdiis, cr, res);
    if ( err < model->tol ) break;
    if ( verbose )
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d\n", it, errp, err, ibp, ib);
    ibp = ib;
    errp = err;
  }
  *niter = it;
  mdiis_close(mdiis);
  return err;
}



#endif /* MDIIS_H__ */

