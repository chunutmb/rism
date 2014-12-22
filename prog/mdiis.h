#ifndef MDIIS_H__
#define MDIIS_H__



/* modified direct inversion of the iterative subspace (MDIIS) method */



typedef struct {
  int ns;
  int npt;
  int mnb; /* maximal number of bases */
  int nb; /* number of functions in the basis */
  double **cr;  /* basis */
  double **res; /* residues */
  double *mat; /* correlations of residues */
  double *mat2; /* temporary matrix for LU decomposition */
  double *coef; /* coefficients */
  double *crbest;
  double errmin;
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
  newarr(m->crbest, ns2 * npt);
  m->errmin = errinf;
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
  delarr(m->crbest);
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
  if ( lusolve(m->mat2, m->coef, nb1, 1e-20) != 0 ) {
    fprintf(stderr, "MDIIS lusolve failed\n");
    exit(1);
  }
  return 0;
}



/* copy out cr */
static void mdiis_copyout(double **cr, double *src,
    int ns, int npt, int *prmask)
{
  int i, j, ij, ipr, l;

  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      if ( !prmask[ij = i*ns + j] ) continue;
      for ( l = 0; l < npt; l++ )
        cr[ij][l] = src[ipr*npt + l];
      ipr++;
      cparr(cr[j*ns + i], cr[ij], npt);
    }
}



/* construct the new c(r) */
static void mdiis_gencr(mdiis_t *m, double **cr, double damp,
    const uv_t *uv)
{
  int ib, il, npt = m->npt, nb = m->nb, npr = uv->npr;

  for ( il = 0; il < npr * npt; il++ )
    m->cr[nb][il] = 0;
  for ( ib = 0; ib < nb; ib++ ) {
    double coef = m->coef[ib];
    for ( il = 0; il < npr * npt; il++ )
      m->cr[nb][il] += coef * (m->cr[ib][il] + damp * m->res[ib][il]);
  }

  /* cr = m->cr[nb] */
  mdiis_copyout(cr, m->cr[nb], m->ns, m->npt, uv->prmask);
}



/* compute the dot product */
static double mdiis_getdot(double *a, double *b, int n)
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

  m->mat[0] = mdiis_getdot(m->res[0], m->res[0], uv->npr * npt);
  for ( ib = 0; ib < mnb; ib++ )
    m->mat[ib*mnb1 + mnb] = m->mat[mnb*mnb1 + ib] = -1;
  m->mat[mnb*mnb1 + mnb] = 0;
  return 0;
}



/* replace base ib by cr */
static int mdiis_update(mdiis_t *m, double **cr, double *res,
    double err, const uv_t *uv)
{
  int i, j, l, id, ib, nb, mnb1, ns = m->ns, npt = m->npt;
  double dot, max;

  nb = m->nb;
  mnb1 = m->mnb + 1;

  /* save this function if it achieves the minimal error so far */
  if ( err < m->errmin ) {
    cparr(m->crbest, m->cr[nb], uv->npr * npt);
    m->errmin = err;
  }

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  } else {
    /* choose the base with the largest residue */
    ib = 0;
    for ( i = 1; i < nb; i++ )
      /* the diagonal represents the error */
      if ( m->mat[i*mnb1+i] > m->mat[ib*mnb1 + ib] )
        ib = i;
    max = m->mat[ib*mnb1 + ib];

    dot = mdiis_getdot(res, res, uv->npr * npt);
    if ( dot > max ) {
#ifndef MDIIS_THRESHOLD
#define MDIIS_THRESHOLD 1.0
#endif
      int reset = ( sqrt(dot) < MDIIS_THRESHOLD );
      if ( verbose ) {
        fprintf(stderr, "MDIIS: bad basis, %g is greater than %g, %s, error:",
          dot, max, reset ? "reset" : "accept");
        for ( i = 0; i < nb; i++ )
          fprintf(stderr, " %g", m->mat[i*mnb1+i]);
        fprintf(stderr, "\n");
      }
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
      = mdiis_getdot(m->res[i], res, uv->npr * npt);
  return ib;
}



static double iter_mdiis(model_t *model,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, double **Qrx,
    uv_t *uv, int *niter)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib, ns = model->ns, npt = model->npt;
  double err, errp, damp = model->mdiis.damp, *res;

  /* open the mdiis object if needed */
  mdiis = mdiis_open(ns, npt, model->mdiis.nbases);
  /* use the space of the last array for `res' */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial base set */
  step_picard(model, res, vrsr, wk,
      cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
  mdiis_build(mdiis, cr, res, uv);

  err = errp = errinf;
  for ( it = 0; it <= model->itmax; it++ ) {
    /* obtain a set of optimal coefficients of combination */
    mdiis_solve(mdiis);
    /* generate a new cr from the set of coefficients */
    mdiis_gencr(mdiis, cr, damp, uv);
    err = step_picard(model, res, vrsr, wk,
        cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
    /* add the new cr into the basis */
    ib = mdiis_update(mdiis, cr, res, err, uv);

    if ( verbose )
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d\n",
          it, errp, err, ibp, ib);
    if ( err < model->tol || it == model->itmax ) {
      mdiis_copyout(cr, mdiis->crbest, ns, npt, uv->prmask);
      err = step_picard(model, res, vrsr, wk,
          cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
      if ( it >= model->itmax && verbose )
        fprintf(stderr, "MDIIS: failure, stage %d, use the best c(r) with residue %g\n",
           uv->stage, mdiis->errmin); //getchar();
      if ( uv_switch(uv) != 0 ) {
        if ( uv->uu1step )
          err = step_uu_infdil_atomicsolute(model, vrsr, wk,
              cr, ck, vklr, tr, tk, Qrx, uv->prmask);
        break;
      }
      /* reset the basis */
      err = step_picard(model, res, vrsr, wk,
          cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
      mdiis_build(mdiis, cr, res, uv);
      mdiis->errmin = err;
      it = -1;
    }
    ibp = ib;
    errp = err;
  }
  *niter = it;
  mdiis_close(mdiis);
  return err;
}



#endif /* MDIIS_H__ */

