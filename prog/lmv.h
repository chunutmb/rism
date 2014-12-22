#ifndef LMV_H__
#define LMV_H__



/* Labik-Malijevsky-Vonka (LMV) solver
 * Reference:
 * Stanislav Labik, Anatol Malijevsky, Petr Vonka
 * A rapidly convergent method of solving the OZ equation
 * Molecular Physics, 1985, Vol. 56, No. 3, 709-715 */



typedef struct {
  int ns;
  int npr;
  int npt;
  int M;
  double *ki;
  double **Cjk; /* d ck / d tk */
  double *mat; /* M x M matrix */
  double *a; /* M array */
  double *dp; /* 3*M array */
  double **costab; /* [3*M][npt] */
  double **tr1;
  double **tk1;
  double **invwc1w;
  double **der;
  double **crbest;
  double errmin;
  double err1; /* error of tk[i] for i < M */
  double err2; /* error of tk[i] for i >= M */
} lmv_t;




/* open an lmv object */
static lmv_t *lmv_open(int ns, int npt, int M, double *ki)
{
  lmv_t *lmv;
  int i, j, ns2, Mp;

  xnew(lmv, 1);
  lmv->ns = ns;
  ns2 = ns * ns;
  lmv->npr = ns * (ns + 1)/2;
  lmv->npt = npt;
  if ( M >= npt ) M = npt;
  if ( verbose ) fprintf(stderr, "select M = %d\n", M);
  lmv->M = M;
  lmv->ki = ki;
  newarr2d(lmv->tr1,     ns2, npt);
  newarr2d(lmv->tk1,     ns2, npt);
  newarr2d(lmv->invwc1w, ns2, npt);
  newarr2d(lmv->der,     ns2, npt);
  newarr2d(lmv->crbest,  ns2, npt);
  lmv->errmin = errinf;

  if ( M > 0 ) {
    Mp = M * lmv->npr;
    newarr2d(lmv->Cjk, ns2, M * M);
    newarr(lmv->mat, Mp * Mp);
    newarr(lmv->a,  Mp);
    newarr(lmv->dp, 3*M);
    newarr2d(lmv->costab, 3*M, npt);
    for ( j = 0; j < 3*M; j++ )
      for ( i = 0; i < npt; i++ )
        lmv->costab[j][i] = cos(PI*(i*2+1)*(j-M)/npt/2);
  }

  return lmv;
}



static void lmv_close(lmv_t *lmv)
{
  int ns2;

  if ( lmv == NULL ) return;
  ns2 = lmv->ns * lmv->ns;
  delarr2d(lmv->tr1,     ns2);
  delarr2d(lmv->tk1,     ns2);
  delarr2d(lmv->invwc1w, ns2);
  delarr2d(lmv->der,     ns2);
  delarr2d(lmv->crbest,  ns2);
  if ( lmv->M > 0 ) {
    delarr2d(lmv->Cjk, ns2);
    delarr(lmv->mat);
    delarr(lmv->a);
    delarr(lmv->dp);
    delarr(lmv->costab);
  }
  free(lmv);
}



/* register a good cr */
static void lmv_savebest(lmv_t *lmv, double **cr, double err)
{
  if ( err < lmv->errmin ) {
    cparr2d(lmv->crbest, cr, lmv->ns * lmv->ns, lmv->npt);
    lmv->errmin = err;
  }
}



/* compute Cjk = d ck / d tk */
static void lmv_getCjk(lmv_t *lmv)
{
  int ns = lmv->ns, npt = lmv->npt, M = lmv->M;
  int i, j, ij, ji, m, k, l;

  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;

      for ( m = 1; m < 3*M - 1; m++ ) {
        for ( lmv->dp[m] = 0, l = 0; l < npt; l++ )
          lmv->dp[m] += lmv->der[ij][l] * lmv->costab[m][l];
        lmv->dp[m] /= npt;
      }

      for ( m = 0; m < M; m++ )
        for ( k = 0; k < M; k++ )
          lmv->Cjk[ij][m*M+k] = lmv->dp[k-m+M] - lmv->dp[k+m+M];

      if ( j == i ) continue;

      ji = j*ns + i;
      for ( m = 0; m < M*M; m++ )
        lmv->Cjk[ji][m] = lmv->Cjk[ij][m];
    }
  }
}



/* compute the Jacobian matrix for the Newton-Raphson method */
static void lmv_getjacob(lmv_t *lmv, double **tk, uv_t *uv)
{
  int ns = lmv->ns, npr = lmv->npr, M = lmv->M, Mp;
  int i1, j1, ij1, ipr1, m1, id1, i2, j2, ipr2, m2, id2;

  Mp = M * npr;
  for ( m1 = 0; m1 < M; m1++ ) {
    for ( ipr1 = 0, i1 = 0; i1 < ns; i1++ ) {
      for ( j1 = i1; j1 < ns; j1++, ipr1++ ) {
        ij1 = i1*ns + j1;
        id1 = m1*npr + ipr1;
        lmv->a[id1] = uv->prmask[ij1] ? lmv->ki[m1] * (lmv->tk1[ij1][m1] - tk[ij1][m1]) : 0;

        for ( m2 = 0; m2 < M; m2++ ) {
          for ( ipr2 = 0, i2 = 0; i2 < ns; i2++ ) {
            for ( j2 = i2; j2 < ns; j2++, ipr2++ ) {
              id2 = m2*npr + ipr2;
              //if (id1 >= Mp || id2 >= Mp) {
              //  fprintf(stderr, "id1 %d, id2 %d, m2 %d, ipr2 %d/%d\n", id1, id2, m2, ipr2, npr);
              //  exit(1);
              //}
              /*
               * h = w c (1 - rho w c)^-1 w
               *   = w c w + w c rho w c w + w c rho w c rho w c w + ...
               * thus
               * dh = (1 - w c rho)^(-1) w dc (1 - rho w c)^(-1) w
               * let z = (1 - rho w c)^(-1) w
               * then
               * z^T = w (1 - c w rho)^(-1)
               *     = w (1 - c rho w)^(-1)   (since w rho = rho w)
               *     = w + w c rho w + w c rho w c rho w + ...
               *     = (1 - w c rho)^-1 w
               * and
               * dh = z^T dc z */
              lmv->mat[id1*Mp + id2] =
                (ipr1 == ipr2 ? (m1 == m2) + lmv->Cjk[ij1][m1*M+m2] : 0)
                - lmv->invwc1w[i2*ns+i1][m1] * lmv->Cjk[i2*ns+j2][m1*M+m2]
                * lmv->invwc1w[j2*ns+j1][m2];
            }
          }
        }
      }
    }
  }
}



/* update tk */
static int lmv_update(lmv_t *lmv, double **tk, double dmp, uv_t *uv)
{
  int ns = lmv->ns, npr = lmv->npr, npt = lmv->npt, M = lmv->M;
  int i, j, ij, l, ipr;
  double del;

  /* compute d ck / d tk */
  lmv_getCjk(lmv);

  /* compute the Jacobian matrix for the Newton-Raphson method */
  lmv_getjacob(lmv, tk, uv);

  if ( lusolve(lmv->mat, lmv->a, npr * M, 1e-10) != 0 ) {
    fprintf(stderr, "LU solve failed: stage %d\n", uv->stage);
    return -1;
  }

  /* compute the new t(k) */
  lmv->err1 = lmv->err2 = 0;
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++, ipr++ ) {
      ij = i*ns + j;
      if ( !uv->prmask[ij] ) continue;

      for ( l = 0; l < npt; l++ ) {
        if ( l < M ) {
          /* use the Newton-Raphson method to solve for t(k) of small k */
          del = lmv->a[l*npr+ipr] / lmv->ki[l];
          if ( fabs(del) > lmv->err1 ) lmv->err1 = fabs(del);
        } else {
          /* use the OZ relation to solve for t(k) of large k */
          del = lmv->tk1[ij][l] - tk[ij][l];
          if ( fabs(del) > lmv->err2 ) lmv->err2 = fabs(del);
        }

        tk[ij][l] += dmp * del;
      }
      if ( j > i ) cparr(tk[j*ns + i], tk[ij], npt);
    }

  return 0;
}



/* LMV solver */
static double iter_lmv(model_t *model,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, double **Qrx,
    uv_t *uv, int *niter)
{
  int it, M, ns = model->ns, npt = model->npt;
  double err = 0, errp = errinf, dmp;
  lmv_t *lmv;

  /* set the optimal M */
  M = (model->lmv.M > 0) ? model->lmv.M :
       (int) (2 * model->rmax/model->sigma[ns-1]);

  /* open an lmv object */
  lmv = lmv_open(ns, npt, model->lmv.M, fft_ki);
  cparr2d(lmv->crbest, cr, ns * ns, npt);

  /* initialize t(k) and t(r) */
  step_picard(model, NULL, vrsr, wk, cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
  cparr2d(lmv->tk1, tk, ns * ns, npt);

  /* set the damping factor */
  dmp = (model->lmv.damp > 0) ? model->lmv.damp : 1;

  for ( it = 0; it <= model->itmax; it++ ) {
    /* compute the error of the current c(r) and c(k) */
    sphr_k2r(lmv->tk1, lmv->tr1, ns, uv->prmask);
    err = closure(model, NULL, NULL, vrsr, cr, lmv->tr1, Qrx, uv->prmask, 0.);
    lmv_savebest(lmv, cr, err);

    /* compute c(r) and c(k) from the closure */
    closure(model, NULL, lmv->der, vrsr, cr, tr, Qrx, uv->prmask, 1.);
    sphr_r2k(cr, ck, ns, uv->prmask);
    oz(model, ck, vklr, lmv->tk1, wk, lmv->invwc1w);

    /* compute the new tk */
    if ( lmv_update(lmv, tk, dmp, uv) != 0 ) break;
    sphr_k2r(tk, tr, ns, uv->prmask);

    if ( verbose )
      fprintf(stderr, "it %d: M %d, err %g -> %g, tk_err %g/%g, damp %g\n",
          it, M, errp, err, lmv->err1, lmv->err2, dmp);

    if ( err < model->tol || it == model->itmax ) {
      /* use the best cr discovered so far */
      cparr2d(cr, lmv->crbest, ns * ns, npt);
      /* update the corresponding ck, tr, tk, and the error */
      err = step_picard(model, NULL, vrsr, wk,
          cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
      //fprintf(stderr, "switching stage %d, it %d, err %g, tol %g\n", uv->stage, it, err, model->tol); getchar();
      /* switch between stages */
      if ( uv_switch(uv) != 0 ) {
        if ( uv->uu1step )
          err = step_uu_infdil_atomicsolute(model, vrsr, wk,
              cr, ck, vklr, tr, tk, Qrx, uv->prmask);
        break;
      }
      err = step_picard(model, NULL, vrsr, wk,
          cr, ck, vklr, tr, tk, Qrx, uv->prmask, 0.);
      cparr2d(lmv->tk1, tk, ns * ns, npt);
      lmv->errmin = err;
      it = -1;
    }
    errp = err;
  }
  *niter = it;
  lmv_close(lmv);
  return err;
}



#endif /* LMV_H__ */

