#ifndef LMV_H__
#define LMV_H__



/* Labik-Malijevsky-Vonka solver */



/* compute Cjk */
static void getCjk(double **Cjk, int npt, int M, int ns,
    double **der, double **costab)
{
  int i, j, ij, ji, m, k, l;
  double y, *dp;

  newarr(dp, 3*M);
  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;

      for ( m = 1; m < 3*M - 1; m++ ) {
        for ( y = 0, l = 0; l < npt; l++ )
          y += der[ij][l] * costab[m][l];
        dp[m] = y / npt;
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
  delarr(dp);
}



/* compute the Jacobian matrix for the Newton-Raphson method */
static void getjacob(double *mat, double *b, int M, int npr, int ns,
    int *prmask, double **ntk, double **tk, double **Cjk, double **invwc1w)
{
  int i1, j1, ipr1, m1, id1, i2, j2, ipr2, m2, id2, Mp;
  double y;

  Mp = M * npr;
  for ( m1 = 0; m1 < M; m1++ ) {
    for ( ipr1 = 0, i1 = 0; i1 < ns; i1++ ) {
      for ( j1 = i1; j1 < ns; j1++, ipr1++ ) {
        id1 = m1*npr + ipr1;
        if ( prmask[i1*ns + j1] )
          b[id1] = fft_ki[m1] * (ntk[i1*ns+j1][m1] - tk[i1*ns+j1][m1]);
        else
          b[id1] = 0;

        for ( m2 = 0; m2 < M; m2++ ) {
          for ( ipr2 = 0, i2 = 0; i2 < ns; i2++ ) {
            for ( j2 = i2; j2 < ns; j2++, ipr2++ ) {
              id2 = m2*npr + ipr2;
              //if (id1 >= Mp || id2 >= Mp) {
              //  fprintf(stderr, "id1 %d, id2 %d, m2 %d, ipr2 %d/%d\n", id1, id2, m2, ipr2, npr);
              //  exit(1);
              //}
              y = (ipr1 == ipr2 ? (m1 == m2) + Cjk[i1*ns+j1][m1*M+m2] : 0)
                - invwc1w[i1*ns+i2][m1] * Cjk[i2*ns+j2][m1*M+m2]
                * invwc1w[j2*ns+j1][m2];
              mat[id1*Mp + id2] = y;
            }
          }
        }
      }
    }
  }
}



/* Reference:
 * Stanislav Labik, Anatol Malijevsky, Petr Vonka
 * A rapidly convergent method of solving the OZ equation
 * Molecular Physics, 1985, Vol. 56, No. 3, 709-715 */
static double iter_lmv(model_t *model,
    double **fr, double **wk,
    double **cr, double **der, double **ck, double **vklr,
    double **tr, double **tk, double **ntk,
    double **invwc1w, int *niter)
{
  int i, j, l, ij, it, M, npr, ipr, Mp;
  int ns = model->ns, npt = model->npt;
  double **Cjk = NULL, *mat = NULL, *b = NULL, **costab = NULL;
  double y, err1 = 0, err2 = 0, err = 0, errp = errinf, dmp;
  uv_t *uv;

  /* initialize t(k) and t(r) */
  sphr_r2k(cr, ck, ns, NULL);
  oz(model, ck, vklr, tk, wk, NULL); /* c(k) --> t(k) */
  sphr_k2r(tk, tr, ns, NULL);

  /* set the optimal M */
  M = (model->Mpt > 0) ? model->Mpt :
       (int) (2 * model->rmax/model->sigma[ns-1]);
  if ( M >= npt ) M = npt;
  if ( verbose ) fprintf(stderr, "select M = %d\n", M);

  /* set the damping factor */
  dmp = (model->lmvdamp > 0) ? model->lmvdamp : 1;

  npr = ns * (ns + 1) / 2;
  Mp = M * npr;

  /* initialize the cosine table */
  if ( M > 0 ) {
    newarr2d(Cjk, ns*ns, M*M);
    newarr(mat, Mp*Mp);
    newarr(b, Mp);
    newarr2d(costab, 3*M, npt);
    for ( j = 0; j < 3*M; j++ )
      for ( i = 0; i < npt; i++ )
        costab[j][i] = cos(PI*(i+.5)*(j-M)/npt);
  }

  /* initialize the manager for solvent-solvent iteraction */
  uv = uv_open(model);

  for ( it = 0; it < model->itmax; it++ ) {
    /* compute c(r) and c(k) from the closure */
    err = 0;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        for ( ij = i*ns + j, l = 0; l < npt; l++ ) {
          y = getcr(tr[ij][l], fr[ij][l], &der[ij][l], model->ietype) - cr[ij][l];
          if ( fabs(y) > err ) err = fabs(y);
          cr[ij][l] += y;
        }
    sphr_r2k(cr, ck, ns, NULL);

    /* compute Cjk */
    getCjk(Cjk, npt, M, ns, der, costab);

    /* compute the new t(k) */
    oz(model, ck, vklr, ntk, wk, invwc1w);

    /* compute the Jacobian matrix for the Newton-Raphson method */
    getjacob(mat, b, M, npr, ns, uv->prmask, ntk, tk, Cjk, invwc1w);

    if ( lusolve(mat, b, Mp, DBL_MIN) != 0 ) break;

    /* compute the new t(k) */
    err1 = err2 = 0;
    for ( ipr = 0, i = 0; i < ns; i++ )
      for ( j = i; j < ns; j++, ipr++ ) {
        ij = i*ns + j;
        if ( !uv->prmask[ij] ) continue;

        for ( l = 0; l < npt; l++ ) {
          if ( l < M ) {
            /* use the Newton-Raphson method to solve for t(k) of small k */
            y = b[l*npr+ipr] / fft_ki[l];
            if ( fabs(y) > err1 ) err1 = fabs(y);
          } else {
            /* use the OZ relation to solve for t(k) of large k */
            y = ntk[ij][l] - tk[ij][l];
            if ( fabs(y) > err2 ) err2 = fabs(y);
          }

          tk[ij][l] += dmp * y;
          if ( j > i ) tk[j*ns+i][l] = tk[ij][l];
        }
      }
    /* note that error of t(k) is different from that of c(r) in scale */
    //err = (err1 > err2) ? err1 : err2;

    sphr_k2r(tk, tr, ns, NULL);

    if ( verbose )
      fprintf(stderr, "it %d: M %d, err %g -> %g, tk_err %g/%g, damp %g\n",
          it, M, errp, err, err1, err2, dmp);

    if ( it > 0 && err < model->tol ) {
      /* switch between stages */
      int brk = uv_switch(uv);
      if ( uv->stage == 2 ) {
        oz(model, ck, vklr, tk, wk, invwc1w);
        sphr_k2r(tk, tr, ns, NULL);
        /* NOTE currently the we cannot continue from here
         * thus, we stop at the case of infinite dilution */
        //fprintf(stderr, "turning on solute-solute interaction\n"); getchar();
      }
      if ( brk ) break;
      it = -1;
      err = errinf;
    }
/*
    // adaptively adjust the damping factor
    if ( err > errp ) {
      dmp *= 0.8;
    } else {
      dmp = dmp * 0.9 + 0.1;
    }
*/
    errp = err;
  }
  *niter = it;
  if ( M > 0 ) {
    delarr2d(Cjk, ns*ns);
    delarr(mat);
    delarr(b);
    delarr2d(costab, 3*M);
  }
  uv_close(uv);
  return err;
}



#endif /* LMV_H__ */

