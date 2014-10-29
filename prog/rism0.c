/* basic RISM, assume a single solvent molecule
 *
 *  gcc rism0.c -lfftw3 -lm
 *
 * the output can be found in "out.dat"
 * to view g(r) in gnuplot, type
 *
 *  plot "out.dat" u 1:($2+$3+1)
 *
 * */



#include "util.h"



#define MAXATOM 5
enum { HARD_SPHERE, LJ_FULL, LJ_REPULSIVE };
enum { IE_PY, IE_HNC };

typedef struct {
  int ns;
  double sigma[MAXATOM];
  double eps6_12[MAXATOM][2];
  double rho[MAXATOM];
  double Lpm[MAXATOM*(MAXATOM-1)/2];
  double beta;
  int ljtype;

  int ietype;
  double rmax; /* radius cutoff */
  int npt; /* number of sampling points along r */
  int Mpt; /* number of points for the Newton-Raphson method (LMV) */
  double lmvdamp; /* damping factor for the LMV solver */

  int nlambda; /* number of intermediate stages */

  double charge[MAXATOM];
  double ampch;
  double radch; /* screening distance */

  double C6_12[MAXATOM*(MAXATOM+1)/2][2]; /* alternative to sigma/epsilon */
} model_t;

#include "models.h" /* built-in models from literature */

int model_id = 15;

enum { SOLVER_PICARD, SOLVER_LMV };

int solver = SOLVER_LMV;
double damp = 0.01;
double errinf = 1e9;
int itmax = 100000;
double tol = 1e-7;
int verbose = 1;



/* Lennard-Jones potential in terms of sigma and epsilon */
static double ljpot(double r, double sig, double eps6, double eps12,
    double lam)
{
  double ir, ir6, u1, u2, eps;

  ir = sig/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  if ( eps6 > 0 ) { /* LJ with an attractive tail */
    if ( ir6 < eps6/(2*eps12) ) { /* r > rm */
      u1 = 0; /* repulsive */
      u2 = 4*ir6*(ir6*eps12 - eps6); /* attractive */
    } else { /* r < rm */
      eps = eps6 * eps6 / eps12;
      u1 = 4*ir6*(ir6*eps12 - eps6) + eps; /* repulsive */
      u2 = -eps; /* attractive */
    }
  } else { /* purely repulsive */
    u1 = 4*ir6*(ir6*eps12 - eps6); /* repulsive */
    u2 = 0; /* no attractive tail */
  }
  return u1 + u2 * lam;
}



/* Lennard-Jones potential in terms of C6 and C12 */
static double ljpot6_12(double r, double c6, double c12, double lam)
{
  double sig, eps, ir, ir6;

  if ( c6 < 0 ) { /* attractive */
    sig = pow(-c12/c6, 1./6);
    eps = c6*c6/4/c12;
    return ljpot(r, sig, eps, eps, lam);
  }
  ir = 1/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  return ir6 * (c12 * ir6 + c6);
}




/* repulsive Lennard-Jones potential */
static double ljrpot(double r, double sig, double eps6, double eps12)
{
  double ir, ir6, eps = 0;

  ir = sig/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  if ( ir6 < eps6/(2*eps12) ) return 0;
  if ( eps6 > 0 && eps12 > 0 )
    eps = eps6 * eps6 / eps12;
  return 4*ir6*(ir6*eps12 - eps6) + eps;
}



/* initialize f(r) */
static void initfr(model_t *m, double **fr, double **vrlr, double lam)
{
  int i, j, ij, ji, ipr, l, ns = m->ns, use_c6_12;
  double beta = m->beta, z, u, uelec, ulr;
  double sig, eps6, eps12, c6, c12;

  for ( ipr = 0, i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++, ipr++ ) {
      ij = i*ns + j;
      ji = j*ns + i;
      sig = .5 * (m->sigma[i] + m->sigma[j]);
      eps6 = sqrt(m->eps6_12[i][0] * m->eps6_12[j][0]);
      eps12 = sqrt(m->eps6_12[i][1] * m->eps6_12[j][1]);

      c6 = m->C6_12[ipr][0];
      c12 = m->C6_12[ipr][1];
      use_c6_12 = (fabs(eps6) < DBL_MIN && fabs(eps12) < DBL_MIN);

      for ( l = 0; l < m->npt; l++ ) {
        if ( m->ljtype  == HARD_SPHERE) {
          z = (fft_ri[l] < sig) ? -1 : 0;
        } else { /* Lennard-Jones */
          if ( use_c6_12 ) {
            u = ljpot6_12(fft_ri[l], c6, c12, lam);
          } else if ( m->ljtype == LJ_REPULSIVE ) {
            u = ljrpot(fft_ri[l], sig, eps6, eps12);
          } else {
            u = ljpot(fft_ri[l], sig, eps6, eps12, lam);
          }
          uelec = lam * m->ampch * m->charge[i] * m->charge[j] / fft_ri[l];
          /* set the screen length as sig */
          ulr = uelec * erf( fft_ri[l]/sqrt(2)/m->radch );
          vrlr[ij][l] = beta * ulr;
          if ( j > i ) vrlr[ji][l] = vrlr[ij][l];

          u += uelec - ulr; /* add the short-ranged component */
          z = exp(-beta*u) - 1;
        }
        fr[ij][l] = z;
        if ( j > i ) fr[ji][l] = fr[ij][l];
      }
    }
  }
}



/* initialize the w matrix for intra-molecular covalence bonds */
static void initwk(model_t *m, double **wk)
{
  int i, j, u, ipr, ns = m->ns;

  for ( u = 0; u < m->npt; u++ ) {
    double k = fft_ki[u];
    for ( ipr = 0, i = 0; i < ns; i++ ) {
      wk[i*ns + i][u] = 1; /* for j == i */
      for ( j = i + 1; j < ns; j++, ipr++ ) { /* for j > i */
        double l = m->Lpm[ipr];
        wk[j*ns + i][u] = wk[i*ns + j][u] = (l > 0) ? sin(k*l)/(k*l) : 0;
      }
    }
  }
}



/* Ornstein-Zernike relation: c(k) --> t(k) */
static void oz(model_t *m, double **ck, double **vklr,
    double **tk, double **wk, double **invwc1w)
{
  int i, j, ij, l, u, ns = m->ns;
  double *wc, *wc1, *invwc1, *iwc1w;
  double hk;

  newarr(wc,     ns * ns);
  newarr(wc1,    ns * ns);
  newarr(invwc1, ns * ns);
  newarr(iwc1w,  ns * ns);

  for ( l = 0; l < m->npt; l++ ) {
    /* compute w c
     * note that w c is not symmetric w.r.t. i and j */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i*ns + j;
        wc[ij] = 0;
        for ( u = 0; u < ns; u++ )
          wc[ij] += wk[i*ns + u][l] * (ck[u*ns + j][l] - vklr[u*ns + j][l]);
        wc1[ij] = (i == j) - m->rho[i] * wc[ij];
      }

    /* compute the inverse matrix */
    invmat(ns, wc1, invwc1);

    /* compute (1 - rho w * c)^(-1) w */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i * ns + j;
        iwc1w[ij] = 0;
        for ( u = 0; u < ns; u++ )
          iwc1w[ij] += invwc1[i*ns + u] * wk[u*ns + j][l];
        if ( invwc1w != NULL )
          invwc1w[ij][l] = iwc1w[ij];
      }

    /* compute w c (1 - rho w * c)^(-1) w */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i * ns + j;
        hk = 0;
        for ( u = 0; u < ns; u++ )
          hk += wc[i*ns + u] * iwc1w[u*ns + j];
        tk[ij][l] = hk - ck[ij][l];
      }
  }

  delarr(wc);
  delarr(wc1);
  delarr(invwc1);
  delarr(iwc1w);
}



/* return the cavity distribution function */
static double getyr(double tr, double *dy, int ietype)
{
  double xp;
  if ( ietype == IE_HNC ) {
    xp = exp(tr);
    if ( dy ) *dy = xp;
    return xp;
  } else {
    if ( dy ) *dy = 1;
    return 1 + tr;
  }
}



/* direct Picard iteration
 * do not use this unless for simple models */
static double iter_picard(model_t *model,
    double **fr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, int *niter)
{
  int it, ns = model->ns, i, j, ij, l;
  double y, dcr, err, errp = errinf;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(cr, ck, ns, NULL);
    oz(model, ck, vklr, tk, wk, NULL);
    sphr_k2r(tk, tr, ns, NULL);

    /* solve the closure */
    err = 0;
    for ( l = 0; l < model->npt; l++ ) {
      for ( i = 0; i < ns; i++ ) {
        for ( j = 0; j < ns; j++ ) {
          ij = i * ns + j;
          y = getyr(tr[ij][l], NULL, model->ietype);
          dcr = (fr[ij][l] + 1) * y - 1 - tr[ij][l] - cr[ij][l];
          if ( fabs(dcr) > err ) err = fabs(dcr);
          cr[ij][l] += dcr * damp;
        }
      }
    }

    if ( err < tol ) break;
    if ( verbose )
      fprintf(stderr, "it %d err %g -> %g\n", it, errp, err); //getchar();
    if ( err > errp ) break;
    errp = err;
  }
  *niter = it;
  return err;
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
  int i, j, k, l, m, ij, it, M, npr, ipr, Mp;
  int ns = model->ns, npt = model->npt;
  double **Cjk = NULL, *mat = NULL, *a = NULL, *b = NULL, *costab = NULL;
  double y, dy, err1 = 0, err2 = 0, err, errp = errinf, dmp;
  int *prmask, nsv;

  /* initialize the prmask for solvent-solvent iteraction */
  xnew(prmask, ns * ns);
  for ( i = 1; i < ns; i++ )
    if ( model->Lpm[i-1] < DBL_MIN ) break;
  nsv = i;
  for ( i = 0; i < ns; i++ )
    for ( j = 0; j < ns; j++ )
      prmask[i*ns + j] = (i < nsv && j < nsv);

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
  Mp = M*npr;

  /* initialize the cosine table */
  if ( M > 0 ) {
    newarr2d(Cjk, ns*ns, M*M);
    newarr(mat, Mp*Mp);
    newarr(a, Mp);
    newarr(b, Mp);
    newarr(costab, 4*M*npt);
    for ( i = 0; i < npt; i++ )
      for ( j = 0; j < 4*M; j++ )
        costab[i*4*M + j] = cos(PI*(i+.5)*(j-2*M)/npt);
  }

  for ( it = 0; it < itmax; it++ ) {
    /* compute c(r) and c(k) from the closure */
    for ( i = 0; i < ns; i++ ) {
      for ( j = 0; j < ns; j++ ) {
        for ( l = 0; l < npt; l++ ) {
          ij = i * ns + j;
          y = getyr(tr[ij][l], &dy, model->ietype);
          cr[ij][l] = (fr[ij][l] + 1) * y - tr[ij][l] - 1;
          der[ij][l] = (fr[ij][l] + 1) * dy - 1;
        }
      }
    }
    sphr_r2k(cr, ck, ns, NULL);

    /* compute Cjk */
    for ( i = 0; i < ns; i++ ) {
      for ( j = i; j < ns; j++ ) {
        ij = i * ns + j;

        for ( m = 0; m < M; m++ ) {
          for ( k = 0; k < M; k++ ) {
            y = 0;
            for ( l = 0; l < npt; l++ )
              y += der[ij][l] * (costab[l*4*M+k-m+2*M]
                               - costab[l*4*M+k+m+2*M]);
            y /= npt;
            Cjk[ij][m*M+k] = y;
            if ( j > i ) Cjk[j*ns+i][m*M+k] = y;
          }
        }

      }
    }

    /* compute the new t(k) */
    oz(model, ck, vklr, ntk, wk, invwc1w);

    /* compute the Jacobian for the Newton-Raphson method */
    for ( m = 0; m < M; m++ ) {
      for ( ipr = 0, i = 0; i < ns; i++ ) {
        for ( j = i; j < ns; j++, ipr++ ) {
          int id1, id2, i2, j2, ipr2;

          id1 = m*npr + ipr;
          if ( prmask[i*ns + j] )
            b[id1] = fft_ki[m] * (ntk[i*ns+j][m] - tk[i*ns+j][m]);
          else
            b[id1] = 0;

          for ( k = 0; k < M; k++ ) {
            for ( ipr2 = 0, i2 = 0; i2 < ns; i2++ ) {
              for ( j2 = i2; j2 < ns; j2++, ipr2++ ) {
                id2 = k*npr + ipr2;
                if (id1 >= Mp || id2 >= Mp) {
                  fprintf(stderr, "id1 %d, id2 %d, k %d, ipr2 %d/%d\n", id1, id2, k, ipr2, npr);
                  exit(1);
                }
                y = (ipr == ipr2 ? (m == k) + Cjk[i*ns+j][m*M+k] : 0)
                  - invwc1w[i*ns+i2][m] * Cjk[i2*ns+j2][m*M+k]
                  * invwc1w[j2*ns+j][k];
                mat[id1*Mp + id2] = y;
              }
            }
          }
        }
      }
    }

    if ( linsolve(Mp, mat, a, b) != 0 )
      break;

    /* use the Newton-Raphson method to solve for t(k) of small k */
    err1 = 0;
    for ( m = 0; m < M; m++ ) {
      for ( ipr = 0, i = 0; i < ns; i++ ) {
        for ( j = i; j < ns; j++, ipr++ ) {
          y = a[m*npr+ipr] / fft_ki[m];
          if ( fabs(y) > err1 ) err1 = fabs(y);
          tk[i*ns+j][m] += dmp * y;
          if ( j > i ) tk[j*ns+i][m] = tk[i*ns + j][m];
        }
      }
    }

    /* use the OZ relation to solve for t(k) of large k */
    err2 = 0;
    for ( l = M; l < npt; l++ ) {
      for ( i = 0; i < ns; i++ ) {
        for ( j = 0; j < ns; j++ ) {
          ij = i * ns + j;
          y = prmask[ij] ? ntk[ij][l] - tk[ij][l] : 0;
          if ( fabs(y) > err2 ) err2 = fabs(y);
          tk[ij][l] += dmp * y;
        }
      }
    }

    sphr_k2r(tk, tr, ns, NULL);

    if ( verbose )
      fprintf(stderr, "it %d: M %d, errp %g, err1 %g, err2 %g, damp %g\n",
          it, M, errp, err1, err2, dmp);
    err = err1 > err2 ? err1 : err2;
    if ( err < tol ) {
      if ( prmask[ns*ns-1] ) break;
      if ( prmask[nsv] == 0 ) { /* turn on solute-solvent interaction */
        for ( i = 0; i < ns; i++ )
          for ( j = 0; j < ns; j++ )
            prmask[i*ns + j] = (i < nsv && j >= nsv || j < nsv && i >= nsv);
        fprintf(stderr, "turning on solute-solvent interaction\n"); //getchar();
      } else { /* turn on solute-solute interaction */
        for ( i = 0; i < ns; i++ )
          for ( j = 0; j < ns; j++ )
            prmask[i*ns + j] = (i >= nsv && j >= nsv);
        oz(model, ck, vklr, tk, wk, invwc1w);
        sphr_k2r(tk, tr, ns, NULL);
        /* NOTE currently the we cannot continue from here
         * thus, we stop at the case of infinite dilution */
        //fprintf(stderr, "turning on solute-solute interaction\n"); getchar();
        break;
      }
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
    delarr(a);
    delarr(b);
  }
  free(prmask);
  return err;
}



/* save correlation functions to file `fn' */
static int output(model_t *model,
    double **cr, double **vrlr, double **tr, double **fr,
    double **ck, double **vklr, double **tk, double **wk,
    const char *fn)
{
  int i, j, ij, l, ns = model->ns;
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i * ns + j;
      for ( l = 0; l < model->npt; l++ ) {
        double vr = vrlr[ij][l], vk = vklr[ij][l];
        fprintf(fp, "%g %g %g %g %g %d %d %g %g %g %g %g\n",
            fft_ri[l], cr[ij][l], tr[ij][l], fr[ij][l],
            vr, i, j,
            fft_ki[l], ck[ij][l], tk[ij][l], wk[ij][l],
            vk);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return 0;
}



/* compute the equivalent diameter of the molecule
 * this value is only used for comparison
 * and it does not affect the solver */
static double getdiameter(model_t *m)
{
  int i, j, ipr, ns = m->ns;
  double vol = 0, si, sj, l;

  /* determine the number of sites of the solvent */
  for ( i = 1; i < ns; i++ )
    if ( m->Lpm[i-1] <= 0 ) {
      ns = i;
      break;
    }

  for ( i = 0; i < ns; i++ )
    vol += pow( m->sigma[i], 3 );

  /* deduct the overlap */
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns; j++, ipr++ ) {
      si = m->sigma[i];
      sj = m->sigma[j];
      l = m->Lpm[ipr];
      vol -= (si*si*si + sj*sj*sj)/2 - (si*si + sj*sj)*l*3./4
           - pow(si*si - sj*sj, 2)*3./32/l + l*l*l/2;
    }

  return pow(vol, 1./3);
}



static void dorism(void)
{
  int it, ns, npt, ilam, nlam;
  double err, dia, lam;
  double **fr, **wk, **cr, **cp, **ck, **tr, **tk;
  double **der, **ntk, **vrlr, **vklr;
  model_t *model = models + model_id;

  /* equivalent diameter of the solvent molecule */
  dia = getdiameter(model);

  ns = model->ns;
  npt = model->npt;
  newarr2d(fr,    ns * ns, npt);
  newarr2d(wk,    ns * ns, npt);
  newarr2d(cr,    ns * ns, npt);
  newarr2d(ck,    ns * ns, npt);
  newarr2d(cp,    ns * ns, npt);
  newarr2d(tr,    ns * ns, npt);
  newarr2d(tk,    ns * ns, npt);
  newarr2d(der,   ns * ns, npt);
  newarr2d(ntk,   ns * ns, npt);
  newarr2d(vrlr,  ns * ns, npt);
  newarr2d(vklr,  ns * ns, npt);

  initfftw(model->rmax, npt);
  initwk(model, wk);

  /* lambda is used to gradually switch on long-range interaction */
  nlam = model->nlambda;
  if ( nlam < 1 ) nlam = 1;

  for ( ilam = 1; ilam <= nlam; ilam++ ) {
    lam = 1.*ilam/nlam;

    initfr(model, fr, vrlr, lam);
    sphr_r2k(vrlr, vklr, ns, NULL);

    /* use f(r) as the initial c(r) for the lowest lambda */
    if ( ilam == 1)
      cparr2d(cr, fr, ns * ns, npt);

    if ( solver == SOLVER_LMV ) {
      err = iter_lmv(model, fr, wk, cr, der, ck, vklr, tr, tk, ntk, cp, &it);
    } else {
      err = iter_picard(model, fr, wk, cr, ck, vklr, tr, tk, &it);
    }
    output(model, cr, vrlr, tr, fr, ck, vklr, tk, wk, "out.dat");
    fprintf(stderr, "lambda %g, %d iterations, err %g, d %g, rho*d^3 %g\n",
        lam, it, err, dia, model->rho[0]*dia*dia*dia);
  }

  delarr2d(fr,    ns * ns);
  delarr2d(wk,    ns * ns);
  delarr2d(cr,    ns * ns);
  delarr2d(ck,    ns * ns);
  delarr2d(cp,    ns * ns);
  delarr2d(tr,    ns * ns);
  delarr2d(tk,    ns * ns);
  delarr2d(der,   ns * ns);
  delarr2d(ntk,   ns * ns);
  delarr2d(vrlr,  ns * ns);
  delarr2d(vklr,  ns * ns);
  donefftw();
}



int main(void)
{
  dorism();
  return 0;
}
