/* basic RISM, assume a single solvent molecule
 *
 *  gcc rism0.c -lfftw3 -lm
 *
 * the output can be found in "out.dat"
 * to view g(r) in gnuplot, type
 *
 *  plot "out.dat" u 1:(1+$2+$3)
 *
 * */



#include "debug.h"
#include "util.h"
#include "model.h"
#include "calctd.h"



int model_id = 16;
int skip_uu = 0;
int verbose = 0;
const char *fncrtr = "out.dat";
int sepout = 0;
const char *fncrdnum = "crdnum.dat";
int printk = 0;



/* print help message and die */
static void help(const char *prog)
{
  fprintf(stderr,
      "Reference site interaction model (RISM)\n"
      "Usage:\n"
      "  %s [Options] [input.cfg|model_id]\n"
      "\n"
      "Options:\n"
      "  -o:    followed by the output file, default: %s\n"
      "  -k:    print k-space correlation functions, default %d\n"
      "  -!:    skip the solute-solute stage calculation, default: %d\n"
      "  -$:    separately output file for each lambda, default %d\n"
      "  -#:    followed by the coordination number file, default: %s\n"
      "  -v:    be verbose\n"
      "  -vv:   be more verbose\n"
      "  -h:    display this message\n",
      prog, fncrtr, printk, skip_uu, sepout, fncrdnum);
  exit(1);
}



/* handle command line arguments */
static model_t *doargs(int argc, char **argv)
{
  model_t *m;
  int i, j, ch;
  const char *p, *q, *fncfg = NULL;

  for ( i = 1; i < argc; i++ ) {
    /* it's an argument */
    if ( argv[i][0] != '-' ) {
      if ( striscnum(argv[i]) ) { /* use the stock model */
        model_id = atoi(argv[i]);
      } else { /* load a configuration file */
        model_id = 0;
        fncfg = argv[i];
      }
      continue;
    }

    /* it is an option
     * loop over characters in the options
     * in this way, `-vo' is understood as `-v -o' */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( ch == 'o' || ch == '#' ) {
        /* handle options that require an argument */
        q = p = argv[i] + j + 1;
        if ( *p != '\0' ) {
          /* the argument follows immediately after the option
           * e.g., -oa.dat */
          q = p;
        } else if ( ++i < argc ) {
          /* the option and argument are separated by a space
           * then the argument belongs to the next argv[] element,
           * hence ++i
           * e.g., -o a.dat */
          q = argv[i];
        } else {
          fprintf(stderr, "-%c requires an argument!\n", ch);
          help(argv[0]);
        }
        if ( ch == 'o' ) {
          fncrtr = q;
        } else if ( ch == '#' ) {
          fncrdnum = q;
        }
        break; /* skip the rest of the characters in the option */
      } else if ( ch == '!' ) {
        skip_uu = 1;
      } else if ( ch == '$' ) {
        sepout = 1;
      } else if ( ch == 'v' ) {
        verbose++;
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        help(argv[0]);
      }
    }
  }

  m = models + model_id;
  if ( fncfg != NULL )
    if ( model_load(m, fncfg, verbose) != 0 ) {
      fprintf(stderr, "failed to load %s\n", fncfg);
      help(argv[0]);
    }

  return m;
}



/* Lennard-Jones potential in terms of sigma and epsilon */
static double ljpot(double r, double sig, double eps6, double eps12,
    double lam, double *nrdu)
{
  double ir, ir6, u1, u2, eps;

  ir = sig/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  if ( eps6 > 0 ) { /* LJ with an attractive tail */
    if ( ir6 < eps6/(2*eps12) ) { /* r > rm */
      u1 = 0; /* repulsive */
      u2 = 4*ir6*(ir6*eps12 - eps6); /* attractive */
      *nrdu = ir6 * (48 * ir6 * eps12 - 24 * eps6) * lam;
    } else { /* r < rm */
      eps = eps6 * eps6 / eps12;
      u1 = 4*ir6*(ir6*eps12 - eps6) + eps; /* repulsive */
      u2 = -eps; /* attractive */
      *nrdu = ir6 * (48 * ir6 * eps12 - 24 * eps6);
    }
  } else { /* purely repulsive */
    u1 = 4*ir6*(ir6*eps12 - eps6); /* repulsive */
    u2 = 0; /* no attractive tail */
    *nrdu = ir6 * (48 * ir6 * eps12 - 24 * eps6);
  }
  return u1 + u2 * lam;
}



/* Lennard-Jones potential in terms of C6 and C12 */
static double ljpot6_12(double r, double c6, double c12,
    double lam, double *nrdu)
{
  double sig, eps, ir, ir6;

  if ( c6 < 0 ) { /* attractive */
    sig = pow(-c12/c6, 1./6);
    eps = c6*c6/4/c12;
    return ljpot(r, sig, eps, eps, lam, nrdu);
  }
  ir = 1/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  *nrdu = ir6 * (12 * c12 * ir6 + 6 * c6);
  return ir6 * (c12 * ir6 + c6);
}



/* repulsive Lennard-Jones potential */
static double ljrpot(double r, double sig, double eps6, double eps12,
    double *nrdu)
{
  double ir, ir6, eps = 0;

  ir = sig/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  if ( ir6 < eps6/(2*eps12) ) return 0;
  if ( eps6 > 0 && eps12 > 0 )
    eps = eps6 * eps6 / eps12;
  *nrdu = ir6*(48*ir6*eps12 - 24*eps6);
  return 4*ir6*(ir6*eps12 - eps6) + eps;
}



/* exponential potential: B * exp(-r/rho) */
static double exppot(double r, double B, double rho, double *nrdu)
{
  r /= rho;
  B *= exp(-r);
  *nrdu = B*r;
  return B;
}



/* initialize f(r) */
static void initfr(model_t *m, double **ur, double **nrdur,
    double **fr, double **vrqq, double **vrlr, double **vrsr,
    double lam)
{
  int i, j, ij, ji, ipr, l, ns = m->ns, use_pairpot;
  double beta = m->beta, r, z, u, uelec, nrdu, nrdu2, ulr;
  double sig, eps6, eps12, c6, c12, Bij, rhoij;

  for ( ipr = 0, i = 0; i < ns; i++ ) { /* the first site */
    for ( j = i; j < ns; j++, ipr++ ) { /* the second site */
      ij = i*ns + j;
      ji = j*ns + i;
      sig = .5 * (m->sigma[i] + m->sigma[j]);
      eps6 = sqrt(m->eps6_12[i][0] * m->eps6_12[j][0]);
      eps12 = sqrt(m->eps6_12[i][1] * m->eps6_12[j][1]);

      c6 = m->pairpot[ipr].C6;
      c12 = m->pairpot[ipr].C12;
      Bij = m->pairpot[ipr].B;
      rhoij = m->pairpot[ipr].rho;
      use_pairpot = (fabs(eps6) < DBL_MIN && fabs(eps12) < DBL_MIN);
      if ( use_pairpot ) {
        sig = m->pairpot[ipr].sigma;
        eps6 = m->pairpot[ipr].eps6;
        eps12 = m->pairpot[ipr].eps12;
      }

      for ( l = 0; l < m->npt; l++ ) { /* the radius */
        r = fft_ri[l];
        if ( m->ljtype  == HARD_SPHERE) {
          if (r < sig) {
            vrsr[ij][l] = 2*INFTY;
            z = -1;
          } else {
            vrsr[ij][l] = 0;
            z = 0;
          }
          nrdur[ij][l] = ur[ij][l] = 0;
          vrlr[ij][l] = vrqq[ij][l] = 0;
        } else { /* Lennard-Jones */
          if ( use_pairpot ) {
            if ( fabs(c6) > DBL_MIN || fabs(c12) > DBL_MIN ) {
              u = ljpot6_12(r, c6, c12, lam, &nrdu);
            } else {
              u = ljpot(r, sig, eps6, eps12, lam, &nrdu);
            }
            u += exppot(r, Bij, rhoij, &nrdu2);
            nrdu += nrdu2;
          } else if ( m->ljtype == LJ_REPULSIVE ) {
            u = ljrpot(r, sig, eps6, eps12, &nrdu);
          } else {
            u = ljpot(r, sig, eps6, eps12, lam, &nrdu);
          }
          uelec = lam * m->ampch * m->charge[i] * m->charge[j] / r;
          ur[ij][l] = u + uelec;
          nrdur[ij][l] = nrdu + uelec;
          ulr = uelec * erf( r/sqrt(2)/m->rscreen );
          vrqq[ij][l] = beta * uelec;
          vrlr[ij][l] = beta * ulr;
          vrsr[ij][l] = beta * (u + uelec - ulr);

          z = exp(-vrsr[ij][l]) - 1;
        }
        fr[ij][l] = z;
        if ( j > i ) {
          ur[ji][l] = ur[ij][l];
          nrdur[ji][l] = nrdur[ij][l];
          vrqq[ji][l] = vrqq[ij][l];
          vrlr[ji][l] = vrlr[ij][l];
          vrsr[ji][l] = vrsr[ij][l];
          fr[ji][l] = fr[ij][l];
        }
      } /* loop over l, the radius */
    } /* loop over j, the second site */
  } /* loop over i, the first site */
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
        double l = m->dis[ipr];
        wk[j*ns + i][u] = wk[i*ns + j][u] = (l > 0) ? sin(k*l)/(k*l) : 0;
      }
    }
  }
}



/* Ornstein-Zernike relation: c(k) --> t(k) */
static void oz(model_t *m, double **ck, double **vklr,
    double **tk, double **wk, double **invwc1w)
{
  int i, j, ij, l, ns = m->ns;
  double w[NS2MAX], dw[NS2MAX], c[NS2MAX], wc[NS2MAX], invwc1[NS2MAX];
  double tm1[NS2MAX], tm2[NS2MAX], tm3[NS2MAX];

  for ( l = 0; l < m->npt; l++ ) {
    for ( ij = 0; ij < ns * ns; ij++ ) {
      w[ij] = wk[ij][l];
      c[ij] = ck[ij][l] - vklr[ij][l];
    }

    /* note that w c is not symmetric w.r.t. i and j */
    matmul(wc, w, c, ns); /* wc = w c */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i*ns + j;
        tm1[ij] = m->rho[i] * wc[ij]; /* tm1 = rho w c */
        tm2[ij] = (i == j) - tm1[ij]; /* tm2 = 1 - rho w c */
        dw[ij] = wk[ij][l] - (i == j);
      }

    invmat(tm2, invwc1, ns); /* invwc1 = (1 - wc)^(-1) */

    matmul(tm3, invwc1, w, ns); /* tm3 = (1 - rho w c)^(-1) w */
    if ( invwc1w != NULL )
      for ( ij = 0; ij < ns*ns; ij++ )
        invwc1w[ij][l] = tm3[ij];
    matmul(tm2, tm1, tm3, ns); /* tm2 = rho w c (1 - rho w c)^(-1) w */
    matmul(tm1, wc, tm2, ns); /* tm1 = w c rho w c (1 - rho w c)^(-1) w */

    /* w c w - c = w c (1 + dw) - c = dw c + w c dw */
    matmul(tm2, dw, c, ns);
    matmul(tm3, wc, dw, ns);

    /* compute w c (1 - rho w c)^(-1) w - c
     * = [w c rho w c (1 - rho w c)^(-1) w - w c w] + (w c w - c) */
    for ( ij = 0; ij < ns * ns; ij++ )
      tk[ij][l] = tm1[ij] + tm2[ij] + tm3[ij] - vklr[ij][l];
  }
}



/* return the updated cr */
static double getcr(double tr, double vrsr,
                    double *dcr, int ietype)
{
  double xp, del, fr;

  del = -vrsr + tr;
  if ( ietype == IE_HNC ) {
    xp = exp(del);
    if ( dcr != NULL ) *dcr = xp - 1;
    return xp - 1 - tr;
  } else if ( ietype == IE_PY ) {
    fr = exp(-vrsr) - 1;
    if ( dcr != NULL ) *dcr = fr;
    return fr * (1 + tr);
  } else if ( ietype == IE_KH ) {
    if ( del <= 0 ) { /* HNC */
      xp = exp(del);
      if ( dcr != NULL ) *dcr = xp - 1;
      return xp - 1 - tr;
    } else {
      if ( dcr != NULL ) *dcr = 0;
      return -vrsr;
    }
  } else {
    fprintf(stderr, "unknown closure %d\n", ietype);
    exit(1);
  }
  return 0;
}



struct { int l; double err; } errspot[NS2MAX];

/* apply the closure
 * compute residue vector if needed */
static double closure(model_t *model,
    double *res, double **der, double **vrsr,
    double **cr, double **tr,
    int *prmask, int update, double damp)
{
  int ns = model->ns, npt = model->npt, i, j, ij, ji, id, l;
  double y, err, max, errm = 0;

  for ( id = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;
      ji = j*ns + i;
      if ( prmask && !prmask[ij] ) continue;
      err = max = 0;
      for ( l = 0; l < npt; l++, id++ ) {
        y = getcr(tr[ij][l], vrsr[ij][l], der ? &der[ij][l] : NULL,
                  model->ietype) - cr[ij][l];
        if ( res != NULL ) res[id] = y;
        if ( update ) {
          cr[ij][l] += damp * y;
          if ( j > i ) cr[ji][l] = cr[ij][l];
        }
        if ( fabs(y) > err ) {
          err = fabs(y);
          errspot[ij].l = l;
          errspot[ij].err = err;
        }
        if ( fabs(cr[ij][l]) > max ) max = fabs(cr[ij][l]);
      }
      /* the c(r) between two ions can be extremely large
       * so we use the relative error to be compared with the tolerance */
      if ( (err /= (max + 1e-6)) > errm ) errm = err;
    }
  return errm;
}




/* a step of direct iteration (Picard)
 * compute residue vector if needed */
static double step_picard(model_t *model,
    double *res, double **der,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk,
    int *prmask, int update, double damp)
{
  sphr_r2k(cr, ck, model->ns, NULL);
  oz(model, ck, vklr, tk, wk, NULL);
  sphr_k2r(tk, tr, model->ns, NULL);
  return closure(model, res, der, vrsr, cr, tr,
                 prmask, update, damp);
}



const double errinf = 1e20;
#include "uv.h" /* manager for solvent/solute interaction */
#include "lmv.h" /* LMV solver */
#include "mdiis.h" /* MDIIS solver */



/* direct Picard iteration
 * do not use this unless for simple models
 * does not handle solvent-solute interactions */
static double iter_picard(model_t *model,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk,
    int *niter)
{
  int it;
  double err = 0, errp = errinf;

  for ( it = 0; it < model->itmax; it++ ) {
    err = step_picard(model, NULL, NULL, vrsr, wk, cr, ck, vklr,
        tr, tk, NULL, 1, model->picard.damp);
    if ( err < model->tol ) break;
    if ( verbose )
      fprintf(stderr, "it %d err %g -> %g\n", it, errp, err); //getchar();
    if ( err > errp ) break;
    errp = err;
  }
  *niter = it;
  return err;
}



/* save correlation functions to file `fn' */
static int output(model_t *m,
    double **cr, double **vrqq, double **vrlr,
    double **ur, double **tr, double **fr,
    double **ck, double **vklr, double **tk, double **wk,
    const char *fn, int ilam)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt;
  FILE *fp;
  char fnl[80];
  double eps_rism;

  if ( sepout ) {
    sprintf(fnl, "%s%d", fn, ilam);
  } else {
    strcpy(fnl, fn);
  }
  if ((fp = fopen(fnl, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fnl);
    return -1;
  }
  eps_rism = calcdielec(m);
  /* print some basic information on the first line */
  fprintf(fp, "# %g %g %g %d %g\n",
      1/(m->kBT*m->beta), m->kBU, m->ampch, m->ietype, eps_rism);
  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;
      for ( l = 0; l < npt; l++ ) {
        double vrl = vrlr[ij][l], vkl = vklr[ij][l];
        double vrt = m->beta * ur[ij][l], vrq = vrqq[ij][l];
        /* note that cr[ij][l] and tr[ij][l] exclude
         * the long-range component vrl */
        double crt = cr[ij][l] - vrl, trt = tr[ij][l] + vrl;
        double ckt = ck[ij][l] - vkl, tkt = tk[ij][l] + vkl;
        /* pmfs = beta dW: the short-range correction to
         * the continuum primitive model, in which
         *
         *  W_c = u_LJ + u_qq/eps  (with u = u_LJ + u_qq)
         *
         * `pmfs' is computed from Eq. (50) of the following paper:
         * [ ``The interionic potential of mean force in a molecular polar
         *     solvent from an extended RISM equation''
         *    Hirate, Rossky, and Pettitt,
         *    J. Chem. Phys. 78(6) 4133-4144 (1983) ]
         *
         *    beta dW
         *  = beta W_rism_corrected - beta W_s
         *  = (beta W_s - phi_qq/eps) - beta W_c
         *  = beta W_s + phi*                 (here, phi* = -beta u_LJ)
         *  = beta W + phi/eps_rism + phi*    (here, phi  = -beta u_qq)
         *  = beta u - t(r) - beta u_qq/eps_rism - beta u_LJ
         *  = beta u_q - t(r) - beta u_qq/eps_rism
         * */
        double pmfs = vrq - trt - vrq/eps_rism;
        fprintf(fp, "%g %g %g %g %g %g %d %d %g %g ",
            fft_ri[l], crt, trt, fr[ij][l],
            vrl, vrt, i, j, vrq, pmfs);
        if ( printk )
          fprintf(fp, "%g %g %g %g %g",
              fft_ki[l], ckt, tkt, wk[ij][l], vkl);
        fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  fprintf(stderr, "saved result to %s\n", fnl);
  return 0;
}



static int dorism(model_t *model)
{
  int it, ns, npt, ilam, nlam;
  double err = 0, dia, lam;
  double **ur, **nrdur, **fr, **wk;
  double **cr, **cp, **ck, **tr, **tk, **ntk;
  double **der, **vrlr, **vrqq, **vrsr, **vklr;
  double *um, *mum;

  /* equivalent diameter of the solvent molecule */
  dia = getdiameters(model);

  ns = model->ns;
  npt = model->npt;
  newarr2d(ur,    ns * ns, npt);
  newarr2d(nrdur, ns * ns, npt);
  newarr2d(vrqq,  ns * ns, npt);
  newarr2d(vrlr,  ns * ns, npt);
  newarr2d(vrsr,  ns * ns, npt);
  newarr2d(vklr,  ns * ns, npt);
  newarr2d(fr,    ns * ns, npt);
  newarr2d(wk,    ns * ns, npt);
  newarr2d(cr,    ns * ns, npt);
  newarr2d(ck,    ns * ns, npt);
  newarr2d(cp,    ns * ns, npt);
  newarr2d(tr,    ns * ns, npt);
  newarr2d(tk,    ns * ns, npt);
  newarr2d(ntk,   ns * ns, npt);
  newarr2d(der,   ns * ns, npt);
  xnew(um, ns);
  xnew(mum, ns);

  initfftw(model->rmax, npt);
  initwk(model, wk);

  /* lambda is used to gradually switch on long-range interaction */
  nlam = model->nlambdas;
  if ( nlam < 1 ) nlam = 1;

  for ( ilam = 1; ilam <= nlam; ilam++ ) {
    lam = 1.*ilam/nlam;

    initfr(model, ur, nrdur, fr, vrqq, vrlr, vrsr, lam);
    sphr_r2k(vrlr, vklr, ns, NULL);

    /* use f(r) as the initial c(r) for the lowest lambda */
    if ( ilam == 1 )
      cparr2d(cr, fr, ns * ns, npt);

    if ( model->solver == SOLVER_LMV ) {
      err = iter_lmv(model, vrsr, wk, cr, der, ck, vklr,
          tr, tk, ntk, cp, skip_uu, &it);
    } else if ( model->solver == SOLVER_MDIIS ) {
      err = iter_mdiis(model, vrsr, wk, cr, ck, vklr,
          tr, tk, skip_uu, &it);
    } else {
      err = iter_picard(model, vrsr, wk, cr, ck, vklr,
          tr, tk, &it);
    }

    output(model, cr, vrqq, vrlr, ur, tr, fr, ck, vklr, tk, wk, fncrtr, ilam);
    fprintf(stderr, "lambda %g, %d iterations, err %g, d %g, rho*d^3 %g\n",
        lam, it, err, dia, model->rho[0]*dia*dia*dia);
  }

  calcU(model, ur, cr, tr, fr, um);
  calcchempot(model, cr, tr, vrsr, vrlr, mum);
  calckirk(model, cr, tr, NULL);
  calccrdnum(model, cr, tr, fr, fncrdnum);
  calcdielec(model);

  delarr2d(ur,    ns * ns);
  delarr2d(nrdur, ns * ns);
  delarr2d(vrqq,  ns * ns);
  delarr2d(vrlr,  ns * ns);
  delarr2d(vrsr,  ns * ns);
  delarr2d(vklr,  ns * ns);
  delarr2d(fr,    ns * ns);
  delarr2d(wk,    ns * ns);
  delarr2d(cr,    ns * ns);
  delarr2d(ck,    ns * ns);
  delarr2d(cp,    ns * ns);
  delarr2d(tr,    ns * ns);
  delarr2d(tk,    ns * ns);
  delarr2d(ntk,   ns * ns);
  delarr2d(der,   ns * ns);
  free(um);
  free(mum);
  donefftw();
  return 0;
}



int main(int argc, char **argv)
{
  model_t *m = doargs(argc, argv);
  return dorism(m);
}

