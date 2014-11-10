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



#include "util.h"
#include "model.h"
#include "calctd.h"



int model_id = 16;
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
      "  -$:    separately output file for each lambda, default %d\n"
      "  -#:    followed by the coordination number file, default: %s\n"
      "  -v:    be verbose\n"
      "  -vv:   be more verbose\n"
      "  -h:    display this message\n",
      prog, fncrtr, printk, sepout, fncrdnum);
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



/* initialize f(r) */
static void initfr(model_t *m, double **ur, double **nrdur,
    double **fr, double **vrlr, double **vrsr, double lam)
{
  int i, j, ij, ji, ipr, l, ns = m->ns, use_c6_12;
  double beta = m->beta, z, u, uelec, nrdu, ulr;
  double sig, eps6, eps12, c6, c12;

  for ( ipr = 0, i = 0; i < ns; i++ ) { /* the first site */
    for ( j = i; j < ns; j++, ipr++ ) { /* the second site */
      ij = i*ns + j;
      ji = j*ns + i;
      sig = .5 * (m->sigma[i] + m->sigma[j]);
      eps6 = sqrt(m->eps6_12[i][0] * m->eps6_12[j][0]);
      eps12 = sqrt(m->eps6_12[i][1] * m->eps6_12[j][1]);

      c6 = m->C6_12[ipr][0];
      c12 = m->C6_12[ipr][1];
      use_c6_12 = (fabs(eps6) < DBL_MIN && fabs(eps12) < DBL_MIN);

      for ( l = 0; l < m->npt; l++ ) { /* the radius */
        if ( m->ljtype  == HARD_SPHERE) {
          if (fft_ri[l] < sig) {
            vrsr[ij][l] = 2*INFTY;
            z = -1;
          } else {
            vrsr[ij][l] = 0;
            z = 0;
          }
          nrdur[ij][l] = ur[ij][l] = 0;
          vrlr[ij][l] = 0;
        } else { /* Lennard-Jones */
          if ( use_c6_12 ) {
            u = ljpot6_12(fft_ri[l], c6, c12, lam, &nrdu);
          } else if ( m->ljtype == LJ_REPULSIVE ) {
            u = ljrpot(fft_ri[l], sig, eps6, eps12, &nrdu);
          } else {
            u = ljpot(fft_ri[l], sig, eps6, eps12, lam, &nrdu);
          }
          uelec = lam * m->ampch * m->charge[i] * m->charge[j] / fft_ri[l];
          ur[ij][l] = u + uelec;
          nrdur[ij][l] = nrdu + uelec;
          /* set the screen length as sig */
          ulr = uelec * erf( fft_ri[l]/sqrt(2)/m->rscreen );
          vrlr[ij][l] = beta * ulr;
          vrsr[ij][l] = beta * (u + uelec - ulr);

          z = exp(-vrsr[ij][l]) - 1;
        }
        fr[ij][l] = z;
        if ( j > i ) {
          ur[ji][l] = ur[ij][l];
          nrdur[ji][l] = nrdur[ij][l];
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



/* compute the number solvents
 * here, solvent are defined as molecules with nonzero density
 * zero-density molecules must appear at the end of model */
static int getnsv(model_t *m)
{
  int i;

  for ( i = 0; i < m->ns; i++ )
    if ( m->rho[i] < DBL_MIN ) break;
  return i;
}



/* Ornstein-Zernike relation: c(k) --> t(k) */
static void oz(model_t *m, double **ck, double **vklr,
    double **tk, double **wk, double **invwc1w)
{
  int i, j, ij, l, u, ns = m->ns;
  double *wkl, *ckl, *cvkl, *wc, *wc1, *invwc1, *iwc1w;
  double hk;

  newarr(wkl,    ns * ns);
  newarr(ckl,    ns * ns);
  newarr(cvkl,   ns * ns);
  newarr(wc,     ns * ns);
  newarr(wc1,    ns * ns);
  newarr(invwc1, ns * ns);
  newarr(iwc1w,  ns * ns);

  for ( l = 0; l < m->npt; l++ ) {
    for ( i = 0; i < ns * ns; i++ ) {
      wkl[i] = wk[i][l];
      ckl[i] = ck[i][l];
      cvkl[i] = ckl[i] - vklr[i][l];
    }

    /* compute w c
     * note that w c is not symmetric w.r.t. i and j */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i*ns + j;
        wc[ij] = 0;
        for ( u = 0; u < ns; u++ )
          wc[ij] += wkl[i*ns + u] * cvkl[u*ns + j];
        wc1[ij] = (i == j) - m->rho[i] * wc[ij];
      }

    /* compute the inverse matrix */
    invmat(wc1, invwc1, ns);

    /* compute (1 - rho w * c)^(-1) w */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i*ns + j;
        iwc1w[ij] = 0;
        for ( u = 0; u < ns; u++ )
          iwc1w[ij] += invwc1[i*ns + u] * wkl[u*ns + j];
        if ( invwc1w != NULL )
          invwc1w[ij][l] = iwc1w[ij];
      }

    /* compute w c (1 - rho w * c)^(-1) w */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i*ns + j;
        hk = 0;
        for ( u = 0; u < ns; u++ )
          hk += wc[i*ns + u] * iwc1w[u*ns + j];
        tk[ij][l] = hk - ckl[ij];
      }
  }

  delarr(wkl);
  delarr(ckl);
  delarr(cvkl);
  delarr(wc);
  delarr(wc1);
  delarr(invwc1);
  delarr(iwc1w);
}



/* return the update of cr */
static double getcr(double tr, double vrsr,
                    double *dcr, int ietype)
{
  double xp, del, fr;

  del = -vrsr + tr;
  if ( ietype == IE_HNC ) {
    xp = ( del <= -INFTY ) ? 0 : exp(del);
    if ( dcr != NULL ) *dcr = xp - 1;
    return xp - 1 - tr;
  } else if ( ietype == IE_PY ) {
    fr = ( vrsr >= INFTY ) ? -1 : exp(-vrsr) - 1;
    if ( dcr != NULL ) *dcr = fr;
    return fr * (1 + tr);
  } else if ( ietype == IE_KH ) {
    if ( del <= 0 ) { /* HNC */
      xp = ( del <= -INFTY ) ? 0 : exp(del);
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



/* apply the closure
 * compute residue vector if needed */
static double closure(model_t *model,
    double *res, double **der, double **vrsr,
    double **cr, double **tr, int *prmask,
    int update, double damp)
{
  int ns = model->ns, npt = model->npt, i, j, ij, ji, id, l;
  double y, err = 0, max = 0;

  for ( id = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;
      ji = j*ns + i;
      if ( prmask && !prmask[ij] ) continue;
      for ( l = 0; l < npt; l++, id++ ) {
        y = getcr(tr[ij][l], vrsr[ij][l], der ? &der[ij][l] : NULL,
                  model->ietype) - cr[ij][l];
        if ( res != NULL ) res[id] = y;
        if ( update ) {
          cr[ij][l] += damp * y;
          if ( j > i ) cr[ji][l] = cr[ij][l];
        }
        if ( fabs(y) > err ) err = fabs(y);
        if ( fabs(cr[ij][l]) > max ) max = fabs(cr[ij][l]);
      }
    }
  /* the c(r) between two ions can be extremely large
   * so we use the relative error to be compared with the tolerance */
  return err / (max + 1e-3);
}




/* a step of direct iteration (Picard)
 * compute residue vector if needed */
static double step_picard(model_t *model,
    double *res, double **der,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, int *prmask,
    int update, double damp)
{
  sphr_r2k(cr, ck, model->ns, NULL);
  oz(model, ck, vklr, tk, wk, NULL);
  sphr_k2r(tk, tr, model->ns, NULL);
  return closure(model, res, der, vrsr, cr, tr, prmask, update, damp);
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
    double **tr, double **tk, int *niter)
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
    double **cr, double **vrlr, double **ur, double **tr, double **fr,
    double **ck, double **vklr, double **tk, double **wk,
    const char *fn, int ilam)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt;
  FILE *fp;
  char fnl[80];

  if ( sepout ) {
    sprintf(fnl, "%s%d", fn, ilam);
  } else {
    strcpy(fnl, fn);
  }
  if ((fp = fopen(fnl, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fnl);
    return -1;
  }
  /* print some basic information on the first line */
  fprintf(fp, "# %g %g %g %d\n", 1/(m->kBT*m->beta), m->kBU, m->ampch, m->ietype);
  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;
      for ( l = 0; l < npt; l++ ) {
        double vrl = vrlr[ij][l], vkl = vklr[ij][l];
        fprintf(fp, "%g %g %g %g %g %g %d %d ",
            fft_ri[l], cr[ij][l] - vrl, tr[ij][l] + vrl, fr[ij][l],
            vrl, m->beta * ur[ij][l], i, j);
        if ( printk )
          fprintf(fp, "%g %g %g %g %g",
              fft_ki[l], ck[ij][l] - vkl, tk[ij][l] + vkl, wk[ij][l],
              vkl);
        fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  fprintf(stderr, "saved result to %s\n", fnl);
  return 0;
}



/* compute the equivalent diameter of the molecule
 * this value is only used for comparison
 * and it does not affect the solver */
static double getdiameter(model_t *m)
{
  int i, j, ipr, ns = m->ns, nsv;
  double vol = 0, si, sj, l;

  nsv = getnsv(m);
  for ( i = 0; i < nsv; i++ )
    vol += pow( m->sigma[i], 3 );

  /* deduct the overlap */
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns; j++, ipr++ ) {
      si = m->sigma[i];
      sj = m->sigma[j];
      l = m->dis[ipr];
      if ( l < DBL_MIN ) continue;
      vol -= (si*si*si + sj*sj*sj)/2 - (si*si + sj*sj)*l*3./4
           - pow(si*si - sj*sj, 2)*3./32/l + l*l*l/2;
    }

  return pow(vol, 1./3);
}



static void dorism(model_t *model)
{
  int it, ns, npt, ilam, nlam;
  double err = 0, dia, lam;
  double **ur, **nrdur, **fr, **wk;
  double **cr, **cp, **ck, **tr, **tk;
  double **der, **ntk, **vrlr, **vrsr, **vklr;
  double *um, *mum;

  /* equivalent diameter of the solvent molecule */
  dia = getdiameter(model);

  ns = model->ns;
  npt = model->npt;
  newarr2d(ur,    ns * ns, npt);
  newarr2d(nrdur, ns * ns, npt);
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
  newarr2d(der,   ns * ns, npt);
  newarr2d(ntk,   ns * ns, npt);
  xnew(um, ns);
  xnew(mum, ns);

  initfftw(model->rmax, npt);
  initwk(model, wk);

  /* lambda is used to gradually switch on long-range interaction */
  nlam = model->nlambdas;
  if ( nlam < 1 ) nlam = 1;

  for ( ilam = 1; ilam <= nlam; ilam++ ) {
    lam = 1.*ilam/nlam;

    initfr(model, ur, nrdur, fr, vrlr, vrsr, lam);
    sphr_r2k(vrlr, vklr, ns, NULL);

    /* use f(r) as the initial c(r) for the lowest lambda */
    if ( ilam == 1 )
      cparr2d(cr, fr, ns * ns, npt);

    if ( model->solver == SOLVER_LMV ) {
      err = iter_lmv(model, vrsr, wk, cr, der, ck, vklr, tr, tk, ntk, cp, &it);
    } else if ( model->solver == SOLVER_MDIIS ) {
      err = iter_mdiis(model, vrsr, wk, cr, ck, vklr, tr, tk, &it);
    } else {
      err = iter_picard(model, vrsr, wk, cr, ck, vklr, tr, tk, &it);
    }
    output(model, cr, vrlr, ur, tr, fr, ck, vklr, tk, wk, fncrtr, ilam);
    fprintf(stderr, "lambda %g, %d iterations, err %g, d %g, rho*d^3 %g\n",
        lam, it, err, dia, model->rho[0]*dia*dia*dia);
  }

  calcU(model, ur, cr, tr, fr, um);
  calcmu(model, cr, tr, mum);
  calckirk(model, cr, tr, NULL);
  calccrdnum(model, cr, tr, fr, fncrdnum);

  delarr2d(ur,    ns * ns);
  delarr2d(nrdur, ns * ns);
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
  delarr2d(der,   ns * ns);
  delarr2d(ntk,   ns * ns);
  free(um);
  free(mum);
  donefftw();
}



int main(int argc, char **argv)
{
  model_t *m = doargs(argc, argv);
  dorism(m);
  return 0;
}

