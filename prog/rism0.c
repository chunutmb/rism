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



#define IDBASE 0  /* starting index of site or molecules, 0 or 1 */



#include "debug.h"
#include "util.h"
#include "model.h"
#include "calctd.h"



int model_id = 16;
int verbose = 0;
const char *fncrtr = "out.dat";
const char *fncrtrinp = NULL;
int sepout = 0;
const char *fncrdnum = NULL;
int printk = 0;


/* print help message and die */
static void help(const char *prog)
{
  fprintf(stderr, "Reference site interaction model (RISM)\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options] [input.cfg|model_id]\n\n", prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -o:    followed by the output crtr file, default: %s\n", fncrtr);
  fprintf(stderr, "  -i:    followed by the input crtr file, default: NULL\n");
  fprintf(stderr, "  -k:    print k-space correlation functions, default %d\n", printk);
  fprintf(stderr, "  -$:    separately output file for each lambda, default %d\n", sepout);
  fprintf(stderr, "  -#:    followed by the output coordination number file, default: NULL\n");
  fprintf(stderr, "  -8:    set the magnitude of c(r), default %g\n", CRMAX);
  fprintf(stderr, "  -!:    skip the solute-solute stage calculation, default: false\n");
  fprintf(stderr, "  -u:    always do the solute-solute stage calculation, default: false\n");
  fprintf(stderr, "  -V:    treat all sites as solvent, default: false\n");
  fprintf(stderr, "  -r:    override the density, must be nonnegative, e.g. -d1,0.03 sets the density of the first atom to 0.03\n");
  fprintf(stderr, "  -q:    override the charge, must be nonnegative, e.g. -q2,0.25 sets the charge of the second atom to 0.25\n");
  fprintf(stderr, "  -d:    override the distance, e.g. -d1,2,1.5 sets the distance between the first two atoms to 1.5\n");
  fprintf(stderr, "  -T:    override the temperature, e.g. -T298 sets the temperature to 298K\n");
  fprintf(stderr, "  -C:    override the closure, options: PY, HNC, or KH\n");
  fprintf(stderr, "  -S:    override the solver, options: Picard, LMV, or MDIIS\n");
  fprintf(stderr, "  -^:    override the number of lambda stages\n");
  fprintf(stderr, "  -v:    be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h:    display this message\n");
  exit(1);
}



/* handle command line arguments */
static model_t *doargs(int argc, char **argv)
{
  model_t *m;
  int i, j, ch;
  const char *fncfg = NULL, *prog = argv[0];
  char *p, *q;

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
      if ( strchr("oi#rqdTCS^8", ch) != NULL ) {
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
          help(prog);
        }
        if ( ch == 'o' ) {
          fncrtr = q;
        } else if ( ch == 'i' ) {
          fncrtrinp = q;
        } else if ( ch == '#' ) {
          fncrdnum = q;
        } else if ( ch == 'r' ) { /* override the density */
          if ( model_register_arr(model_usr->rho, q) != 0 )
            help(prog);
        } else if ( ch == 'q' ) { /* override the charge */
          if ( model_register_arr(model_usr->charge, q) != 0 )
            help(prog);
        } else if ( ch == 'd' ) { /* override the distance of a bond */
          if ( model_register_disij(model_usr, q) != 0 )
            help(prog);
        } else if ( ch == 'T' ) {
          model_usr->beta = 1/atof(q);
        } else if ( ch == 'C' ) { /* override the closure */
          model_usr->ietype = model_select(q, IE_COUNT, ietypes);
        } else if ( ch == 'S' ) { /* override the solver */
          model_usr->solver = model_select(q, SOLVER_COUNT, solvers);
        } else if ( ch == '^' ) { /* override the number of lambdas */
          model_usr->nlambdas = atoi(q);
        } else if ( ch == '8' ) { /* override the maximal value of c(r) */
          model_usr->crmax = atof(q);
        }
        break; /* skip the rest of the characters in the option */
      } else if ( ch == '!' ) {
        model_usr->douu = DOUU_NEVER;
      } else if ( ch == 'u' ) {
        model_usr->douu = DOUU_ALWAYS;
      } else if ( ch == 'V' ) {
        model_usr->douu = DOUU_ALLSOLVENT;
      } else if ( ch == 'k' ) {
        printk = 1;
      } else if ( ch == '$' ) {
        sepout = 1;
      } else if ( ch == 'v' ) {
        verbose++;
      } else if ( ch == 'h' ) {
        help(prog);
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        help(prog);
      }
    }
  }

  m = models + model_id;
  if ( fncfg != NULL )
    if ( model_load(m, fncfg, verbose) != 0 ) {
      fprintf(stderr, "failed to load %s\n", fncfg);
      help(prog);
    }
  /* override options from the command-line */
  model_override(m, model_usr);
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

  if ( c6 < 0 && c12 > 0 ) { /* attractive */
    sig = pow(-c12/c6, 1./6);
    eps = c6*c6/4/c12;
    return ljpot(r, sig, eps, eps, lam, nrdu);
  }
  if ( c6 > 0 ) lam = 1;
  ir = 1/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  *nrdu = ir6 * (12 * c12 * ir6 + 6 * c6 * lam);
  return ir6 * (c12 * ir6 + c6 * lam);
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



/* Huggins-Mayer potential: B * exp(-r/rho) - C/r^6
 * which is used in
 * Alkali halides in water: Ion-solvent correlations and ion-ion potentials of
 * mean force at infinite dilution
 * B. Montgomery Pettitt and Peter J. Rossky
 * J. Chem. Phys. 84(10) 5836-5844 (1986) */
static double HMpot(double r, double B, double C, double rho,
    double lam, double *nrdu)
{
  double ir, ir6, r1;

  ir = 1/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;

  r1 = r/rho;
  B *= exp(-r1);

  *nrdu = B * r1 - 6 * C * ir6 * lam;
  return B - C*ir6;
}



/* find the peak radius of the short-range part of the Huggins-Mayer potential:
 *  B exp(-r/rho) - C/r^6 */
static double solveHMrmin(double B, double C, double rho)
{
  double r, rp;
  int i;
  const int itmax = 1000;

  if ( B <= 0 || C <= 0 ) return 0;
  rp = rho;
  for ( i = 0; i < itmax; i++ ) {
    r = pow(6*C*rho/B*exp(rp/rho), 1./7);
    if (fabs(r - rp) <  1e-10) break;
    rp = r;
  }
  return r;
}



/* initialize f(r) */
static void initfr(model_t *m, double **ur, double **nrdur,
    double **fr, double **vrlj, double **vrqq, double **vrpp,
    double **vrlr, double **vrsr, double lam)
{
  const double vrmin = -30;
  int i, j, ij, ji, ipr, l, ns = m->ns, npt = m->npt, use_pairpot;
  double beta = m->beta, r, z, u, uelec, nrdu = 0, ulr, rscreen;
  double sig, eps6, eps12, c6, c12, Bij, rhoij, rminij;
  static int once;

  getmols(m); /* parse the system into molecules */
  for ( ipr = 0, i = 0; i < ns; i++ ) { /* the first site */
    for ( j = i; j < ns; j++, ipr++ ) { /* the second site */
      ij = i*ns + j;
      ji = j*ns + i;

      sig = m->pairpot[ipr].sigma;
      eps6 = m->pairpot[ipr].eps6;
      eps12 = m->pairpot[ipr].eps12;
      c6 = m->pairpot[ipr].C6;
      c12 = m->pairpot[ipr].C12;
      Bij = m->pairpot[ipr].B;
      rhoij = m->pairpot[ipr].rho;
      rminij = solveHMrmin(Bij, -c6, rhoij);
      /* if any of the pair energy is set, we use the pair potential */
      use_pairpot = ( fabs(eps6) > DBL_MIN || fabs(eps12) > DBL_MIN
          || fabs(c6) > DBL_MIN || fabs(c12) > DBL_MIN || fabs(Bij) > DBL_MIN );
      if ( !use_pairpot ) { /* otherwise, we use the starndard LJ form */
        sig = .5 * (m->sigma[i] + m->sigma[j]);
        eps6 = sqrt(m->eps6_12[i][0] * m->eps6_12[j][0]);
        eps12 = sqrt(m->eps6_12[i][1] * m->eps6_12[j][1]);
      } else {
        if ( fabs(Bij) > 0 ) sig = rhoij;
      }
      rscreen = m->pairpot[ipr].rscreen;
      if ( rscreen <= 0 ) rscreen = m->rscreen * sqrt(2);
      if ( !once )
        fprintf(stderr, "i %d, j %d, sig %4.2f, eps6 %10.2e, eps12 %10.2e, "
            "c6 %10.2e c12 %10.2e, B %10.2e, rho %4.2f, rmin %g, rscreen %g\n",
          i + IDBASE, j + IDBASE, sig, eps6, eps12, c6, c12, Bij, rhoij, rminij, rscreen);

      for ( l = 0; l < npt; l++ ) { /* the radius */
        r = fft_ri[l];
        if ( m->ljtype  == HARD_SPHERE) {
          if (r < sig) {
            vrsr[ij][l] = 2*INFTY;
            z = -1;
          } else {
            vrsr[ij][l] = 0;
            z = 0;
          }
          vrlj[ij][l] = vrsr[ij][l];
          nrdur[ij][l] = ur[ij][l] = 0;
          vrlr[ij][l] = vrqq[ij][l] = vrpp[ij][l] = 0;
        } else { /* Lennard-Jones */
          if ( use_pairpot ) {
            if ( Bij > 0 ) {
              /* to avoid blow-up, the minimal r is rminij */
              u = HMpot((r < rminij ? rminij : r), Bij, -c6, rhoij, lam, &nrdu);
            } else if ( fabs(c6) > DBL_MIN || fabs(c12) > DBL_MIN ) {
              u = ljpot6_12(r, c6, c12, lam, &nrdu);
            } else {
              u = ljpot(r, sig, eps6, eps12, lam, &nrdu);
            }
          } else if ( m->ljtype == LJ_REPULSIVE ) {
            u = ljrpot(r, sig, eps6, eps12, &nrdu);
          } else {
            u = ljpot(r, sig, eps6, eps12, lam, &nrdu);
          }
          uelec = lam * m->ampch * m->charge[i] * m->charge[j] / r;
          ur[ij][l] = u + uelec;
          nrdur[ij][l] = nrdu + uelec;
          ulr = uelec * erf( r/rscreen );
          vrlj[ij][l] = beta * u;
          vrqq[ij][l] = beta * uelec;
          /* vrpp concerns only the net charge on the molecule */
          vrpp[ij][l] = beta * lam * m->ampch / r
                      * m->chargemol[ m->mol[i] ] * m->chargemol[ m->mol[j] ];
          vrlr[ij][l] = beta * ulr;
          vrsr[ij][l] = beta * (u + uelec - ulr);
          if ( vrsr[ij][l] < vrmin ) {
            fprintf(stderr, "vrsr(%d-%d, %g) = %g, too negative, may decrease rscreen\n",
                i, j, r, vrsr[ij][l]);
            vrsr[ij][l] = vrmin;
          }

          z = exp(-vrsr[ij][l]) - 1;
        }
        fr[ij][l] = z;
      } /* loop over l, the radius */
      if ( j > i ) {
        cparr(ur[ji],     ur[ij],     npt);
        cparr(nrdur[ji],  nrdur[ij],  npt);
        cparr(vrlj[ji],   vrlj[ij],   npt);
        cparr(vrqq[ji],   vrqq[ij],   npt);
        cparr(vrpp[ji],   vrpp[ij],   npt);
        cparr(vrlr[ji],   vrlr[ij],   npt);
        cparr(vrsr[ji],   vrsr[ij],   npt);
        cparr(fr[ji],     fr[ij],     npt);
      }
    } /* loop over j, the second site */
  } /* loop over i, the first site */
  once = 1;
}



/* initialize the w matrix for intra-molecular covalent bonds */
static void initwk(model_t *m, double **wk)
{
  int i, j, ipr, u, ns = m->ns, npt = m->npt;

  for ( ipr = 0, i = 0; i < ns; i++ ) {
    for ( u = 0; u < npt; u++ ) // diagonal
      wk[i*ns + i][u] = 1;
    for ( j = i + 1; j < ns; j++, ipr++ ) {
      int ij = i*ns + j;
      double l = m->dis[ipr];
      if ( l <= 0 ) continue;
      for ( u = 0; u < npt; u++ ) {
        double kl = fft_ki[u] * l;
        wk[ij][u] = sin(kl)/kl;
      }
      cparr(wk[j*ns + i], wk[ij], npt);
    }
  }
}



/* Ornstein-Zernike relation: c(k) --> t(k) */
static void oz(model_t *m, double **ck, double **vklr,
    double **tk, double **wk, double **invwc1w)
{
  int i, j, ij, l, ns = m->ns;
  double w[NS2MAX], dw[NS2MAX], c[NS2MAX], wc[NS2MAX], invwc1[NS2MAX], v[NS2MAX] = {0};
  double tm1[NS2MAX], tm2[NS2MAX], tm3[NS2MAX];

  for ( l = 0; l < m->npt; l++ ) {
    for ( ij = 0; ij < ns * ns; ij++ ) {
      if ( vklr != NULL ) v[ij] = vklr[ij][l];
      w[ij] = wk[ij][l];
      c[ij] = ck[ij][l] - v[ij];
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

    invmat(tm2, invwc1, ns); /* invwc1 = (1 - rho w c)^(-1) */

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
      tk[ij][l] = tm1[ij] + tm2[ij] + tm3[ij] - v[ij];
  }
}



/* compute the chain sum Q(r) for the charge interation */
static void getQr(model_t *m, double **vrqq, double **Qr,
    double **Qk, double **wk, double **vkqq,
    double **Qrx, double **vrlj)
{
  int ij, l, ns = m->ns, npt = m->npt;

  /* c(r) --> -vqq(r) for large r */
  sphr_r2k(vrqq, vkqq, ns, NULL);
  for ( l = 0; l < npt; l++ )
    for ( ij = 0; ij < ns * ns; ij++ )
      vkqq[ij][l] = -vkqq[ij][l];
  oz(m, vkqq, NULL, Qk, wk, NULL);
  sphr_k2r(Qk, Qr, ns, NULL);
  for ( ij = 0; ij < ns * ns; ij++ )
    for ( l = 0; l < npt; l++ )
      Qrx[ij][l] = Qr[ij][l] - vrlj[ij][l];
}



/* return the updated cr */
static double getcr(double tr, double vrsr, double Qrx,
    double *dcr, int ietype, double crmax)
{
  double xp, del, cr = 0;

  del = -vrsr + tr;
  if ( ietype == IE_HNC ) {
    xp = exp(del);
    cr = xp - 1 - tr;
    if ( dcr != NULL ) *dcr = xp - 1;
  } else if ( ietype == IE_PY ) {
#ifdef PYA
    /* currently, the PYA closure does not work
     * this branch is disabled by default */
    /* del == -vrlj + Qr + taur == Qrx + taur */
    double taur = tr - vrsr - Qrx;
#else
    double taur = tr;
    (void) Qrx;
#endif
    xp = exp(del - taur);
    cr = xp * (1 + taur) - 1 - tr;
    if ( dcr != NULL ) *dcr = xp - 1;
  } else if ( ietype == IE_KH ) {
    if ( del <= 0 ) { /* HNC */
      xp = exp(del);
      cr = xp - 1 - tr;
      if ( dcr != NULL ) *dcr = xp - 1;
    } else {
      cr = -vrsr;
      if ( dcr != NULL ) *dcr = 0;
    }
  } else {
    fprintf(stderr, "unknown closure %d\n", ietype);
    exit(1);
  }
  if ( cr > crmax ) cr = crmax;
  else if ( cr < -crmax ) cr = -crmax;
  return cr;
}



struct { int l; double err; } errspot[NS2MAX];

/* apply the closure
 * compute residue vector if needed */
static double closure(model_t *model,
    double *res, double **der, double **vrsr,
    double **cr, double **tr, double **Qrx,
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
        y = getcr(tr[ij][l], vrsr[ij][l], Qrx[ij][l],
            der ? &der[ij][l] : NULL, model->ietype, model->crmax)
          - cr[ij][l];
        if ( res != NULL ) res[id] = y;
        if ( update ) cr[ij][l] += damp * y;
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
      if ( j > i ) cparr(cr[ji], cr[ij], npt);
    }
  return errm;
}




/* a step of direct iteration (Picard)
 * compute residue vector if needed */
static double step_picard(model_t *model,
    double *res, double **der,
    double **vrsr, double **wk,
    double **cr, double **ck, double **vklr,
    double **tr, double **tk, double **Qrx,
    int *prmask, int update, double damp)
{
  sphr_r2k(cr, ck, model->ns, NULL);
  oz(model, ck, vklr, tk, wk, NULL);
  sphr_k2r(tk, tr, model->ns, NULL);
  return closure(model, res, der, vrsr, cr, tr, Qrx,
                 prmask, update, damp);
}



/* the uu step for infinitely diluted atomic solute */
static double step_uu_infdil_atomicsolute(model_t *model,
    double **vrsr, double **wk, double **cr, double **ck,
    double **vklr, double **tr, double **tk, double **Qrx,
    int *prmask)
{
  return step_picard(model, NULL, NULL, vrsr, wk, cr, ck,
      vklr, tr, tk, Qrx, prmask, 1, 1.0);
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
    double **tr, double **tk, double **Qrx,
    int *niter)
{
  int it;
  double err = 0, errp = errinf;

  for ( it = 0; it < model->itmax; it++ ) {
    err = step_picard(model, NULL, NULL, vrsr, wk, cr, ck, vklr,
        tr, tk, Qrx, NULL, 1, model->picard.damp);
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
static char *output(model_t *m,
    double **cr, double **vrqq, double **vrpp, double **vrlr,
    double **ur, double **tr, double **fr,
    double **ck, double **vklr, double **tk, double **wk,
    double **Qr, double **Qrx, double **Qk,
    const char *fn, int ilam)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt;
  FILE *fp;
  static char fnl[FILENAME_MAX];
  double eps_rism;

  if ( sepout ) {
    sprintf(fnl, "%s%d", fn, ilam);
  } else {
    strcpy(fnl, fn);
  }
  if ((fp = fopen(fnl, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fnl);
    return NULL;
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
        double vrt = m->beta * ur[ij][l], vrq = vrqq[ij][l], vrp = vrpp[ij][l];
        /* note that cr[ij][l] and tr[ij][l] exclude
         * the long-range component vrl */
        double crt = cr[ij][l] - vrl, trt = tr[ij][l] + vrl;
        double ckt = ck[ij][l] - vkl, tkt = tk[ij][l] + vkl;
        /* pmfs = beta dW: the short-range correction to
         * the continuum primitive model, in which
         *
         *  W_c = u_LJ + u_pp/eps  (with u = u_LJ + u_pp)
         *
         *  Here u_pp is net molecular electrostatic interaction
         *
         * `pmfs' is computed from Eq. (50) of the following paper:
         * [ ``The interionic potential of mean force in a molecular polar
         *     solvent from an extended RISM equation''
         *    Hirate, Rossky, and Pettitt,
         *    J. Chem. Phys. 78(6) 4133-4144 (1983) ]
         *
         *    beta dW
         *  = beta W_rism_corrected - beta W_c
         *  = (beta W_s - phi_pp/eps) - beta W_c
         *  = beta W_s + phi*                 (here, phi* = -beta u_LJ)
         *  = beta W + phi/eps_rism + phi*    (here, phi  = -beta u_pp)
         *  = beta u - t(r) - beta u_pp/eps_rism - beta u_LJ
         *  = beta u_q - t(r) - beta u_pp/eps_rism
         * */
        double pmfs = -getlngr(m, crt, trt, vrt) - vrt + vrq - vrp/eps_rism;
        fprintf(fp, "%g %g %g %g %g %g %d %d %g %g %g %g %g ",
            fft_ri[l], crt, trt, fr[ij][l],
            vrl, vrt, i + IDBASE, j + IDBASE,
            vrq, pmfs, vrp, Qr[ij][l], Qrx[ij][l]);
        if ( printk )
          fprintf(fp, "%g %g %g %g %g %g",
              fft_ki[l], ckt, tkt, wk[ij][l], vkl, Qk[ij][l]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return fnl;
}



/* load correlation functions from the file `fn' */
static int loadcrtr(model_t *m, double **cr, double **tr,
    const char *fn)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt;
  double ri, crt, trt, fr, vrl;
  char ln[40960];
  FILE *fp;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }
  if ( fgets(ln, sizeof ln, fp) == NULL
    || ln[0] != '#' ) {
    fprintf(stderr, "invalid file %s\n", fn);
    return -1;
  }
  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;
      for ( l = 0; l < npt; l++ ) {
        do {
          if ( fgets(ln, sizeof ln, fp) == NULL ) {
            fprintf(stderr, "%s stops at i %d, j %d, l %d\n",
                fn, i, j, l);
            return -1;
          }
          strstrip(ln);
        } while ( ln[0] == '\0' ); /* skip blank lines */

        if ( 5 != sscanf(ln, "%lf %lf %lf %lf %lf",
              &ri, &crt, &trt, &fr, &vrl) ) {
          fprintf(stderr, "%s: corrupted on line i %d j %d l %d\n%s",
              fn, i, j, l, ln);
          return -1;
        }
        cr[ij][l] = crt + vrl;
        tr[ij][l] = trt - vrl;
      }
      if ( j > i ) {
        int ji = j*ns + i;
        cparr(cr[ji], cr[ij], npt);
        cparr(tr[ji], tr[ij], npt);
      }
    }
  }
  fclose(fp);
  fprintf(stderr, "successfully loaded %s\n", fn);
  return 0;
}



static int dorism(model_t *model)
{
  int it, ns, npt, ilam, nlam;
  double err = 0, dia, eps, lam;
  double **ur, **nrdur, **fr, **wk;
  double **cr, **cp, **ck, **tr, **tk, **ntk;
  double **der, **vrlr, **vrlj, **vrqq, **vrpp, **vrsr, **vklr;
  double **Qr, **Qk, **Qrx;
  double *um, *mum;
  const char *fnout;
  uv_t *uv;

  /* equivalent diameter of the solvent molecule */
  dia = getdiameters(model);

  ns = model->ns;
  npt = model->npt;
  newarr2d(ur,    ns * ns, npt);
  newarr2d(nrdur, ns * ns, npt);
  newarr2d(vrlj,  ns * ns, npt);
  newarr2d(vrqq,  ns * ns, npt);
  newarr2d(vrpp,  ns * ns, npt);
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
  newarr2d(Qr,    ns * ns, npt);
  newarr2d(Qk,    ns * ns, npt);
  newarr2d(Qrx,   ns * ns, npt);
  xnew(um, ns);
  xnew(mum, ns);

  initfft(model->rmax, npt);
  initwk(model, wk);

  /* lambda is used to gradually switch on long-range interaction */
  nlam = model->nlambdas;
  if ( nlam < 1 ) nlam = 1;
  /* if we use an input file, don't use lambda */
  if ( fncrtrinp != NULL ) nlam = 1;

  for ( ilam = 1; ilam <= nlam; ilam++ ) {
    lam = 1.*ilam/nlam;

    initfr(model, ur, nrdur, fr, vrlj, vrqq, vrpp, vrlr, vrsr, lam);
    sphr_r2k(vrlr, vklr, ns, NULL);
    getQr(model, vrqq, Qr, Qk, wk, ntk, Qrx, vrlj);

    if ( ilam == 1 ) {
      if ( fncrtrinp != NULL ) {
        /* load previous c(r) as the initial condition */
        if ( loadcrtr(model, cr, tr, fncrtrinp) != 0 ) {
          fprintf(stderr, "error loading %s\n", fncrtrinp);
          goto EXIT;
        }
        output(model, cr, vrqq, vrpp, vrlr, ur, tr, fr,
          ck, vklr, tk, wk, Qr, Qrx, Qk, "test.dat", ilam);
      } else {
        /* use f(r) as the initial c(r) for the lowest lambda */
        cparr2d(cr, fr, ns * ns, npt);
      }
    }

    /* initialize the manager for solute-solvent iteraction */
    uv = uv_open(model);

    if ( model->solver == SOLVER_LMV ) {
      err = iter_lmv(model, vrsr, wk, cr, der, ck, vklr,
          tr, tk, ntk, Qrx, cp, uv, &it);
    } else if ( model->solver == SOLVER_MDIIS ) {
      err = iter_mdiis(model, vrsr, wk, cr, ck, vklr,
          tr, tk, Qrx, uv, &it);
    } else {
      err = iter_picard(model, vrsr, wk, cr, ck, vklr,
          tr, tk, Qrx, &it);
    }

    uv_close(uv);

    fnout = output(model, cr, vrqq, vrpp, vrlr, ur, tr, fr,
        ck, vklr, tk, wk, Qr, Qrx, Qk, fncrtr, ilam);
    fprintf(stderr, "lambda %4.2f, %4d iterations, err %10.4e, output %s\n",
        lam, it, err, fnout);
  }

  calcU(model, ur, cr, tr, vrsr, um);
  calcchempot(model, cr, tr, vrsr, vrlr, mum, verbose);
  calckirk(model, cr, tr, vrsr, NULL);
  calccrdnum(model, cr, tr, vrsr, fncrdnum);
  eps = calcdielec(model);
  printf("dielectric constant %g, d %g, rho*d^3 %g\n",
      eps, dia, model->rho[0]*dia*dia*dia);

EXIT:
  delarr2d(ur,    ns * ns);
  delarr2d(nrdur, ns * ns);
  delarr2d(vrlj,  ns * ns);
  delarr2d(vrqq,  ns * ns);
  delarr2d(vrpp,  ns * ns);
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
  delarr2d(Qr,    ns * ns);
  delarr2d(Qk,    ns * ns);
  delarr2d(Qrx,   ns * ns);
  free(um);
  free(mum);
  donefft();
  return 0;
}



int main(int argc, char **argv)
{
  model_t *m = doargs(argc, argv);
  return dorism(m);
}

