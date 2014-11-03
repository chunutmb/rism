#ifndef MODEL_H__
#define MODEL_H__



/* input parameters */



#define MAXATOM 6 /* maximal number of atoms/sites */

enum { HARD_SPHERE, LJ_FULL, LJ_REPULSIVE };
enum { IE_PY, IE_HNC };
enum { SOLVER_PICARD, SOLVER_LMV, SOLVER_MDIIS };

typedef struct {
  double damp;
} picard_params_t;

typedef struct {
  int M; /* number of points for the Newton-Raphson method */
  double damp;
} lmv_params_t;

typedef struct {
  int nbases; /* number of bases */
  double damp;
} mdiis_params_t;



typedef struct {
  int ns;
  double sigma[MAXATOM];
  double eps6_12[MAXATOM][2];
  double C6_12[MAXATOM*(MAXATOM+1)/2][2]; /* alternative to sigma/epsilon */
  double rho[MAXATOM];
  double dis[MAXATOM*(MAXATOM-1)/2];
  double beta;
  int ljtype;

  double charge[MAXATOM];
  double ampch;
  double rscreen; /* screening distance */

  int ietype;
  double rmax; /* radius cutoff */
  int npt; /* number of sampling points along r */

  int nlambdas; /* number of intermediate stages */
  int itmax; /* maximial number of iterations in each stage */
  double tol; /* tolerance of error */
  int solver; /* solver */

  picard_params_t picard;
  lmv_params_t    lmv;
  mdiis_params_t  mdiis;
} model_t;



/* return the index of string from a predefined array */
static int model_select(const char *s, int n, const char **arr)
{
  int i;

  for ( i = 0; i < n; i++ )
    if ( strcmp(arr[i], s) == 0 )
      return i;
  fprintf(stderr, "cannot find %s\n", s);
  return 0;
}



/* return the index of an array index */
static int model_getidx(char *s, int n)
{
  char *p = strchr(s, '(');
  int i;

  if ( p == NULL ) p = strchr(s, '[');
  if ( p == NULL ) {
    fprintf(stderr, "cannot find the index of %s\n", s);
    exit(1);
  }
  i = atoi(p + 1) - 1;
  if ( i >= n ) {
    fprintf(stderr, "bad index for %s, i %d >= %d\n", s, i, n);
    exit(1);
  }
  return i;
}



/* return the index of an array index */
static int model_getidx2(char *s, int *j, int n)
{
  char *p = strchr(s, '('), *q;
  int i, k;

  if ( p == NULL ) p = strchr(s, '[');
  if ( p == NULL ) {
    fprintf(stderr, "cannot find the index of %s\n", s);
    exit(1);
  }
  p++;
  q = strchr(p, ',');
  if ( q == NULL ) {
    fprintf(stderr, "cannot find the second index of %s\n", s);
    exit(1);
  }
  *q++ = '\0';
  i = atoi(strstrip(p)) - 1;
  *j = atoi(strstrip(q)) - 1;
  if ( i > *j ) k = i, i = *j, *j = k;
  if ( i >= n || *j >= n ) {
    fprintf(stderr, "bad index for %s, i %d or j %d >= %d\n",
        s, i, *j, n);
    exit(1);
  }
  return i;
}



/* load model from file `fn' */
static int model_load(model_t *m, const char *fn, int verbose)
{
  FILE *fp;
  char buf[800], *p, *key, *val;
  int i, j, ipr, ns = -1, inpar;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  while ( fgets(buf, sizeof buf, fp) ) {
    strstrip(buf); /* remove trailing spaces */
    if ( buf[0] == '\0' || buf[0] == '#' ) continue;

    /* parse the line to a key and a value */
    /* find the end of the key */
    inpar = 0; /* within (...) */
    for ( p = buf; *p; p++ ) {
      if ( (!inpar && isspace(*p)) || *p == '=' ) {
        *p = '\0'; /* end the key part */
        break;
      }
      if ( !inpar && (*p == '(' || *p == '[') )
        inpar = 1; /* enter a parentheses block */
      else if ( inpar && (*p == ')' || *p == ']') )
        inpar = 0; /* leave a parentheses block */
      *p = (char) tolower(*p);
    }
    key = buf;

    /* find the beginning of the value */
    for ( p++; isspace(*p) || *p == '=' ; ) p++;
    val = p;
    for ( ; *p; p++ ) *p = (char) tolower(*p);

    //printf("key: %-40s val: %s\n", key, val);

#define CHECK_NS(key) if ( ns <= 0 ) { \
      fprintf(stderr, "must set `ns' before `%s'\n", key); \
      exit(1); }

    if ( strcmp(key, "ns") == 0 ) {
      m->ns = ns = atoi(val);
    } else if ( strncmp(key, "sigma", 5) == 0 ) {
      CHECK_NS(key);
      i = model_getidx(key, ns);
      m->sigma[i] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "sigma(%d)     = %g\n", i, m->sigma[i]);
    } else if ( strncmp(key, "eps6", 4) == 0 ) {
      CHECK_NS(key);
      i = model_getidx(key, ns);
      m->eps6_12[i][0] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "eps6(%d)      = %g\n", i, m->eps6_12[i][0]);
    } else if ( strncmp(key, "eps12", 5) == 0 ) {
      CHECK_NS(key);
      i = model_getidx(key, ns);
      m->eps6_12[i][1] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "eps12(%d)     = %g\n", i, m->eps6_12[i][1]);
    } else if ( strncmp(key, "rho", 3) == 0 ) {
      CHECK_NS(key);
      i = model_getidx(key, ns);
      m->rho[i] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "rho(%d)       = %g\n", i, m->rho[i]);
    } else if ( strncmp(key, "charge", 6) == 0 ) {
      CHECK_NS(key);
      i = model_getidx(key, ns);
      m->charge[i] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "charge(%d)    = %g\n", i, m->charge[i]);
    } else if ( strncmp(key, "dis", 3) == 0 ) {
      CHECK_NS(key);
      i = model_getidx2(key, &j, ns);
      ipr = ns*i - (i+1)*(i+2)/2 + j;
      m->dis[ipr] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "dis(%d, %d)    = %g (pair %d)\n", i, j, m->dis[ipr], ipr);
    } else if ( strncmp(key, "c6", 2) == 0 ) {
      CHECK_NS(key);
      i = model_getidx2(key, &j, ns);
      ipr = ns*i - i*(i+1)/2 + j;
      m->C6_12[ipr][0] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "c6(%d, %d)     = %g (pair %d)\n", i, j, m->C6_12[ipr][0], ipr);
    } else if ( strncmp(key, "c12", 3) == 0 ) {
      CHECK_NS(key);
      i = model_getidx2(key, &j, ns);
      ipr = ns*i - i*(i+1)/2 + j;
      m->C6_12[ipr][1] = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "c12(%d, %d)    = %g (pair %d)\n", i, j, m->C6_12[ipr][1], ipr);
    } else if ( strcmp(key, "t") == 0 || strcmp(key, "temp") == 0 ) {
      double temp = atof(val);
      m->beta = 1/temp;
      if ( verbose >= 2 )
        fprintf(stderr, "T            = %g\n", temp);
    } else if ( strcmp(key, "ljtype") == 0 ) {
      const char *ljtypes[3];
      ljtypes[HARD_SPHERE] = "hard-sphere";
      ljtypes[LJ_FULL] = "lj-full";
      ljtypes[LJ_REPULSIVE] = "lj-repulsive";
      m->ljtype = model_select(val, 3, ljtypes);
      if ( verbose >= 2 )
        fprintf(stderr, "LJ type      = %s\n", ljtypes[m->ljtype]);
    } else if ( strcmp(key, "rscreen") == 0 ) {
      m->rscreen = atof(val);
    } else if ( strcmp(key, "closure") == 0
             || strcmp(key, "ietype") == 0 ) {
      const char *ietypes[2];
      ietypes[IE_PY] = "py";
      ietypes[IE_HNC] = "hnc";
      m->ietype = model_select(val, 2, ietypes);
      if ( verbose >= 2 )
        fprintf(stderr, "closure      = %s\n", ietypes[m->ietype]);
    } else if ( strcmp(key, "rmax") == 0 ) {
      m->rmax = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "rmax         = %g\n", m->rmax);
    } else if ( strcmp(key, "npt") == 0
             || strcmp(key, "n-pts") == 0 ) {
      m->npt = atoi(val);
      if ( verbose >= 2 )
        fprintf(stderr, "npt          = %d\n", m->npt);
    } else if ( strncmp(key, "nlambda", 6) == 0 ) {
      m->nlambdas = atoi(val);
      if ( verbose >= 2 )
        fprintf(stderr, "nlambdas     = %d\n", m->nlambdas);
    } else if ( strcmp(key, "tol") == 0 ) {
      m->tol = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "tol          = %g\n", m->tol);
    } else if ( strcmp(key, "solver") == 0 ) {
      const char *solvers[3];
      solvers[SOLVER_PICARD] = "picard";
      solvers[SOLVER_LMV] = "lmv";
      solvers[SOLVER_MDIIS] = "mdiis";
      m->solver = model_select(val, 3, solvers);
      if ( verbose >= 2 )
        fprintf(stderr, "solver       = %s\n", solvers[m->solver]);
    } else if ( strcmp(key, "picard_damp") == 0 ) {
      m->picard.damp = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "Picard_damp  = %g\n", m->picard.damp);
    } else if ( strcmp(key, "lmv_m") == 0 ) {
      m->lmv.M = atoi(val);
      if ( verbose >= 2 )
        fprintf(stderr, "LMV_M        = %d\n", m->lmv.M);
    } else if ( strcmp(key, "lmv_damp") == 0 ) {
      m->lmv.damp = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "LMV_damp     = %g\n", m->lmv.damp);
    } else if ( strcmp(key, "mdiis_nbases") == 0 ) {
      m->mdiis.nbases = atoi(val);
      if ( verbose >= 2 )
        fprintf(stderr, "MDIIS_nbases = %d\n", m->mdiis.nbases);
    } else if ( strcmp(key, "mdiis_damp") == 0 ) {
      m->mdiis.damp = atof(val);
      if ( verbose >= 2 )
        fprintf(stderr, "MDIIS_damp   = %g\n", m->mdiis.damp);
    }
  }
  if ( verbose >= 2 ) {
    fprintf(stderr, "press Enter to start...");
    getchar();
  }

  fclose(fp);

  if ( m->ns < 1 ) {
    fprintf(stderr, "invalid number of sites %d\n", m->ns);
    return -1;
  }
  return 0;
}



/* built-in models from literature
 *
 * References:
 *
 * Solution of a new integral equation for pair correlation function in molecular liquids
 * Lawrence J. Lowden and David Chandler
 * J. Chem. Phys. 59 (12) 6587 (1973)
 *
 * Applications of the RISM equation to diatomic fluids:
 * the liquids nitrogen, oxygen and bromine
 * C.S. Hsu, David Chandler and L.J. Lowden
 * Chem. Phys. 14 213-228 (1976)
 *
 * Comparisons of Monte Carlo and RISM calculations of pair correlation functions
 * David Chandler, C. S. Hsu and William B. Streett
 * J. Chem. Phys. 66 (11) 5231 (1977)
 *
 * Computation of the correlation functions for fluids composed of
 * diatomic molecules by means of the method of integration equations
 * Kazumitsu Kojima and Kiyoshi Arakawa
 * Bulletin of the Chemical Society of Japan 51(7) 1977-1981 (1978)
 *
 * An extended RISM equation for molecular polar fluids
 * Fumio Hirata and Peter J. Rossky
 * Chem. Phys. Lett. 83(2) 329-334 (1981)
 *
 * Application of an extended RISM equation to dipolar and quadrupolar fluids
 * Fumio Hirata, B. Montgomery Pettitt, Peter J. Rossky
 * J. Chem. Phys. 77(1) 509-520 (1982)
 *
 * Integral equation prediction of liquid state structure for
 * waterlike intermolecular potentials
 * B. Montgomery Pettitt and Peter J. Rossky
 * J. Chem. Phys. 77(3) 1451-1457 (1982)
 *
 * The interionic potential of mean force in a molecular polar
 * solvent from an extended RISM equation
 * Fumio Hirata, Peter J. Rossky, and B. Montgomery Pettitt
 * J. Chem. Phys. 78(6) 4133-4144 (1983) */
model_t models[] =
{
  /* 0. default model */
  {0, {0}, {{0}}, {{0}},
    {0}, {0}, 1./300, LJ_FULL,
    {0}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {20, 0.5},
    {5, 0.5}
  },
  /* 1. LC 1973, and Model I of CHS 1977 */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.500, 0.500}, {0.600}, 1.000, HARD_SPHERE,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_MDIIS,
    1.0,
    0, 1.0,
    5, 1.0
  },
  /* 2. Model II of CHS 1977 */
  {2, {0.790, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.686, 0.686}, {0.490}, 1.000, HARD_SPHERE,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.5},
    {0, 1.0},
    {5, 1.0}
  },
  /* 3. Model III of CHS 1977 */
  {2, {0.675, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.825, 0.825}, {0.346}, 1.000, HARD_SPHERE,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.3},
    {0, 0.5},
    {5, 1.0}
  },
  /* 4. LC1973, liquid nitrogen */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.696, 0.696}, {1.1/3.341}, 1/1.83, LJ_REPULSIVE,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 10000, 1e-7,
    SOLVER_LMV,
    {0.15},
    {0, 0.5},
    {3, 0.02}
  },
  /* 5. KA1978, liquid nitrogen */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.6964, 0.6964}, {1.1/3.341}, 1/1.61, LJ_FULL,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    5, 10000, 1e-7,
    SOLVER_MDIIS,
    {0.01}, /* does not work */
    {0, 0.5},
    {5, 0.2}
  },
  /* 6. HR1981, liquid nitrogen, neutral */
  {2, {3.341, 3.341}, {{1, 1}, {1, 1}}, {{0}},
    {0.01867, 0.01867}, {1.1}, 1/1.636, LJ_FULL,
    {0}, 0, 0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {10, 0.3},
    {5, 0.1}
  },
  /* 7. HR1981, liquid nitrogen, charged, also HPR1982, model I */
  {2, {3.341, 3.341}, {{44.0, 44.0}, {44.0, 44.0}}, {{0}},
    {0.01867, 0.01867}, {1.1}, 1./72, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {10, 0.5},
    {5, 0.1} /* does not work! */
  },
  /* 8. HPR1982, HCl, model II */
  {2, {2.735, 3.353}, {{20.0, 20.0}, {259.0, 259.0}}, {{0}},
    {0.018, 0.018}, {1.257}, 1./210, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {10, 0.5},
    {5, 0.1}
  },
  /* 9. HPR1982, HCl, model III */
  {2, {0.4, 3.353}, {{20.0, 20.0}, {259.0, 259.0}}, {{0}},
    {0.018, 0.018}, {1.3}, 1./210, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {10, 0.5},
    {5, 0.1}
  },
  /* 10. PR1982, H2O, model I
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * d(H1, H2) = 1.51369612 (104.5 degree) */
  {3, {2.8, 0.4, 0.4}, {{0}},
    { {262.566, 309408} /* O-O */,
      {0, 689.348} /* O-H1 */, {0, 689.348} /* O-H2 */,
      {0, 610.455} /* H1-H1 */, {0, 610.455} /* H1-H2 */, {0, 610.455} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.513696}, 1./(KBNAC*300), LJ_FULL,
    {-0.866088, 0.433044, 0.433044}, KE2C, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {20, 0.5},
    {5, 0.5}
  },
  /* 11. PR1982, H2O, model II (SPC)
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * d(H1, H2) = 1.633081624 (109.48 degree) */
  {3, {3.166, 0.4, 0.4}, {{0}},
    { {-625.731, 629624} /* O-O */,
      {0, 225.180} /* O-H1 */, {0, 225.180} /* O-H2 */,
      {0, 0} /* H1-H1 */, {0, 0} /* H1-H2 */, {0, 0} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {1.0, 1.0, 1.633081624}, 1./(KBNAC*300), LJ_FULL,
    {-0.82, 0.41, 0.41}, KE2C, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {20, 0.5},
    {5, 0.5}
  },
  /* 12. PR1982, H2O, model III (TIPS)
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * d(H1, H2) = 1.51369612 (104.5 degree) */
  {3, {3.215, 0.4, 0.4}, {{0}},
    { {-525.000, 580000} /* O-O */,
      {0, 225.180} /* O-H1 */, {0, 225.180} /* O-H2 */,
      {0, 0} /* H1-H1 */, {0, 0} /* H1-H2 */, {0, 0} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.513696}, 1./(KBNAC*300), LJ_FULL,
    {-0.8, 0.4, 0.4}, KE2C, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {20, 0.5},
    {5, 0.5}
  },
  /* 13. SPCE, H2O
   * atom 0: O, atom 1: H1, atom 2: H2
   * the following data are copied from /Bossman/Software/3Drism/h2o_lib/spce */
  {3, {3.1666, 0.4, 0.4},
    {{78.2083543, 78.2083543}, {0, 23.150478}, {0, 23.150478}}, {{0}},
    {0.033314, 0.033314, 0.033314},
    {1.0, 1.0, 1.633}, 1./300, LJ_FULL,
    {-0.8476, 0.4238, 0.4238}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {25, 0.5},
    {5, 0.5}
  },
  /* 14. TIP3, H2O
   * atom 0: O, atom 1: H1, atom 2: H2
   * the following data are copied from /Bossman/Software/3Drism/h2o_lib/tip3 */
  {3, {3.15, 0.4, 0.4},
    { {76.5364, 76.5364}, {0, 23.1509}, {0, 23.1509} }, {{0}},
    {0.033314, 0.033314, 0.033314},
    {0.95719835, 0.95719835, 1.5139}, 1./300, LJ_FULL,
    {-0.834, 0.417, 0.417}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {25, 0.5},
    {5, 0.5}
  },
  /* 15. HRP1983, solvent + solute (+/- ion pair)
   * Cf. model 7 */
  {4, {3.341, 3.341, 3.341, 3.341},
    { {44.0, 44.0}, {44.0, 44.0}, {44.0, 44.0}, {44.0, 44.0} }, {{0}},
    {0.01867, 0.01867, 0, 0},
    {1.1}, 1./200, LJ_FULL,
    {0.200, -0.200, 1.0, -1.0}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {15, 0.5},
    {5, 0.5}
  },
  /* 16. SPCE, H2O, Na+, Cl-
   * atom 0: O, atom 1: H1, atom 2: H2, atom 3: Na, atom 4: Cl-
   * the following data are copied from /Bossman/Software/3Drism/h2o_lib/spce */
  {5, {3.1666, 0.4, 0.4, 2.35, 4.4},
    { {78.2083543, 78.2083543}, {0, 23.150478}, {0, 23.150478},
      {65.42, 65.42}, {50.32, 50.32} }, {{0}},
    {0.033314, 0.033314, 0.033314, 0, 0},
    {1.0, 1.0, 0, 0, 1.633}, 1./300, LJ_FULL,
    {-0.8476, 0.4238, 0.4238, 1, -1}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {15, 0.5},
    {5, 0.5}
  },
  {0} /* empty model, place holder */
};




#endif /* MODEL_H__ */
