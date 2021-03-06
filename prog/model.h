#ifndef MODEL_H__
#define MODEL_H__



/* input parameters */



/* to use more atoms, define MAXATOM when compiling, e.g.,
 *  icc -DMAXATOM=64 rism0.c -lfftw3 */
#ifndef MAXATOM
#define MAXATOM   32 /* maximal number of atoms/sites */
#endif
#define NSMAX     (MAXATOM)
#define NS2MAX    (NSMAX * NSMAX)



#define CRMAX 1e30



enum { FALSE, TRUE, BOOLEAN_COUNT };
const char *booleans[] = {"False", "True", "BOOLEAN_COUNT"};

enum { HARD_SPHERE, LJ_FULL, LJ_REPULSIVE, LJ_COUNT };
const char *ljtypes[] = {"Hard-sphere", "LJ-full", "LJ-repulsive", "LJ_COUNT"};

enum { IE_PY, IE_HNC, IE_KH, IE_COUNT };
const char *ietypes[] = {"PY", "HNC", "KH", "IE_COUNT"};

enum { SOLVER_PICARD, SOLVER_LMV, SOLVER_MDIIS, SOLVER_COUNT };
const char *solvers[] = {"Picard", "LMV", "MDIIS", "SOLVER_COUNT"};

enum { DOUU_ATOMIC, DOUU_NEVER, DOUU_ALWAYS, DOUU_ALLSOLVENT, DOUU_COUNT };
const char *douutypes[] = {"Atomic", "Never", "Always", "Allsolvent", "DOUU_COUNT"};




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
  double C6;    /* C6/r^6 */
  double C12;   /* C12/r^12 */
  double B;     /* B exp(-r/rho) */
  double rho;   /* B exp(-r/rho), this is not the density */
  double sigma; /* eps6 (sigma/r)^6 + eps12 (sigma/r)^12 */
  double eps6;  /* eps6 (sigma/r)^6 */
  double eps12; /* eps12 (sigma/r)^12 */
  double rscreen; /* screening radius */
} pairpot_t;

typedef struct {
  int ns;
  double sigma[MAXATOM];
  double eps6_12[MAXATOM][2];

  pairpot_t pairpot[MAXATOM*(MAXATOM+1)/2];
  /* alternative to sigma/epsilon, pairwise entries */

  double rho[MAXATOM];
  double dis[MAXATOM*(MAXATOM-1)/2];
  double beta; /* 1/(kBT*T), sometimes given as 1/T without kB */
  double kBT; /* Boltzmann constant, to be multiplied on T */
  double kBU; /* Boltzmann constant, only used if the unit
                 of LJ energy is already divided by kB */
  int ljtype;

  double charge[MAXATOM];
  double ampch;
  double rscreen; /* screening distance */

  int ietype;
  double rmax; /* radius cutoff */
  int npt; /* number of sampling points along r */

  int nlambdas; /* number of intermediate stages */
  int itmax; /* maximal number of iterations in each stage */
  double tol; /* tolerance of error */
  int solver; /* solver */

  picard_params_t picard;
  lmv_params_t    lmv;
  mdiis_params_t  mdiis;

  int douu; /* one of DOUU_XXX values */
  double crmax; /* maximal magnitude of c(r), can be used to
                   improve convergence for highly charged systems */

  /* the following elements are to be computed by the program */
  double disij[MAXATOM][MAXATOM]; /* matrix form of dis[] */
  int nmol, mol[MAXATOM];
  double chargemol[MAXATOM]; /* net charge of each molecule */
  double diameter[MAXATOM];
} model_t;



/* return the index of string from a predefined array */
static int model_select(const char *s, int n, const char **arr)
{
  int i;

  for ( i = 0; i < n; i++ )
    if ( strcmpfuzzy(arr[i], s) == 0 )
      return i;
  if ( isdigit(s[0]) ) {
    i = atoi(s);
    if ( i >= 0 && i < n ) return i;
  }
  fprintf(stderr, "Error: cannot select %s from the array of %d items:", s, n);
  for ( i = 0; i < n; i++ )
    fprintf(stderr, " %s", arr[i]);
  fprintf(stderr, "\n");
  exit(1);
  return 0;
}



typedef struct {
  const char *key;
  double val;
} constmap_t;



/* map the string value to the numerical value */
static double model_map(const char *s, const constmap_t *arr)
{
  int i;

  for ( i = 0; arr[i].key != NULL; i++ )
    if ( strcmpfuzzy(arr[i].key, s) == 0 )
      return arr[i].val;
  fprintf(stderr, "Error: cannot find %s in", s);
  for ( i = 0; arr[i].key != NULL; i++ )
    fprintf(stderr, " %s", arr[i].key);
  fprintf(stderr, "\n");
  exit(1);
  return 0;
}



/* return the index of an array index */
static int model_getidx(char *s, int n)
{
  char *p = strchr(s, '(');
  int i;

  if ( n <= 0 ) {
    fprintf(stderr, "Error: getidx has array size %d of %s, have you set `ns'?\n", n, s);
    exit(1);
  }
  if ( p == NULL ) p = strchr(s, '[');
  if ( p == NULL ) {
    fprintf(stderr, "Error: getidx cannot find the index of [%s]\n", s);
    exit(1);
  }
  i = atoi(p + 1) - 1;
  if ( i >= n ) {
    fprintf(stderr, "Error: getidx has bad index for %s, i %d > %d\n", s, i + 1, n);
    exit(1);
  }
  return i;
}



/* return the index of a pair
 * `hasii' is 1 if the two indices can be the same */
static int model_getidx2(char *s, int *i, int *j, int n, int hasii)
{
  char *p = strchr(s, '('), *q;
  int k;

  if ( n <= 0 ) {
    fprintf(stderr, "Error: getidx2 has array size %d for %s, have you set `ns'?\n", n, s);
    exit(1);
  }
  if ( p == NULL ) p = strchr(s, '[');
  if ( p == NULL ) {
    fprintf(stderr, "Error: getidx2 cannot find the index of [%s]\n", s);
    exit(1);
  }
  p++;
  q = strchr(p, ',');
  if ( q == NULL ) {
    fprintf(stderr, "Error: getidx2 cannot find the second index of %s\n", s);
    exit(1);
  }
  *q++ = '\0';
  *i = atoi(strstrip(p)) - 1;
  *j = atoi(strstrip(q)) - 1;
  if ( *i > *j ) k = *i, *i = *j, *j = k;
  if ( *i >= n || *j >= n ) {
    fprintf(stderr, "Warning: getidx2 bad index for %s, i %d or j %d > %d\n",
        s, *i + 1, *j + 1, n);
    exit(1);
  }
  if ( hasii ) { /* i == j is allowed */
    return n*(*i) - *i*(*i+1)/2 + *j;
  } else { /* i != j */
    return n*(*i) - (*i+1)*(*i+2)/2 + *j;
  }
}



/* load model from file `fn' */
static int model_load(model_t *m, const char *fn, int verbose)
{
  FILE *fp;
  char buf[800], *p, *key, *val;
  int i, j, ipr, ns = -1, inpar;
  double temp = 300;
  const constmap_t constants[] = {
    {"cal_to_j",        CAL_TO_J},
    {"j_to_cal",        J_TO_CAL},
    {"cal_to_kj",       CAL_TO_KJ},
    {"j_to_kcal",       J_TO_KCAL},
    {"kcal_to_j",       KCAL_TO_J},
    {"kj_to_cal",       KJ_TO_CAL},
    {"na",              NA},
    {"ec",              EC},
    {"eps0_si",         EPS0_SI},
    {"kb_si",           KB_SI},
    {"kb_kj",           KB_KJ},
    {"kb_j",            KB_J},
    {"kb_erg",          KB_ERG},
    {"kb_kcal",         KB_KCAL},
    {"kb_cal",          KB_CAL},
    {"kb",              KB},
    {"kbna_si",         KBNA_SI},
    {"kbna_kj",         KBNA_KJ},
    {"kbna_j",          KBNA_J},
    {"kbna_erg",        KBNA_ERG},
    {"kbna_kcal",       KBNA_KCAL},
    {"kbna_cal",        KBNA_CAL},
    {"kbna",            KBNA},
    {"kbnac",           KBNAC},
    {"ke2_si",          KE2_SI},
    {"ke2_akj",         KE2_AKJ},
    {"ke2_aj",          KE2_AJ},
    {"ke2_aerg",        KE2_AERG},
    {"ke2_akcal",       KE2_AKCAL},
    {"ke2_acal",        KE2_ACAL},
    {"ke2",             KE2},
    {"ke2na_si",        KE2NA_SI},
    {"ke2na_akj",       KE2NA_AKJ},
    {"ke2na_aj",        KE2NA_AJ},
    {"ke2na_aerg",      KE2NA_AERG},
    {"ke2na_akcal",     KE2NA_AKCAL},
    {"ke2na_acal",      KE2NA_ACAL},
    {"ke2na",           KE2NA},
    {"ke2nac",          KE2NAC},
    {"ke2pk_si",        KE2PK_SI},
    {"ke2pk_a",         KE2PK_A},
    {"ke2pk",           KE2PK},
    {"erg_to_j",        ERG_TO_J},
    {"erg_to_cal",      ERG_TO_CAL},
    {"erg_to_kj",       ERG_TO_KJ},
    {"erg_to_kcal",     ERG_TO_KCAL},
    {"erg_to_jpmol",    ERG_TO_JPMOL},
    {"erg_to_calpmol",  ERG_TO_CALPMOL},
    {"erg_to_kjpmol",   ERG_TO_KJPMOL},
    {"erg_to_kcalpmol", ERG_TO_KCALPMOL},
    {"j_to_erg",        J_TO_ERG},
    {"cal_to_erg",      CAL_TO_ERG},
    {"kj_to_erg",       KJ_TO_ERG},
    {"kcal_to_erg",     KCAL_TO_ERG},
    {"jpmol_to_erg",    JPMOL_TO_ERG},
    {"calpmol_to_erg",  CALPMOL_TO_ERG},
    {"kjpmol_to_erg",   KJPMOL_TO_ERG},
    {"kcalpmol_to_erg", KCALPMOL_TO_ERG},
    {NULL,              0},
  };

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

#define ECHO_(name, val, fmt) if ( verbose >= 2 ) { \
      fprintf(stderr, "%-16s= " fmt "\n", name, val); }

/* print a real number */
#define ECHO(name, val) ECHO_(name, val, "%g")

/* print an integer */
#define ECHO_INT(name, val) ECHO_(name, val, "%d")

/* print a string */
#define ECHO_STR(name, val) ECHO_(name, val, "%s")

/* print an element of a 1D array */
#define ECHO_ARR1(name, val) { char key_[32]; \
      sprintf(key_, "%s(%d)", name, i); \
      ECHO(key_, val); }

/* print an element of a 2D array */
#define ECHO_ARR2(name, val) { char key_[32], val_[32]; \
      sprintf(key_, "%s(%d, %d)", name, i, j); \
      sprintf(val_, "%-12g (pair %d)", val, ipr); \
      ECHO_STR(key_, val_); }

    if ( strcmpfuzzy(key, "ns") == 0 ) {
      m->ns = ns = atoi(val);
      if ( ns > MAXATOM ) {
        fprintf(stderr, "too many sites %d > %d, recompile the program with increased MAXATOM\n",
            ns, MAXATOM);
        exit(1);
      }
      ECHO_INT("ns", m->ns);
    } else if ( strstartswith(key, "sigma(") ) {
      i = model_getidx(key, ns);
      m->sigma[i] = atof(val);
      ECHO_ARR1("sigma", m->sigma[i]);
    } else if ( strstartswith(key, "eps(") ) {
      i = model_getidx(key, ns);
      m->eps6_12[i][0] = m->eps6_12[i][1] = atof(val);
      ECHO_ARR1("eps6/12", m->eps6_12[i][0]);
    } else if ( strstartswith(key, "eps6(") ) {
      i = model_getidx(key, ns);
      m->eps6_12[i][0] = atof(val);
      ECHO_ARR1("eps6", m->eps6_12[i][0]);
    } else if ( strstartswith(key, "eps12(") ) {
      i = model_getidx(key, ns);
      m->eps6_12[i][1] = atof(val);
      ECHO_ARR1("eps12", m->eps6_12[i][0]);
    } else if ( strstartswith(key, "rho(") ) {
      i = model_getidx(key, ns);
      m->rho[i] = atof(val);
      ECHO_ARR1("rho", m->rho[i]);
    } else if ( strstartswith(key, "charge") ) {
      i = model_getidx(key, ns);
      m->charge[i] = atof(val);
      ECHO_ARR1("charge", m->charge[i]);
    } else if ( strstartswith(key, "dis(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 0);
      m->dis[ipr] = atof(val);
      ECHO_ARR2("dis", m->dis[ipr]);
    } else if ( strstartswith(key, "c6ij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].C6 = atof(val);
      ECHO_ARR2("C6", m->pairpot[ipr].C6);
    } else if ( strstartswith(key, "c12ij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].C12 = atof(val);
      ECHO_ARR2("C12", m->pairpot[ipr].C12);
    } else if ( strstartswith(key, "sigmaij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].sigma = atof(val);
      ECHO_ARR2("sigma", m->pairpot[ipr].sigma);
    } else if ( strstartswith(key, "epsij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].eps6 = m->pairpot[ipr].eps12 = atof(val);
      ECHO_ARR2("eps6/12", m->pairpot[ipr].eps6);
    } else if ( strstartswith(key, "eps6ij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].eps6 = atof(val);
      ECHO_ARR2("eps6", m->pairpot[ipr].eps6);
    } else if ( strstartswith(key, "eps12ij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].eps12 = atof(val);
      ECHO_ARR2("eps12", m->pairpot[ipr].eps12);
    } else if ( strstartswith(key, "bij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].B = atof(val);
      ECHO_ARR2("B", m->pairpot[ipr].B);
    } else if ( strstartswith(key, "rhoij(") ) { /* radius, not the density */
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].rho = atof(val);
      ECHO_ARR2("rho", m->pairpot[ipr].rho);
    } else if ( strstartswith(key, "rscreenij(") ) {
      ipr = model_getidx2(key, &i, &j, ns, 1);
      m->pairpot[ipr].rscreen = atof(val);
      ECHO_ARR2("rscreen", m->pairpot[ipr].rscreen);
    } else if ( strcmpfuzzy(key, "T") == 0
             || strcmpfuzzy(key, "temp") == 0 ) {
      temp = atof(val);
      ECHO("T", temp);
    } else if ( strcmpfuzzy(key, "kBT") == 0 ) {
      if ( isalpha(val[0]) ) {
        m->kBT = model_map(val, constants);
      } else {
        m->kBT = atof(val);
      }
      ECHO("kBT", m->kBT);
    } else if ( strcmpfuzzy(key, "kB") == 0 ||
                strcmpfuzzy(key, "kBU") == 0 ||
                strcmpfuzzy(key, "kBE") == 0 ) {
      if ( isalpha(val[0]) ) {
        m->kBU = model_map(val, constants);
      } else {
        m->kBU = atof(val);
      }
      ECHO("kBU", m->kBU);
    } else if ( strcmpfuzzy(key, "LJ-type") == 0 ) {
      m->ljtype = model_select(val, LJ_COUNT, ljtypes);
      ECHO_STR("LJ-type", ljtypes[m->ljtype]);
    } else if ( strcmpfuzzy(key, "ampch") == 0 ) {
      if ( isalpha(val[0]) ) {
        m->ampch = model_map(val, constants);
      } else {
        m->ampch = atof(val);
      }
      ECHO("unit of e^2", m->ampch);
    } else if ( strcmpfuzzy(key, "rscreen") == 0 ) {
      m->rscreen = atof(val);
      ECHO("rscreen", m->rscreen);
    } else if ( strcmpfuzzy(key, "closure") == 0
             || strcmpfuzzy(key, "ietype") == 0 ) {
      m->ietype = model_select(val, IE_COUNT, ietypes);
      ECHO_STR("closure", ietypes[m->ietype]);
    } else if ( strcmpfuzzy(key, "rmax") == 0 ) {
      m->rmax = atof(val);
      ECHO("rmax", m->rmax);
    } else if ( strcmpfuzzy(key, "npt") == 0
             || strcmpfuzzy(key, "n-pts") == 0 ) {
      m->npt = atoi(val);
      ECHO_INT("npt", m->npt);
    } else if ( strstartswith(key, "nlambda") ) {
      m->nlambdas = atoi(val);
      ECHO_INT("nlambdas", m->nlambdas);
    } else if ( strcmpfuzzy(key, "itmax") == 0 ) {
      m->itmax = atoi(val);
      ECHO_INT("itmax", m->itmax);
    } else if ( strcmpfuzzy(key, "tol") == 0 ) {
      m->tol = atof(val);
      ECHO("tol", m->tol);
    } else if ( strcmpfuzzy(key, "solver") == 0 ) {
      m->solver = model_select(val, SOLVER_COUNT, solvers);
      ECHO_STR("solver", solvers[m->solver]);
    } else if ( strcmpfuzzy(key, "Picard-damp") == 0 ) {
      m->picard.damp = atof(val);
      ECHO("Picard-damp", m->picard.damp);
    } else if ( strcmpfuzzy(key, "LMV-M") == 0 ) {
      m->lmv.M = atoi(val);
      ECHO_INT("LMV-M", m->lmv.M);
    } else if ( strcmpfuzzy(key, "LMV-damp") == 0 ) {
      m->lmv.damp = atof(val);
      ECHO("LMV-damp", m->lmv.damp);
    } else if ( strcmpfuzzy(key, "MDIIS-nbases") == 0 ) {
      m->mdiis.nbases = atoi(val);
      ECHO_INT("MDIIS-nbases", m->mdiis.nbases);
    } else if ( strcmpfuzzy(key, "MDIIS-damp") == 0 ) {
      m->mdiis.damp = atof(val);
      ECHO("MDIIS-damp", m->mdiis.damp);
    } else if ( strcmpfuzzy(key, "douu") == 0 ) {
      m->douu = model_select(val, DOUU_COUNT, douutypes);
      ECHO_STR("douu", douutypes[m->douu]);
    } else if ( strcmpfuzzy(key, "crmax") == 0 ) {
      m->crmax = atof(val);
      ECHO("crmax", m->crmax);
    } else {
      fprintf(stderr, "Warning: unknown option %s = %s\n", key, val);
      getchar();
    }
  }
  if ( verbose >= 2 ) {
    fprintf(stderr, "press Enter to start...");
    getchar();
  }

  fclose(fp);
  m->beta = 1./(m->kBT * temp);

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
    {0}, {0}, 1./300,
    /* the unit of the Boltzmann constant is kJ/mol/K
     * however, the parameter is often already multiplied
     * in the Lennard-Jones parameters */
    1.0, KBNA, LJ_FULL,
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
    {0.500, 0.500}, {0.600},
    1.000, 1.0, 1.0, HARD_SPHERE,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_PICARD,
    {1.0},
    {0, 1.0},
    {5, 1.0}
  },
  /* 2. Model II of CHS 1977 */
  {2, {0.790, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.686, 0.686}, {0.490},
    1.000, 1.0, 1.0, HARD_SPHERE,
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
    {0.825, 0.825}, {0.346},
    1.000, 1.0, 1.0, HARD_SPHERE,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.3},
    {20, 0.5},
    {5, 1.0}
  },
  /* 4. LC1973, liquid nitrogen */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.696, 0.696}, {1.1/3.341},
    1/1.46, 1.0, 1.0, LJ_REPULSIVE,
    {0}, 0, 0,
    IE_PY, 10.24, 1024,
    1, 10000, 1e-7,
    SOLVER_MDIIS,
    {0.15},
    {20, 0.5},
    {5, 0.5}
  },
  /* 5. KA1978, liquid nitrogen */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.6964, 0.6964}, {1.1/3.341},
    1/1.61, 1.0, 1.0, LJ_FULL,
    {0}, 0, 0,
    IE_PY, 20.48, 1024,
    5, 10000, 1e-7,
    SOLVER_MDIIS,
    {0.01}, /* does not work */
    {0, 0.4},
    {5, 0.7}
  },
  /* 6. HR1981, liquid nitrogen, neutral */
  {2, {3.341, 3.341},
    {{1, 1}, {1, 1}}, {{0}},
    {0.01867, 0.01867}, {1.1},
    1/1.636, 1.0, 1.0, LJ_FULL,
    {0}, 0, 0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {10, 0.3},
    {5, 0.7}
  },
  /* 7. HR1981, liquid nitrogen, charged, also HPR1982, model I */
  {2, {3.341, 3.341},
    {{44.0, 44.0}, {44.0, 44.0}}, /* in Kelvins */ {{0}},
    {0.01867, 0.01867}, {1.1},
    1./72, 1.0, KBNA, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 1000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {10, 0.5},
    {5, 0.7}
  },
  /* 8. HPR1982, HCl, model II */
  {2, {2.735, 3.353},
    {{20.0, 20.0}, {259.0, 259.0}}, /* in Kelvins */ {{0}},
    {0.018, 0.018}, {1.257},
    1./210, 1.0, KBNA, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 1000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {10, 0.5},
    {5, 1.0}
  },
  /* 9. HPR1982, HCl, model III */
  {2, {0.4, 3.353},
    {{20.0, 20.0}, {259.0, 259.0}}, /* in Kelvins */ {{0}},
    {0.018, 0.018}, {1.3},
    1./210, 1.0, KBNA, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {10, 0.5},
    {5, 1.0}
  },
  /* 10. PR1982, H2O, model 1
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * the unit of C6 is kcal A^6 / mol
   * the unit of C12 is kcal A^12 / mol
   * d(H1, H2) = 1.51369612 (104.5 degree) */
  {3, {2.8, 0.4, 0.4}, {{0}},
    { {262.566, 309408} /* O-O */,
      {0, 689.348} /* O-H1 */, {0, 689.348} /* O-H2 */,
      {0, 610.455} /* H1-H1 */, {0, 610.455} /* H1-H2 */,
      {0, 610.455} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.513696},
    1./(KBNAC*300), KBNAC, 1.0, LJ_FULL,
    {-0.866088, 0.433044, 0.433044}, KE2NAC, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {20, 0.5},
    {5, 1.0}
  },
  /* 11. PR1982, H2O, model 2 (SPC)
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * the unit of C6 is kcal A^6 / mol
   * the unit of C12 is kcal A^12 / mol
   * d(H1, H2) = 1.633081624 (109.48 degree) */
  {3, {3.166, 0.4, 0.4}, {{0}},
    { {-625.731, 629624} /* O-O */,
      {0, 225.180} /* O-H1 */, {0, 225.180} /* O-H2 */,
      {0, 0} /* H1-H1 */, {0, 0} /* H1-H2 */,
      {0, 0} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {1.0, 1.0, 1.633081624},
    1./(KBNAC*300), KBNAC, 1.0, LJ_FULL,
    {-0.82, 0.41, 0.41}, KE2NAC, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    {0.01},
    {20, 0.5},
    {5, 1.0}
  },
  /* 12. PR1982, H2O, model 3 (TIPS)
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * the unit of C6 is kcal A^6 / mol
   * the unit of C12 is kcal A^12 / mol
   * d(H1, H2) = 1.51369612 (104.5 degree) */
  {3, {3.215, 0.4, 0.4}, {{0}},
    { {-525.000, 580000} /* O-O */,
      {0, 225.180} /* O-H1 */, {0, 225.180} /* O-H2 */,
      {0, 0} /* H1-H1 */, {0, 0} /* H1-H2 */,
      {0, 0} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.513696},
    1./(KBNAC*300), KBNAC, 1.0, LJ_FULL,
    {-0.8, 0.4, 0.4}, KE2NAC, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {20, 0.5},
    {5, 1.0}
  },
  /* 13. SPC/E, H2O
   * atom 0: O, atom 1: H1, atom 2: H2
   * the following data are copied from
   *  /Bossman/Software/3Drism/h2o_lib/spce
   * the unit of LJ energy is Kelvin
   * cf. https://en.wikipedia.org/wiki/Water_model */
  {3, {3.1666, 0.4, 0.4},
    { {78.21, 78.21}, /* O */
      {0, 23.15}, /* H1 */ {0, 23.15} /* H2 */ }, {{0}},
    {0.03334, 0.03334, 0.03334},
    {1.0, 1.0, 1.633},
    1./300, 1.0, KBNA, LJ_FULL,
    {-0.8476, 0.4238, 0.4238}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {25, 0.5},
    {5, 1.0}
  },
  /* 14. TIP3P, H2O
   * atom 0: O, atom 1: H1, atom 2: H2
   * the following data are copied from
   *  /Bossman/Software/3Drism/h2o_lib/tip3
   * the unit of LJ energy is Kelvin
   * cf. https://en.wikipedia.org/wiki/Water_model */
  {3, {3.15, 0.4, 0.4},
    { {76.5364, 76.5364}, /* O */
      {0, 23.1509}, /* H1 */ {0, 23.1509} /* H2 */ }, {{0}},
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.5139},
    1./300, 1.0, KBNA, LJ_FULL,
    {-0.834, 0.417, 0.417}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {25, 0.5},
    {5, 1.0}
  },
  /* 15. HRP1983, solvent + solute (+/- ion pair)
   * Cf. model 7 */
  {4, {3.341, 3.341, 3.341, 3.341},
    { {44.0, 44.0}, {44.0, 44.0},
      {44.0, 44.0}, {44.0, 44.0} }, /* in Kelvins */
    {{0}},
    {0.01867, 0.01867, 0, 0},
    {1.1},
    1./200, 1.0, KBNA, LJ_FULL,
    {0.200, -0.200, 1.0, -1.0}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {15, 0.5},
    {5, 1.0}
  },
  /* 16. SPC/E, H2O, Na+, Cl-
   * atom 0: O, atom 1: H1, atom 2: H2, atom 3: Na, atom 4: Cl-
   * the following data are copied from
   *  /Bossman/Software/3Drism/h2o_lib/spce
   * the unit of LJ energy is Kelvin */
  {5, {3.1666, 0.4, 0.4, 2.35, 4.4},
    { {78.2083543, 78.2083543}, /* O */
      {0, 23.150478}, /* H1 */ {0, 23.150478}, /* H2 */
      {65.42, 65.42}, {50.32, 50.32} }, {{0}},
    {0.03334, 0.03334, 0.03334, 0.0, 0.0},
    {1.0, 1.0, 0, 0, 1.633},
    1./300, 1.0, KBNA, LJ_FULL,
    {-0.8476, 0.4238, 0.4238, 1, -1}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_MDIIS,
    {0.01},
    {25, 0.5},
    {5, 1.0}
  },
  {0} /* empty model, place holder */
};



/* model that contains user settings */
model_t model_usr[1] = { {0, {0.0}, {{0.0}}, {{0}}, /* pairpot */
  {0.0}, {0.0}, 0.0, 0.0, 0.0, -1, /* ljtype */
  {0.0}, 0.0, 0.0, -1, /* ietype */
  0.0, 0, 0, 0, 0.0, -1, /* solver */
  {0}, {0}, {0},
  -1, /* douu */
  0.0, /* crmax */
  {{0.0}}, 0, {0}, {0.0}, {0.0}
} };



/* replace a zero by DBL_MIN */
static double replzero(double x)
{
  return fabs(x) <= 0 ? DBL_MIN : x;
}



/* override an array element
 * the input string assumes the format of `i,x' */
static int model_register_arr(double *arr, char *s)
{
  char *p;
  int i;

  /* split the string at the two commas */
  if ( (p = strchr(s, ',')) == NULL ) {
    fprintf(stderr, "invalid format for array: %s\n", s);
    return -1;
  }
  *p = '\0';
  i = atoi(s);
  if ( i <= 0 || i > NSMAX ) {
    fprintf(stderr, "invalid index %d\n", i);
    return -1;
  }
  arr[i-1] = replzero( atof(p+1) );
  return 0;
}




/* override the distance of a covalent bond
 * the input string assumes the format of `i,j,dis'
 * the distance is saved only to disij, not dis
 * for the number of sites is not yet known */
static int model_register_disij(model_t *m_usr, char *s)
{
  char *p1, *p2;
  int i, j;
  double dis;

  /* split the string at the two commas */
  if ( (p1 = strchr(s, ',')) == NULL
    || (p2 = strchr(p1 + 1, ',')) == NULL ) {
    fprintf(stderr, "invalid format for the distance %s\n", s);
    return -1;
  }
  *p1 = '\0';
  *p2 = '\0';
  i = atoi(s);
  j = atoi(p1 + 1);
  if ( i <= 0 || j <= 0 || i > NSMAX || j > NSMAX || i == j ) {
    fprintf(stderr, "invalid distance pair %d and %d\n", i, j);
    return -1;
  }
  i -= 1;
  j -= 1;
  p2++;
  if ( strncmpfuzzy(p2, "inf", 3) == 0 )
    dis = -1; /* this means infinity; note, 0 will be ignored */
  else
    dis = replzero( atof(p2) );
  m_usr->disij[i][j] = m_usr->disij[j][i] = dis;
  return 0;
}



/* override settings in `m' according to `m_usr' */
static int model_override(model_t *m, const model_t *m_usr)
{
  int i, j, ipr, ns = m->ns;
  double x;

  /* override the density and charges */
  for ( i = 0; i < ns; i++ ) {
    if ( (x = m_usr->rho[i]) > 0 ) {
      m->rho[i] = m_usr->rho[i];
      fprintf(stderr, "density of %d is set to %g\n", i+IDBASE, m->rho[i]);
    }
    if ( fabs(x = m_usr->charge[i]) > 0 ) {
      m->charge[i] = m_usr->charge[i];
      fprintf(stderr, "charge of %d is set to %g\n", i+IDBASE, m->charge[i]);
    }
  }

  /* override the distances */
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns; j++, ipr++ )
      if ( (x = m_usr->disij[i][j]) > 0 || x < 0 ) {
        if ( x < 0 ) x = 0;
        m->dis[ipr] = m->disij[i][j] = m->disij[j][i] = x;
        fprintf(stderr, "distance between %d and %d is set to %g\n",
            i + IDBASE, j + IDBASE, x);
      }

  /* override the temperature */
  if ( m_usr->beta > 0 )
    m->beta = m_usr->beta / m->kBT;

  /* override the solver of the integral equation */
  if ( m_usr->solver >= 0  && m_usr->solver < SOLVER_COUNT ) {
    m->solver = m_usr->solver;
    fprintf(stderr, "solver is changed to %s\n", solvers[m->solver]);
  }

  /* override the closure */
  if ( m_usr->ietype >= 0 && m_usr->ietype < IE_COUNT ) {
    m->ietype = m_usr->ietype;
    fprintf(stderr, "ietype is changed to %s\n", ietypes[m->ietype]);
  }

  /* override how to do solute-solute interactions */
  if ( m_usr->douu >= 0 && m_usr->douu < DOUU_COUNT
    && m_usr->douu != m->douu ) {
    m->douu = m_usr->douu;
    fprintf(stderr, "douu is changed to %s\n", douutypes[m->douu]);
  }

  /* override the number of lambdas */
  if ( m_usr->nlambdas > 0 )
    m->nlambdas = m_usr->nlambdas;

  /* override the maximal value of c(r) */
  m->crmax = m_usr->crmax;
  if ( m->crmax <= 0 ) m->crmax = CRMAX;
  return 0;
}



#endif /* MODEL_H__ */
