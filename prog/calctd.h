#ifndef TDCALC_H__
#define TDCALC_H__



/* compute thermodynamic quantities */



/* compute the distance matrix, m->disij */
static void getdisij(model_t *m)
{
  int i, j, ipr, ns = m->ns;

  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns; j++, ipr++ )
      m->disij[i][j] = m->disij[j][i] = m->dis[ipr];
}



/* compute an array that specifies the molecule id
 * from the site id
 * return the number of molecules */
static int getmols(model_t *m)
{
  int i, j, ns = m->ns;
  int q[MAXATOM+1], nq = 0, iq, im, *mol = m->mol;
  static int once;

  getdisij(m);

  for ( i = 0; i < ns; i++ ) mol[i] = -1;
  mol[0] = im = 0;
  q[nq++] = 0; /* place the first site in the queue */

  iq = 0;
  while ( 1 ) {
    i = q[iq];
    /* add neighbors of i into the queue */
    for ( j = 0; j < ns; j++ ) {
      /* if the two sites i and j are joined by a covalent bond
       * then disij[i][j] > 0, and they belong to the same molecule */
      if ( m->disij[i][j] > 0 && mol[j] < 0 ) {
        mol[j] = mol[i];
        q[nq++] = j; /* add j into the queue */
      }
    }
    if ( ++iq >= nq ) { /* queue is empty */
      /* search for unreached sites */
      for ( i = 0; i < ns; i++ )
        if ( mol[i] < 0 ) break;
      if ( i >= ns ) break; /* all sites are settled */
      mol[i] = ++im;
      q[nq++] = i;
    }
  }
  /* print out the molecular partition at the first time of calling */
  if ( !once ) {
    for ( i = 0; i < ns; i++ )
      fprintf(stderr, "site %d: molecule %d\n", i, mol[i]);
    once = 1;
  }
  return m->nmol = im + 1;
}



/* compute the Kirkword integrals */
static int calckirk(model_t *m, double **cr, double **tr,
    double *kirk)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt;
  double kij, dr, ri, hr;

  dr = m->rmax / m->npt;
  for ( i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++ ) {
      ij = i * ns + j;
      kij = 0;
      /* integrating over the radius */
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        hr = cr[ij][l] + tr[ij][l];
        kij += hr * ri * ri;
      }
      kij *= 4 * PI * dr;
      if (kirk != NULL)
        kirk[j*ns + i] = kirk[i*ns + j] = kij;
      printf("i %d, j %d, Gij %g\n", i, j, kij);
    }
  }
  return 0;
}



/* return the radial distribution function */
static double getgr(model_t *m, double cr, double tr, double fr)
{
  double gr_tiny = m->tol * 100, gr;
  gr = cr + tr + 1;
  if ( gr < gr_tiny || m->ietype == IE_HNC ) /* only for HNC */
    gr = (fr + 1) * exp( tr );
  return gr;
}



/* compute the running coordination numbers */
static int calccrdnum(model_t *m,
    double **cr, double **tr, double **fr,
    const char *fn)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt;
  double nij, dr, ri, gr;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  dr = m->rmax / m->npt;
  for ( i = 0; i < ns; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      ij = i * ns + j;
      nij = 0;
      /* integrating over other radius */
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        gr = getgr(m, cr[ij][l], tr[ij][l], fr[ij][l]);
        nij += m->rho[j] * 4 * PI * gr * ri * ri * dr;
        fprintf(fp, "%g %g %d %d\n", ri, nij, i, j);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return 0;
}



/* compute the internal energy
 * NOTE: this routine does not work for charged systems */
static int calcU(model_t *m, double **ur,
    double **cr, double **tr, double **fr,
    double *um)
{
  int i, j, ij, im, l, ns = m->ns, npt = m->npt;
  double uij, dr, ri, gr;

  getmols(m);
  dr = m->rmax / m->npt;
  for ( i = 0; i < m->nmol; i++ ) um[i] = 0;
  for ( i = 0; i < ns; i++ ) {
    im = m->mol[i];
    for ( j = 0; j < ns; j++ ) {
      if ( m->rho[j] <= 0 ) continue;
      ij = i * ns + j;
      /* integrating over other radius */
      uij = 0;
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        gr = getgr(m, cr[ij][l], tr[ij][l], fr[ij][l]);
        uij += ur[ij][l] * gr * ri * ri;
      }
      uij *= m->kB * .5 * 4 * PI * m->rho[j] * dr;
      um[im] += uij;
    }
  }
  for ( i = 0; i < m->nmol; i++ )
    printf("mol %d: U %g\n", i, um[i]);
  return m->nmol;
}



/* compute the chemical potential
 * NOTE: this routine does not work for charged systems */
static int calcmu(model_t *m, double **cr, double **tr, double *mum)
{
  int i, j, ij, im, l, ns = m->ns, npt = m->npt;
  double muij, dr, ri;

  if ( m->ietype != IE_HNC )
    fprintf(stderr, "Warning: chemical potential only works for the HNC closure\n");

  getmols(m);
  dr = m->rmax / m->npt;
  for ( i = 0; i < m->nmol; i++ ) mum[i] = 0;
  for ( i = 0; i < ns; i++ ) {
    im = m->mol[i];
    for ( j = 0; j < ns; j++ ) {
      if ( m->rho[j] <= 0 ) continue;
      ij = i * ns + j;
      /* integrating over other radius */
      muij = 0;
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        muij += (0.5*(cr[ij][l] + tr[ij][l])*tr[ij][l] - cr[ij][l]) * ri * ri;
      }
      muij *= (m->kB / m->beta) * 4 * PI * m->rho[j] * dr;
      //printf("i %d, j %d, muij %g\n", i, j, muij);
      mum[im] += muij;
    }
  }
  for ( i = 0; i < m->nmol; i++ )
    printf("mol %d: mu %g\n", i, mum[i]);
  return m->nmol;
}



#endif /* TDCALC_H__ */