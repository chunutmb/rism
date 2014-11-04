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
  if ( !once ) {
    for ( i = 0; i < ns; i++ )
      fprintf(stderr, "site %d: molecule %d\n", i, mol[i]);
    once = 1;
  }

  return m->nmol = im + 1;
}



/* compute the number solvents */
static int getnsv(model_t *m)
{
  int i, c;

  getmols(m);
  for ( c = 0, i = 0; i < m->ns; i++ )
    if ( m->mol[i] == 0 ) c++;
  return c;
}



/* compute the internal energy
 * NOTE: this routine does not work for charged systems */
static int calcU(model_t *m, double **ur,
    double **cr, double **tr, double **fr,
    double *um)
{
  int i, j, ij, im, l, dohnc, ns = m->ns, npt = m->npt;
  double uij, dr, ri, gr;

  dohnc = (m->ietype == IE_HNC);
  getmols(m);
  dr = m->rmax / m->npt;
  for ( i = 0; i < m->nmol; i++ ) um[i] = 0;
  for ( i = 0; i < ns; i++ ) {
    im = m->mol[i];
    //printf("i %d, im %d\n", i, im);
    for ( j = 0; j < ns; j++ ) {
      if ( m->rho[j] <= 0 ) continue;
      ij = i * ns + j;
      /* integrating over other radius */
      uij = 0;
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        /* the following formula may produce negative values */
        gr = cr[ij][l] + tr[ij][l] + 1;
        if ( gr < 0 || dohnc ) /* only for HNC */
          gr = (fr[ij][l] + 1) * exp( tr[ij][l] );
        uij += ur[ij][l] * gr * ri * ri;
        //printf("i %d, j %d, ri %g, ur %g, g(r) %g, u(r) g(r) %g\n", i, j, ri, ur[ij][l], gr, ur[ij][l]*gr);
      }
      uij *= .5 * 4 * PI * m->rho[j] * dr;
      um[im] += uij;
      //printf("i %d, j %d, uij %g\n", i, j, uij);
    }
  }
  for ( i = 0; i < m->nmol; i++ )
    printf("mol %d: U %g\n", i, um[i]);
  return m->nmol;
}

#endif /* TDCALC_H__ */