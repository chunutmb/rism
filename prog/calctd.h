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
      fprintf(stderr, "site %d: molecule %d\n", i + IDBASE, mol[i] + IDBASE);
    once = 1;
  }
  m->nmol = im + 1;

  /* compute the net charge of each molecule */
  for ( im = 0; im < m->nmol; im++ )
    m->chargemol[im] = 0;
  for ( i = 0; i < ns; i++ )
    m->chargemol[ mol[i] ] += m->charge[i];

  return m->nmol;
}



/* compute the equivalent diameter of the molecule
 * this value is only used for comparison
 * and it does not affect the solver */
static double getdiameters(model_t *m)
{
  int i, j, ipr, ns = m->ns, imol, nmol;
  double vol, si, sj, l;

  nmol = getmols(m);
  for ( imol = 0; imol < nmol; imol++ ) {
    vol = 0;
    for ( i = 0; i < ns; i++ )
      if ( m->mol[i] == imol )
        vol += pow( m->sigma[i], 3 );

    /* deduct the overlap */
    for ( ipr = 0, i = 0; i < ns; i++ )
      for ( j = i + 1; j < ns; j++, ipr++ ) {
        if ( m->mol[i] != imol || m->mol[j] != imol )
          continue;
        si = m->sigma[i];
        sj = m->sigma[j];
        l = m->dis[ipr];
        if ( l < DBL_MIN || l > (si + sj)/2 ) continue;
        vol -= (si*si*si + sj*sj*sj)/2 - (si*si + sj*sj)*l*3./4
             - pow(si*si - sj*sj, 2)*3./32/l + l*l*l/2;
      }
    m->diameter[imol] = pow(vol, 1./3);
  }
  return m->diameter[0];
}



/* compute the dielectric constant */
static double calcdielec(model_t *m)
{
  int imol, i, j, k, ic, arr[MAXATOM];
  double mu, y, eps = 0;

  getmols(m);
  for ( imol = 0; imol < m->nmol; imol++ ) {
    /* collect atoms that belong to this molecule */
    for ( ic = 0, i = 0; i < m->ns; i++ )
      if ( m->mol[i] == imol && m->rho[i] > 0 )
        arr[ic++] = i;
    /* skip an atomic solvent or solute of infinite dilution */
    if ( ic <= 1 ) continue;

    /* compute the dipole moment of molecule imol */
    if ( ic == 2 ) {
      i = arr[0];
      j = arr[1];
      mu = fabs(m->charge[j] * m->disij[i][j]);
    } else if ( ic == 3 ) {
      double rij, rjk, rki, cosang, sinang, mux, muy;
      i = arr[0];
      j = arr[1];
      k = arr[2];
      rij = m->disij[i][j];
      rjk = m->disij[j][k];
      rki = m->disij[i][k];
      cosang = (rij*rij + rki*rki - rjk*rjk)/2/(rij*rki);
      if ( cosang < -1 || cosang > 1 ) {
        fprintf(stderr, "Error: invalid geometry of "
            "molecule %d (%d,%d,%d), %g, %g, %g, cos %g\n",
            imol, i, j, k, rij, rki, rki, cosang);
        exit(1);
      }
      sinang = sqrt(1 - cosang*cosang);
      mux = m->charge[j]*rij + m->charge[j]*rki*cosang;
      muy = m->charge[j]*rki*sinang;
      mu = sqrt(mux*mux + muy*muy);
    } else {
      fprintf(stderr, "Error: dipole moment of %d atoms not implemented\n", ic);
      mu = 0;
    }

    y = 4 * PI * m->beta * m->rho[arr[0]] * mu * mu * m->ampch / 9;
    eps += y;
    //printf("imol %d, %d sites, dipole moment %g, y %g\n", imol, ic, mu, y);
  }
  eps = 1 + 3 * eps;
  return eps;
}



/* return the logarithm of the radial distribution function */
static double getlngr(model_t *m, double cr, double tr, double vrsr)
{
#ifdef GR_TINY /* user-specified threshold */
  double gr_tiny = GR_TINY;
#else
  double gr_tiny = m->tol * 100;
#endif
  double lngr = -vrsr + tr; /* the HNC case, stable */
  if ( m->ietype != IE_HNC ) { /* try to use the exact g(r) if it is positive */
    double gr0 = 1 + cr + tr;
    return ( gr0 > gr_tiny ) ? log(gr0) : lngr;
  }
  return lngr;
}



/* return the radial distribution function
 * the result is nonnegative */
static double getgr(model_t *m, double cr, double tr, double vrsr)
{
  return exp( getlngr(m, cr, tr, vrsr) );
}



/* compute the Kirkwood integrals
 * NOTE: this routine may fail for charged systems */
static int calckirk(model_t *m, double **cr, double **tr,
    double **vrsr, double *kirk)
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
        hr = getgr(m, cr[ij][l], tr[ij][l], vrsr[ij][l]) - 1;
        kij += hr * ri * ri;
      }
      kij *= 4 * PI * dr;
      if (kirk != NULL)
        kirk[j*ns + i] = kirk[i*ns + j] = kij;
      printf("i %d, j %d, Gij %g\n", i + IDBASE, j + IDBASE, kij);
    }
  }
  return 0;
}



/* compute the running coordination numbers
 * NOTE: this routine may fail for charged systems */
static int calccrdnum(model_t *m,
    double **cr, double **tr, double **vrsr,
    const char *fn)
{
  int i, j, ij, l, ns = m->ns, npt = m->npt, lmax, lmin;
  double *nij, y, dy, dr, ri, gr, gr1;
  double rimax, grmax, rimin, grmin;
  FILE *fp = NULL;

  if ( fn != NULL && (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  newarr(nij, npt + 1);
  dr = m->rmax / m->npt;
  for ( i = 0; i < ns; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      if ( m->rho[j] <= 0 ) continue;
      ij = i * ns + j;
      nij[0] = y = dy = 0;
      lmax = 0, rimax = grmax = 0;
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        gr = getgr(m, cr[ij][l], tr[ij][l], vrsr[ij][l]);
        if ( gr > grmax ) {
          lmax = l;
          rimax = ri;
          grmax = gr;
        }
        dy = m->rho[j] * 4 * PI * gr * ri * ri * dr;
        nij[l] = y + .5*dy;
        y += dy;
        if ( fp != NULL ) fprintf(fp, "%g %g %d %d\n", ri, nij[l], i, j);
      }
      if ( fp != NULL ) fprintf(fp, "\n");

      /* find the first minimum after the principle peak */
      rimin = grmin = 0, lmin = 0;
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        gr = getgr(m, cr[ij][l], tr[ij][l], vrsr[ij][l]);
        gr1 = getgr(m, cr[ij][l+1], tr[ij][l+1], vrsr[ij][l+1]);
        if ( l >= lmax && gr1 > gr ) {
          rimin = ri;
          grmin = gr;
          lmin = l;
          break;
        }
      }
      printf("i %d, j %d, princ. peak r%8.3f, g(r)%8.3f, "
          "first min.%8.3f, g(r)%8.3f, coord. no.%8.3f\n",
          i + IDBASE, j + IDBASE, rimax, grmax, rimin, grmin, nij[lmin]);
    }
  }
  delarr(nij);
  if ( fp != NULL ) fclose(fp);
  return 0;
}



/* compute the internal energy
 * NOTE: this routine may fail for charged systems */
static int calcU(model_t *m, double **ur,
    double **cr, double **tr, double **vrsr,
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
        gr = getgr(m, cr[ij][l], tr[ij][l], vrsr[ij][l]);
        uij += ur[ij][l] * gr * ri * ri;
      }
      uij *= m->kBU * .5 * 4 * PI * m->rho[j] * dr;
      um[im] += uij;
    }
  }
  for ( i = 0; i < m->nmol; i++ )
    printf("mol %d: U %g\n", i + IDBASE, um[i]);
  return m->nmol;
}



/* compute the chemical potential
 * NOTE: this routine may fail for charged systems */
static int calcchempot(model_t *m, double **cr, double **tr,
    double **vrsr, double **vrlr, double *bmum, int verbose)
{
  int i, j, ij, im, l, ns = m->ns, npt = m->npt;
  double bmuij, dr, ri, c, t, h, vrl, y, B;

  if ( m->ietype != IE_HNC && m->ietype != IE_KH )
    fprintf(stderr, "Warning: chemical potential only exact for the HNC or KH closure\n");

  getmols(m);
  dr = m->rmax / m->npt;
  for ( i = 0; i < m->nmol; i++ )
    bmum[i] = 0;
  for ( i = 0; i < ns; i++ ) {
    im = m->mol[i];
    for ( j = 0; j < ns; j++ ) {
      if ( m->rho[j] <= 0 ) continue;
      ij = i * ns + j;
      /* integrating over other radius */
      bmuij = 0;
      for ( l = 0; l < npt; l++ ) {
        ri = (l + .5) * dr;
        vrl = vrlr[ij][l];
        c = cr[ij][l] - vrl;
        t = tr[ij][l] + vrl;
        h = c + t;
        y = -c;
        if ( m->ietype == IE_KH ) {
          if ( -vrsr[ij][l] + tr[ij][l] > 0 )
            y -= 0.5 * c * h;
          else
            y += 0.5 * t * h;
        } else if ( m->ietype == IE_HNC ) {
          y += 0.5 * t * h;
        } else if ( m->ietype == IE_PY ) {
          y += 0.5 * t * h;
          if ( tr[ij][l] > -1 ) {
            B = log( 1 + tr[ij][l] ) - tr[ij][l]; /* bridge function */
            y += 2./3 * B * h + B;
          }
        }
        bmuij += y * ri * ri;
      }
      bmuij *= 4 * PI * m->rho[j] * dr;
      if ( verbose ) printf("i %d, j %d, bmuij %g\n", i + IDBASE, j + IDBASE, bmuij);
      bmum[im] += bmuij;
    }
  }
  for ( i = 0; i < m->nmol; i++ )
    printf("mol %d: chemical potential %+12g, X beta %+g\n",
        i + IDBASE, bmum[i] * m->kBU / m->beta, bmum[i]);
  return m->nmol;
}



#endif /* TDCALC_H__ */
