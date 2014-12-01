var ns = 3;
var ns2 = 9;
var npt = 1024;
var rmax = 20.48;
var itmax = 10000;
var tol = 1e-7;
var picard_damp = 1.0;

var temp = 1.0;
var beta = 1.0; // to be computed
var kBT = 1, kBU = 1, ampch = 1;

var ljtype;
var rscreen = 1.0;
var ljparams;
var pairparams;
var charge;
var rho; // density
var dis; // chemical bonds

var nlambdas = 1;

var crmax = 1e30;

var verbose = true;

var ur, nrdur, fr, cr, tr, cp;
var vrqq, vrlr, vrsr, vklr;
var ck, tk, wk, ntk, der;



function change_verbose()
{
  verbose = grab("verbose").checked;
}



function read_params()
{
  var i, i1, j;
  var sigma, eps6, eps12;

  ns = get_int("ns", 3);
  ns2 = ns * ns;
  npt = get_int("npt", 1024);
  rmax = get_float("rmax", 20.48);
  ietype = grab("ietype").value;
  itmax = get_int("itmax", 100000);
  tol = get_float("tol", 1e-7);
  picard_damp = get_float("picard_damp", 1.0);
  change_verbose();

  temp = get_float("temp", 1.0);
  var unit_eps = grab("unit_eps6").value;
  if ( unit_eps == "reduced" ) {
    kBT = 1;
    kBU = 1;
    ampch = 1;
  } else if ( unit_eps == "K" ) {
    kBT = 1;
    kBU = 0.00831446214547;
    ampch = 167100.9566;
  } else if ( unit_eps == "erg" ) {
    kBT = 1.3806488e-16;
    kBU = 1;
    ampch = 2.3070773523707e-11;
  } else if ( unit_eps == "kJpermol" ) {
    kBT = 0.00831446214547;
    kBU = 1;
    ampch = 1389.354578;
  } else if ( unit_eps == "kcalpermol" ) {
    kBT = 0.00198720414667;
    kBU = 1;
    ampch = 322.0637137;
  }
  beta = 1/(kBT*temp);

  rscreen = get_float("rscreen", 1.0);

  nlambdas = Math.max( get_int("nlambdas", 1), 1 );

  crmax = get_float("crmax", 1e30);

  // Lennard-Jones parameters
  ljtype = grab("ljtype").value;
  ljparams = newarr(ns);
  charge = newarr(ns);
  rho = newarr(ns);
  for ( i = 0; i < ns; i++ ) {
    i1 = i + 1;
    sigma = get_float("sigma_" + i1, 1.0);
    eps6 = get_float("eps6_" + i1, 1.0);
    if ( grab("sameeps_" + i1).checked ) {
      eps12 = eps6;
    } else {
      eps12 = get_float("eps12_" + i1, 1.0);
    }
    ljparams[i] = { "sigma": sigma, "eps6": eps6, "eps12": eps12 };
    charge[i] = get_float("charge_" + i1, 0.0);
    rho[i] = get_float("rho_" + i1, 0.0);
  }

  // pair parameters
  pairparams = newarr(ns*ns);
  var npr = get_int("npairs", 0), ipr;
  for ( ipr = 0; ipr < npr; ipr++ ) {
    i1 = ipr + 1;
    var pairi = get_int("pairi_" + i1, 0) - 1;
    var pairj = get_int("pairj_" + i1, 0) - 1;
    if ( pairi < 0 || pairj < 0 ) {
      console.log("invalid pair between " + pairi + " and " + pairj);
      continue;
    }
    var C6ij = get_float("pairC6_" + i1, 0);
    var C12ij = get_float("pairC12_" + i1, 0);
    var sigmaij = get_float("pairsigma_" + i1, 0);
    var eps6ij = get_float("paireps6_" + i1, 0);
    var eps12ij = get_float("paireps12_" + i1, 0);
    var Bij = get_float("pairB_" + i1, 0);
    var rhoij = get_float("pairrho_" + i1, 0);
    pairparams[pairi*ns + pairj] = pairparams[pairj*ns + pairi] = {
      C6: C6ij,
      C12: C12ij,
      sigma: sigmaij,
      eps6: eps6ij,
      eps12: eps12ij,
      B: Bij,
      rho: rhoij,
    }
  }

  // distance matrix
  dis = newarr(ns*ns);
  var nbonds = get_int("nbonds", 0), ibond;
  for ( ibond = 0; ibond < nbonds; ibond++ ) {
    i1 = ibond + 1;
    var bondi = get_int("bondi_" + i1, 0) - 1;
    var bondj = get_int("bondj_" + i1, 0) - 1;
    if ( bondi < 0 || bondj < 0 || bondi == bondj ) {
      console.log("invalid pair between " + bondi + " and " + bondj);
      continue;
    }
    var dij = get_float("bondlen_" + i1, 0);
    dis[bondi*ns + bondj] = dis[bondj*ns + bondi] = dij;
  }
}



/* Lennard-Jones potential in terms of sigma and epsilon */
function ljpot(r, sig, eps6, eps12, lam)
{
  var ir, ir6, u1, u2, eps, nrdu;

  ir = sig/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  if ( eps6 > 0 ) { // LJ with an attractive tail
    if ( ir6 < eps6/(2*eps12) ) { // r > rm
      u1 = 0; // repulsive
      u2 = 4*ir6*(ir6*eps12 - eps6); // attractive
      nrdu = ir6 * (48 * ir6 * eps12 - 24 * eps6) * lam;
    } else { // r < rm
      eps = eps6 * eps6 / eps12;
      u1 = 4*ir6*(ir6*eps12 - eps6) + eps; // repulsive
      u2 = -eps; // attractive
      nrdu = ir6 * (48 * ir6 * eps12 - 24 * eps6);
    }
  } else { // purely repulsive
    u1 = 4*ir6*(ir6*eps12 - eps6); // repulsive
    u2 = 0; // no attractive tail
    nrdu = ir6 * (48 * ir6 * eps12 - 24 * eps6);
  }
  return [u1 + u2 * lam, nrdu];
}



/* Lennard-Jones potential in terms of C6 and C12 */
function ljpot6_12(r, c6, c12, lam)
{
  var sig, eps, ir, ir6;

  if ( c6 < 0 && c12 > 0 ) { /* attractive */
    sig = Math.pow(-c12/c6, 1./6);
    eps = c6*c6/4/c12;
    return ljpot(r, sig, eps, eps, lam);
  }
  if ( c6 > 0 ) lam = 1;
  ir = 1/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  return [ir6 * (c12 * ir6 + c6 * lam),
          ir6 * (12 * c12 * ir6 + 6 * c6 * lam)];
}



/* repulsive Lennard-Jones potential */
function ljrpot(r, sig, eps6, eps12)
{
  var ir, ir6, eps = 0, nrdu;

  ir = sig/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;
  if ( ir6 < eps6/(2*eps12) ) return [0, 0];
  if ( eps6 > 0 && eps12 > 0 )
    eps = eps6 * eps6 / eps12;
  return [4*ir6*(ir6*eps12 - eps6) + eps,
          ir6*(48*ir6*eps12 - 24*eps6)];
}



/* Huggins-Mayer potential: B * exp(-r/rho) - C/r^6
 * which is used in
 * Alkali halides in water: Ion-solvent correlations and ion-ion potentials of
 * mean force at infinite dilution
 * B. Montgomery Pettitt and Peter J. Rossky
 * J. Chem. Phys. 84(10) 5836-5844 (1986) */
function HMpot(r, B, C, rho, lam)
{
  var ir, ir6, r1;

  ir = 1/r;
  ir6 = ir * ir * ir;
  ir6 *= ir6;

  r1 = r/rho;
  B *= exp(-r1);

  return [B - C*ir6,
          B * r1 - 6 * C * ir6 * lam];
}



/* find the peak radius of the short-range part of the Huggins-Mayer potential:
 *  B exp(-r/rho) - C/r^6 */
function solveHMrmin(B, C, rho)
{
  var r, rp;

  if ( B <= 0 || C <= 0 ) return 0;
  rp = rho;
  for ( var i = 0; i < 1000; i++ ) {
    r = pow(6*C*rho/B*exp(rp/rho), 1./7);
    if (fabs(r - rp) <  1e-10) break;
    rp = r;
  }
  return r;
}



function initfr(lam)
{
  var vrmin = -30;
  var i, j, ij, ji, ipr, l, use_pairpot;
  var r, z, u, uelec, nrdu, ulr, rscrn, u_du;
  var sig, eps6, eps12, c6, c12, Bij, rhoij, rminij;

  for ( ipr = 0, i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++, ipr++ ) {
      ij = i * ns + j;
      ji = j * ns + i;
      use_pairpot = (pairparams[ij] != 0);
      if ( !use_pairpot ) {
        sig = .5 * (ljparams[i].sigma + ljparams[j].sigma);
        eps6 = Math.sqrt( ljparams[i].eps6 * ljparams[j].eps6 );
        eps12 = Math.sqrt( ljparams[i].eps12 * ljparams[j].eps12 );
      } else {
        sig = pairparams[ij].sigma;
        eps6 = pairparams[ij].eps6;
        eps12 = pairparams[ij].eps12;
        c6 = pairparams[ij].C6;
        c12 = pairparams[ij].C12;
        Bij = pairparams[ij].B;
        rhoij = pairparams[ij].rho;
        rminij = solveHMrmin(Bij, -c6, rhoij);
      }
      rscrn = rscreen * Math.sqrt(2);
      for ( l = 0; l < npt; l++ ) {
        r = fft_ri[l];
        if ( ljtype == "Hard-sphere" ) {
          if ( r < sig ) {
            vrsr[ij][l] = 2*INFTY;
            z = -1;
          } else {
            vrsr[ij][l] = 0;
            z = 0;
          }
          nrdur[ij][l] = ur[ij][l] = 0;
          vrlr[ij][l] = vrqq[ij][l] = 0;
        } else { // Lennard-Jones
          if ( use_pairpot ) {
            if ( Bij > 0 ) {
              u_du = HMpot((r < rminij ? rmin : r), Bij, -c6, rhoij, lam)
            } else if ( Math.abs(c6)  > 0 || Math.abs(c12) > 0 ) {
              u_du = ljpot6_12(r, c6, c12, lam);
            } else {
              u_du = ljpot(r, sig, eps6, eps12, lam);
            }
          } else if ( ljtype == "LJ-repulsive" ) {
            u_du = ljrpot(r, sig, eps6, eps12);
          } else {
            u_du = ljpot(r, sig, eps6, eps12, lam);
          }
          u = u_du[0];
          nrdu = u_du[1];
          uelec = ulr = 0;
          uelec = lam * ampch * charge[i] * charge[j] / r;
          ur[ij][l] = u + uelec;
          nrdur[ij][l] = nrdu + uelec;
          ulr = uelec * erf( r / rscrn );
          vrqq[ij][l] = beta * uelec;
          vrlr[ij][l] = beta * ulr;
          vrsr[ij][l] = beta * (u + uelec - ulr);
          if ( vrsr[ij][l] < vrmin ) {
            console.log("vr truncated: ", i, j, r, vrsr[ij][l]);
            vrsr[ij][l] = vrmin;
          }
          z = Math.exp(-vrsr[ij][l]) - 1;
        }
        fr[ij][l] = z;
      }
      if ( j > i ) {
        cparr(ur[ji],     ur[ij],     npt);
        cparr(nrdur[ji],  nrdur[ij],  npt);
        cparr(vrqq[ji],   vrqq[ij],   npt);
        cparr(vrlr[ji],   vrlr[ij],   npt);
        cparr(vrsr[ji],   vrsr[ij],   npt);
        cparr(fr[ji],     fr[ij],     npt);
      }
    }
  }
}



/* initialize the w matrix for intra-molecular covalent bonds */
function initwk(wk)
{
  var i, j, ij, u, l;

  for ( i = 0; i < ns; i++ ) {
    for ( u = 0; u < npt; u++ ) // diagonal
      wk[i*ns + i][u] = 1.0;
    for ( j = i + 1; j < ns; j++ ) {
      ij = i*ns + j;
      l = dis[ij];
      if ( l <= 0 ) continue;
      for ( u = 0; u < npt; u++ ) {
        var k = fft_ki[u];
        wk[ij][u] = Math.sin(k*l)/(k*l);
      }
      cparr(wk[j*ns + i], wk[ij], npt);
    }
  }
}



function prepare()
{
  var i;

  ur = newarr2d(ns2, npt);
  nrdur = newarr2d(ns2, npt);
  vrqq = newarr2d(ns2, npt);
  vrlr = newarr2d(ns2, npt);
  vrsr = newarr2d(ns2, npt);
  vklr = newarr2d(ns2, npt);
  fr = newarr2d(ns2, npt);
  wk = newarr2d(ns2, npt);
  cr = newarr2d(ns2, npt);
  ck = newarr2d(ns2, npt);
  cp = newarr2d(ns2, npt);
  tr = newarr2d(ns2, npt);
  tk = newarr2d(ns2, npt);
  ntk = newarr2d(ns2, npt);
  der = newarr2d(ns2, npt);

  initfftw(rmax, npt);
  initwk(wk);
}



/* Ornstein-Zernike relation: c(k) --> t(k) */
function oz(ck, vklr, tk, wk, invwc1w)
{
  var i, j, ij, l;
  var w = newarr(ns2);
  var dw = newarr(ns2);
  var c = newarr(ns2);
  var wc = newarr(ns2);
  var invwc1 = newarr(ns2);
  var tm1 = newarr(ns2);
  var tm2 = newarr(ns2);
  var tm3 = newarr(ns2);

  for ( l = 0; l < npt; l++ ) {
    for ( ij = 0; ij < ns * ns; ij++ ) {
      w[ij] = wk[ij][l];
      c[ij] = ck[ij][l] - vklr[ij][l];
    }

    /* note that w c is not symmetric w.r.t. i and j */
    matmul(wc, w, c, ns); /* wc = w c */
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ ) {
        ij = i*ns + j;
        tm1[ij] = rho[i] * wc[ij]; /* tm1 = rho w c */
        tm2[ij] = (i == j) - tm1[ij]; /* tm2 = 1 - rho w c */
        dw[ij] = w[ij] - (i == j);
      }

    invmat(tm2, invwc1, ns); /* invwc1 = (1 - wc)^(-1) */

    matmul(tm3, invwc1, w, ns); /* tm3 = (1 - rho w c)^(-1) w */
    if ( invwc1w != null )
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
function getcr(tr, vrsr, ietype, crmax)
{
  var xp, del, fr, cr, dcr;

  del = -vrsr + tr;
  if ( ietype == "HNC" ) {
    xp = Math.exp(del);
    cr = xp - 1 - tr;
    dcr = xp - 1;
  } else if ( ietype == "PY" ) {
    fr = Math.exp(-vrsr) - 1;
    cr = fr * (1 + tr);
    dcr = fr;
  } else if ( ietype == "KH" ) {
    if ( del <= 0 ) { // HNC
      xp = Math.exp(del);
      cr = xp - 1 - tr;
      dcr = xp - 1;
    } else {
      cr = -vrsr;
      dcr = 0;
    }
  } else {
    throw new Error("unknown closure " + ietype);
  }
  if ( cr > crmax ) cr = crmax;
  else if ( cr < -crmax ) cr = -crmax;
  return [cr, dcr];
}



/* apply the closure
 * compute residue vector if needed */
function closure(res, der, vrsr, cr, tr, prmask, update, damp)
{
  var i, j, ij, ji, id, l;
  var ret, y, err, max, errm = 0;

  for ( id = 0, i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      ij = i*ns + j;
      ji = j*ns + i;
      if ( prmask != null && !prmask[ij] ) continue;
      err = max = 0;
      for ( l = 0; l < npt; l++, id++ ) {
        ret = getcr(tr[ij][l], vrsr[ij][l], ietype, crmax);
        y = ret[0] - cr[ij][l];
        if ( der != null ) der[ij][l] = ret[1];
        if ( res != null ) res[id] = y;
        if ( update )
          cr[ij][l] += damp * y;
        err = Math.max(err, Math.abs(y));
        max = Math.max(max, Math.abs(cr[ij][l]));
      }
      // the c(r) between two ions can be extremely large
      // so we use the relative error to be compared with the tolerance
      if ( (err /= (max + 1e-6)) > errm ) errm = err;
      if ( j > i ) cparr(cr[ji], cr[ij], npt);
    }
  return errm;
}



/* a step of direct iteration (Picard)
 * compute residue vector if needed */
function step_picard(res, der, vrsr, wk, cr, ck, vklr,
    tr, tk, prmask, update, damp)
{
  sphr_r2k(cr, ck, ns, null);
  oz(ck, vklr, tk, wk, null);
  sphr_k2r(tk, tr, ns, null);
  return closure(res, der, vrsr, cr, tr, prmask, update, damp);
}





