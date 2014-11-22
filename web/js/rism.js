var ns = 3;
var ns2 = 9;
var npt = 1024;
var rmax = 20.48;
var itmax = 10000;
var tol = 1e-6;
var picard_damp = 1.0;

var beta = 1.0; // TODO
var ljtype;
var rscreen = 1.0;
var ljparams;

var ur, nrdur, fr, cr, tr;
var vrqq, vrlr, vrsr;
var fk, ck, tk;



function read_params()
{
  var i, i1;
  var sigma, eps6, eps12;

  ns = get_int("ns", 3);
  ns2 = ns * ns;
  npt = get_int("npt", 1024);
  ietype = grab("ietype").value;
  itmax = get_int("itmax", 100000);
  tol = get_float("tol", 1e-6);
  picard_damp = get_float("picard_damp", 1.0);
  rscreen = get_float("rscreen", 1.0);
  ljtype = grab("ljtype").value;
  ljparams = newnumarr(ns);
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



function initfr(lam)
{
  var i, j, ij, ji, ipr, l, use_pairpot;
  var r, z, u, uelec, nrdu, ulr, tmp;
  var sig, eps6, eps12, c6, c12, Bij, rhoij, rminij;

  for ( ipr = 0, i = 0; i < ns; i++ ) {
    for ( j = i; j < ns; j++, ipr++ ) {
      ij = i * ns + j;
      ji = j * ns + i;
      use_pairpot = false;
      sig = .5 * (ljparams[i].sigma + ljparams[j].sigma);
      eps6 = Math.sqrt( ljparams[i].eps6 * ljparams[j].eps6 );
      eps12 = Math.sqrt( ljparams[i].eps12 * ljparams[j].eps12 );
      for ( l = 0; l < npt; l++ ) {
        r = fft_ri[l];
        if ( ljtype == "Hard-sphere" ) { // TODO
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
            // TODO
          } else if ( false ) {
            // TODO
          } else {
            tmp = ljpot(r, sig, eps6, eps12, lam);
          }
          u = tmp[0];
          nrdu = tmp[1];
          uelec = ulr = 0;
          //uelec = lam * ampch * charge[i] * charge[j] / r;
          // TODO
          vrsr[ij][l] = beta * (u + uelec - ulr);
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



function prepare()
{
  var i;

  initfftw(rmax, npt);

  ur = newnumarr2d(ns2, npt);
  nrdur = newnumarr2d(ns2, npt);
  fr = newnumarr2d(ns2, npt);
  cr = newnumarr2d(ns2, npt);
  tr = newnumarr2d(ns2, npt);
  vrqq = newnumarr2d(ns2, npt);
  vrlr = newnumarr2d(ns2, npt);
  vrsr = newnumarr2d(ns2, npt);
  fk = newnumarr2d(ns2, npt);
  ck = newnumarr2d(ns2, npt);
  tk = newnumarr2d(ns2, npt);

  var lam = 1.0;
  initfr(lam);
  sphr_r2k(fr, fk);
  cparr2d(cr, fr, ns*ns, npt); // cr = fr
}



function solve()
{
  read_params();
  prepare();
  // TODO
}



function mkplot()
{
  var i, j, ij, l;

  solve();

  var options = {
    //showRoller: true,
    xlabel: '<i>r</i>',
    ylabel: '<i>c</i>(<i>r</i>)',
    yRangePad: 1,
    width: 350,
  };
  dat = "r";
  for ( i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ )
      dat += ",cr(" + i + "-" + j + ")";
  dat += "\n";
  for ( l = 0; l < npt; l++ ) {
    dat += "" + fft_ri[l];
    for ( i = 0; i < ns; i++ ) {
      for ( j = i; j < ns; j++ ) {
        ij = i * ns + j;
        dat += "," + cr[ij][l];
      }
    }
    dat += "\n";
  }
  var crplot = new Dygraph(grab("cr_plot"), dat, options);
}


