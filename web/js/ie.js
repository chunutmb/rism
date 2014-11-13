var rho = 0.7;
var npt = 256;
var rmax = 5.12;
var itmax = 10000;
var tol = 1e-6;
var damp = 0.5;

var ri, ki, dr, dk;
var fr, cr, tr;
var fk, ck, tk;

function sphr_r2k(ar, ak)
{
  var i;

  for ( i = 0; i < npt; i++ ) ak[i] = ar[i];
  fft3dsphr11(ak, npt, dr, dk, 1);
}


function sphr_k2r(ak, ar)
{
  var i;

  for ( i = 0; i < npt; i++ ) ar[i] = ak[i];
  fft3dsphr11(ar, npt, dk, dr, 0);
}


function prepare()
{
  var i;

  dr = rmax / npt;
  dk = Math.PI / rmax;
  ri = new Array(npt);
  ki = new Array(npt);
  fr = new Array(npt);
  cr = new Array(npt);
  tr = new Array(npt);
  fk = new Array(npt);
  ck = new Array(npt);
  tk = new Array(npt);
  for ( i = 0; i < npt; i++ ) {
    ri[i] = (i + .5) * dr;
    ki[i] = (i + .5) * dk;
  }

  for ( i = 0; i < npt; i++ )
    fr[i] = (ri[i] < 1) ? -1 : 0;
  sphr_r2k(fr, fk);

  for ( i = 0; i < npt; i++ )
    cr[i] = fr[i];
}

function picard()
{
  var it;
  var err = 1e9;

  for ( it = 1; it <= itmax && err > tol; it++ ) {
    sphr_r2k(cr, ck);
    // Ornstein-Zernike relation
    for ( i = 0; i < npt; i++ )
      tk[i] = rho * ck[i] * ck[i]/(1 - rho * ck[i]);
    sphr_k2r(tk, tr);
    // closure
    var del;
    err = 0;
    for ( i = 0; i < npt; i++ ) {
      del = fr[i] * (1 + tr[i]) - cr[i];
      err = Math.max( Math.abs(del), err );
      cr[i] += del * damp;
    }
  }
  console.log("Picard finished in " + it + " iterations, error " + err);
}

function solve()
{
  prepare();
  picard();
}

function mkplot()
{
  var i;

  solve();

  var options = {
    //showRoller: true,
    xlabel: '<i>r</i>',
    ylabel: '<i>g</i>(<i>r</i>)',
    yRangePad: 0,
  };
  var dat = "r,gr\n";
  for ( i = 0; i < npt; i++ ) {
    var gr = 1 + tr[i] + cr[i];
    dat += "" + ri[i] + "," + gr + "\n";
  }
  var grplot = new Dygraph(document.getElementById("gr_plot"), dat, options);

  var options = {
    //showRoller: true,
    xlabel: '<i>r</i>',
    ylabel: '<i>c</i>(<i>r</i>)',
    yRangePad: 0,
  };
  var dat = "r,cr\n";
  for ( i = 0; i < npt; i++ ) {
    dat += "" + ri[i] + "," + cr[i] + "\n";
  }
  var crplot = new Dygraph(document.getElementById("cr_plot"), dat, options);
}

$(document).ready(new function() {
});
