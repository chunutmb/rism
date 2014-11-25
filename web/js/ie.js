var rho = 0.7;
var ljtype = "Hard-sphere";
var ietype = "PY";
var npt = 1024;
var rmax = 5.12;
var itmax = 10000;
var tol = 1e-6;
var picard_damp = 1.0;

var ri, ki, dr, dk;
var ur, fr, cr, tr;
var fk, ck, tk;



// read parameters
function read_params()
{
  rho = get_float("rho", 0.5);
  temp = get_float("temp", 1.0);
  ljtype = grab("ljtype").value;
  //sigma = get_float("sigma", 1.0);
  //eps = get_float("eps", 1.0)
  ietype = grab("ietype").value;
  npt = get_int("npt", 256);
  itmax = get_int("itmax", 100000);
  tol = get_float("tol", 1e-7);
  picard_damp = get_float("picard_damp", 1.0);
}



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



function initfr()
{
  var i;

  for ( i = 0; i < npt; i++ ) {
    if ( ljtype == "Hard-sphere" ) {
      if ( ri[i] < 1 ) {
        ur[i] = 10000;
        fr[i] = -1;
      } else {
        ur[i] = 0;
        fr[i] = 0;
      }
    } else {
      var invr = 1.0/ri[i], invr6;
      invr6 = invr * invr * invr;
      invr6 *= invr6;
      if ( ljtype == "LJ-full" ) {
        ur[i] = 4 * invr6 * (invr6 - 1);
      } else if ( ljtype == "LJ-repulsive" ) {
        if ( invr6 > 0.5 ) {
          ur[i] = 4 * invr6 * (invr6 - 1) + 1;
        } else {
          ur[i] = 0;
        }
      }
      fr[i] = Math.exp( -ur[i]/temp ) - 1;
    }
  }
}



function prepare()
{
  var i;

  dr = rmax / npt;
  dk = Math.PI / rmax;
  ri = newarr(npt);
  ki = newarr(npt);
  ur = newarr(npt);
  fr = newarr(npt);
  cr = newarr(npt);
  tr = newarr(npt);
  fk = newarr(npt);
  ck = newarr(npt);
  tk = newarr(npt);
  for ( i = 0; i < npt; i++ ) {
    ri[i] = (i + .5) * dr;
    ki[i] = (i + .5) * dk;
  }

  initfr();
  sphr_r2k(fr, fk);

  for ( i = 0; i < npt; i++ )
    cr[i] = fr[i];
}



function getdcr(cr, tr, fr)
{
  if ( ietype == "PY" ) {
    return fr[i] * (1 + tr[i]) - cr[i];
  } else if ( ietype == "HNC" ) {
    return (1 + fr[i]) * Math.exp( tr[i] ) - 1 - tr[i] - cr[i];
  }
}



function closure(cr, tr, fr, damp)
{
  var del, err = 0;

  for ( i = 0; i < npt; i++ ) {
    del = getdcr(cr, tr, fr);
    err = Math.max( Math.abs(del), err );
    cr[i] += del * damp;
  }
  return err;
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
    err = closure(cr, tr, fr, picard_damp);
  }
  console.log("Picard finished in " + it + " iterations, error " + err);
}



function solve()
{
  read_params();
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
    yRangePad: 1,
    width: 350,
  };
  var dat = [];
  for ( i = 0; i < npt; i++ )
    dat.push( [ ri[i], Math.max(1 + tr[i] + cr[i], 0) ] );
  var grplot = new Dygraph(document.getElementById("gr_plot"), dat, options);

  var options = {
    //showRoller: true,
    xlabel: '<i>r</i>',
    ylabel: '<i>c</i>(<i>r</i>)',
    yRangePad: 1,
    width: 350,
  };
  var dat = [];
  for ( i = 0; i < npt; i++ )
    dat.push( [ ri[i], cr[i] ] );
  var crplot = new Dygraph(grab("cr_plot"), dat, options);
}
