/* direct Picard iteration
 * do not use this unless for simple models
 * does not handle solvent-solute interactions */
function iter_picard(vrsr, wk, cr, ck, vklr, tr, tk)
{
  var it;
  var err = 0, errp = errinf;

  for ( it = 0; it < itmax; it++ ) {
    err = step_picard(null, null, vrsr, wk, cr, ck, vklr,
        tr, tk, null, 1, picard_damp);
    if ( verbose ) console.log("it", it, "err", errp, "->", err);
    if ( err < tol ) break;
    if ( err > errp ) break;
    errp = err;
  }
  return [err, it];
}



function solve()
{
  read_params();
  prepare();

  var ret, err, it, uv;
  var solver = grab("solver").value;

  for ( var ilam = 1; ilam <= nlambdas; ilam++ ) {
    var lam = 1.*ilam/nlambdas;
    initfr(lam);
    sphr_r2k(vrlr, vklr, ns, null);
    if ( ilam == 1 )
      cparr2d(cr, fr, ns2, npt); // cr = fr

    uv = uv_open(DOUU_ATOMIC);

    if ( solver == "LMV" ) {
      ret = iter_lmv(vrsr, wk, cr, der, ck, vklr,
          tr, tk, ntk, cp, uv);
    } else if ( solver == "MDIIS" ) {
      // TODO
      ret = [0, 0]
    } else {
      ret = iter_picard(vrsr, wk, cr, ck, vklr, tr, tk)
    }
    err = ret[0];
    it = ret[1];
    if ( verbose ) console.log("lambda", lam, "error", err, "it", it);
  }
}



function mkplot()
{
  var i, j, ij, l;

  solve();

  var options_gr = {
    xlabel: '<i>r</i>',
    ylabel: '<i>g</i>(<i>r</i>)',
    yRangePad: 1,
    width: 480,
    axisLabelFontSize: 10,
  };

  var options_cr = {
    xlabel: '<i>r</i>',
    ylabel: '<i>c</i>(<i>r</i>)',
    yRangePad: 1,
    width: 480,
    axisLabelFontSize: 10,
  };

  // write the header of the table
  datgr = datcr = "r";
  for ( i = 0; i < ns; i++ )
    for ( j = i; j < ns; j++ ) {
      datcr += ",cr(" + i + "-" + j + ")";
      datgr += ",gr(" + i + "-" + j + ")";
    }
  datcr += "\n";
  datgr += "\n";

  // fill  the table
  for ( l = 0; l < npt; l++ ) {
    datcr += "" + fft_ri[l];
    datgr += "" + fft_ri[l];
    for ( i = 0; i < ns; i++ ) {
      for ( j = i; j < ns; j++ ) {
        ij = i * ns + j;
        datcr += "," + cr[ij][l];
        datgr += "," + (cr[ij][l] + tr[ij][l] + 1);
      }
    }
    datcr += "\n";
    datgr += "\n";
  }
  var crplot = new Dygraph(grab("cr_plot"), datcr, options_cr);
  var grplot = new Dygraph(grab("gr_plot"), datgr, options_gr);
}


