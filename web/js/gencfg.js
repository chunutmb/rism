
function change_unit_eps()
{
  var unit = grab("unit_eps6").value;
  grab("unit_eps12").value = unit;
  // TODO: allow a different unit for pairs
  // and do the conversion later
  grab("pairunit_eps6").value = unit;
  grab("pairunit_eps12").value = unit;
  grab("pairunit_C6").value = unit;
  grab("pairunit_C12").value = unit;
  grab("pairunit_B").value = unit;
}



function change_eps6(id, prefix)
{
  if ( grab(prefix + "sameeps_" + id).checked ) {
    grab(prefix + "eps12_" + id).value
      = grab(prefix + "eps6_" + id).value;
  }
}



function change_ns()
{
  var ns = grab("ns").value;
  if ( !is_int(ns) ) return;
  ns = parseInt(ns);
  var tab = grab("siteTable");
  var tbody = tab.lastChild;
  var rowid, row, td;

  // determine the number of existing rows
  for ( var i = 0; ; i++ ) {
    rowid = "site_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row == null ) break;
  }
  var nrows = i; // number of existing rows

  // remove redundant elements, if any
  for ( var i = ns; i < nrows; i++ ) {
    rowid = "site_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row.parentNode == tbody )
      tbody.removeChild(row);
  }

  // create nonexisting elements
  for ( var i = nrows; i < ns; i++ ) {
    var i1 = "" + (i+1);
    rowid = "site_row_" + i1;
    row = document.createElement("tr");

    row.setAttribute("id", rowid);
    tbody.appendChild(row);

    // site id
    td = document.createElement("td");
    td.innerHTML = i1;
    row.appendChild(td);

    // sigma
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="sigma_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // eps6
    td = document.createElement("td");
    td.innerHTML = ('<input type="text" size="10" value="0" id="eps6_ROW" '
      + 'onchange="change_eps6(ROW, \'\')">').replace(/ROW/g, i1);
    row.appendChild(td);

    // eps12
    td = document.createElement("td");
    td.innerHTML = ('<input type="text" size="10" value="0" id="eps12_ROW">'
      + '<input type="checkbox" id="sameeps_ROW" onchange="change_eps6(ROW, \'\')"> '
      + 'Same as <span class="math"><i>&epsilon;</i><sub>&#8326;</sub></span>').
        replace(/ROW/g, i1);
    row.appendChild(td);

    // density, rho
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="14" value="0" id="rho_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // charge
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="charge_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);
  }
}



function change_nprs()
{
  var nprs = grab("npairs").value;
  if ( !is_int(nprs) ) return;
  nprs = parseInt(nprs);
  var tab = grab("pairTable");
  var tbody = tab.lastChild;
  var rowid, row, td;

  // determine the number of existing rows
  for ( var i = 0; ; i++ ) {
    rowid = "pair_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row == null ) break;
  }
  var nrows = i; // number of existing rows

  // remove redundant elements, if any
  for ( var i = nprs; i < nrows; i++ ) {
    rowid = "pair_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row.parentNode == tbody )
      tbody.removeChild(row);
  }

  // create nonexisting elements
  for ( var i = nrows; i < nprs; i++ ) {
    var i1 = "" + (i+1);

    rowid = "pair_row_" + i1;
    row = document.createElement("tr");

    row.setAttribute("id", rowid);
    tbody.appendChild(row);

    // site i
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="4" value="0" id="pairi_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // site j
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="4" value="0" id="pairj_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // C6
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="13" value="0" id="pairC6_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // C12
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="14" value="0" id="pairC12_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // sigma
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="4" value="0" id="pairsigma_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // eps6
    td = document.createElement("td");
    td.innerHTML = ('<input type="text" size="10" value="0" id="paireps6_ROW" '
      + 'onchange="change_eps6(ROW, \'pair\')">').
        replace(/ROW/g, i1);
    row.appendChild(td);

    // eps12
    td = document.createElement("td");
    td.innerHTML = ('<input type="text" size="10" value="0" id="paireps12_ROW"> <br>'
      + '<input type="checkbox" id="pairsameeps_ROW" onchange="change_eps6(ROW, \'pair\')"> '
      + 'Same as <span class="math"><i>&epsilon;</i><sub>&#8326;</sub></span>').
        replace(/ROW/g, i1);
    row.appendChild(td);

    // rho, in B * exp(-r/rho)
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="4" value="0" id="pairrho_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);

    // B, in B * exp(-r/rho)
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="12" value="0" id="pairB_ROW">'.
      replace(/ROW/g, i1);
    row.appendChild(td);
  }
}



function change_nbonds()
{
  var nbonds = grab("nbonds").value;
  if ( !is_int(nbonds) || nbonds < 0 ) return;
  nbonds = parseInt(nbonds);
  var tab = grab("bondTable");
  var tbody = tab.lastChild;
  var rowid, row, td;

  // determine the number of existing rows
  for ( var i = 0; ; i++ ) {
    rowid = "bond_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row == null ) break;
  }
  var nrows = i; // number of existing rows

  // remove redundant elements, if any
  for ( var i = nbonds; i < nrows; i++ ) {
    rowid = "bond_row_" + (i+1);
    row = document.getElementById(rowid);
    if ( row.parentNode == tbody )
      tbody.removeChild(row);
  }

  // create nonexisting elements
  for ( var i = nrows; i < nbonds; i++ ) {
    var i1 = "" + (i+1);
    rowid = "bond_row_" + i1;
    row = document.createElement("tr");

    row.setAttribute("id", rowid);
    tbody.appendChild(row);

    // site i
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="bondi_ROW">'.
      replace(/ROW/g, i);
    row.appendChild(td);

    // site j
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="bondj_ROW">'.
      replace(/ROW/g, i);
    row.appendChild(td);

    // bond length
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="14" value="0" id="bondlen_ROW">'.
      replace(/ROW/g, i);
    row.appendChild(td);
  }
}



function gencfg()
{
  var s = "# configuration file for rism0\n";

  var ns = grab("ns").value;
  if ( !is_int(ns) || ns <= 0 ) {
    alert("# of sites [" + grab("ns").value + "] is invalid");
    return;
  }
  s += "ns            = " + ns + "\n";

  var ljtype = grab("ljtype").value;
  for ( var i = 0; i < ns; i++ ) {
    var i1 = i + 1;
    s += "\n";
    s += "sigma(" + i1 + ")      = " + grab("sigma_"  + i1).value + "\n";
    if ( ljtype != "Hard-sphere" ) {
      var sameeps = grab("sameeps_" + i1).checked;
      if ( sameeps ) {
        s += "eps(" + i1 + ")        = " + grab("eps6_"   + i1).value + "\n";
      } else {
        s += "eps6(" + i1 + ")       = " + grab("eps6_"   + i1).value + "\n";
        s += "eps12(" + i1 + ")      = " + grab("eps12_"  + i1).value + "\n";
      }
    }
    s += "rho(" + i1 + ")        = " + grab("rho_"    + i1).value + "\n";
    s += "charge(" + i1 + ")     = " + grab("charge_" + i1).value + "\n";
  }
  s += "\n";

  var npairs = get_int("npairs");
  if ( npairs > 0 ) {
    for ( var i = 0; i < npairs; i++ ) {
      var i1 = i + 1;
      var pairi = get_int("pairi_" + i1);
      var pairj = get_int("pairj_" + i1);
      var sid = "ij(" + pairi + ", " + pairj + ")";
      var C6val = grab("pairC6_" + i1).value;
      var C12val = grab("pairC12_" + i1).value;
      var Bval = grab("pairB_" + i1).value;

      if ( parseFloat(C6val) != 0 || parseFloat(C12val) != 0 ) {
        s += "C6" + sid + "    = " + C6val + "\n";
        s += "C12" + sid + "   = " + C12val + "\n";
      } else {
        s += "sigma" + sid + " = " + grab("pairsigma_" + i1).value + "\n";
        var sameeps = grab("pairsameeps_" + i1).checked;
        if ( sameeps ) {
          s += "eps" + sid + "   = " + grab("paireps6_" + i1).value + "\n";
        } else {
          s += "eps6" + sid + "  = " + grab("paireps6_" + i1).value + "\n";
          s += "eps12" + sid + " = " + grab("paireps12_" + i1).value + "\n";
        }
      }
      if ( parseFloat(Bval) != 0 ) {
        s += "rho" + sid + "   = " + grab("pairrho_" + i1).value + "\n";
        s += "B" + sid + "     = " + grab("pairB_" + i1).value + "\n";
      }
      s += "\n";
    }
  }

  var nbonds = get_int("nbonds");
  if ( nbonds > 0 ) {
    for ( var i = 0; i < nbonds; i++ ) {
      var i1 = i + 1;
      var bondi = get_int("bondi_" + i1);
      var bondj = get_int("bondj_" + i1);
      var bondlen = get_float("bondlen_" + i1);
      s += "dis(" + bondi + ", " + bondj + ")     = " + bondlen + "\n";
    }
    s += "\n";
  }

  var unit_eps = grab("unit_eps6").value;
  var kBT, kBU, ampch;
  if ( unit_eps == "reduced" ) {
    kBT = 1;
    kBU = 1;
    ampch = 1;
  } else if ( unit_eps == "K" ) {
    kBT = 1;
    kBU = "KBNA";
    ampch = "KE2PK";
  } else if ( unit_eps == "erg" ) {
    kBT = "KB_ERG";
    kBU = 1;
    ampch = "KE2_AERG";
  } else if ( unit_eps == "kJpermol" ) {
    kBT = "KBNA";
    kBU = 1;
    ampch = "KE2NA";
  } else if ( unit_eps == "kcalpermol" ) {
    kBT = "KBNAC";
    kBU = 1;
    ampch = "KE2NAC";
  }

  s += "T             = " + grab("temp").value + "\n";
  s += "kBT           = " + kBT + "\n";
  s += "kBU           = " + kBU + "\n";
  s += "ljtype        = " + grab("ljtype").value + "\n";
  s += "ampch         = " + ampch + "\n";
  s += "rscreen       = " + grab("rscreen").value + "\n";
  s += "closure       = " + grab("ietype").value + "\n";
  s += "rmax          = " + grab("rmax").value + "\n";
  s += "npt           = " + grab("npt").value + "\n";
  s += "nlambdas      = " + grab("nlambdas").value + "\n";
  s += "itmax         = " + grab("itmax").value + "\n";
  s += "tol           = " + grab("tol").value + "\n";
  s += "solver        = " + grab("solver").value + "\n";
  s += "picard_damp   = " + grab("picard_damp").value + "\n";
  s += "lmv_damp      = " + grab("lmv_damp").value + "\n";
  s += "lmv_M         = " + grab("lmv_M").value + "\n";
  s += "mdiis_damp    = " + grab("mdiis_damp").value + "\n";
  s += "mdiis_nbases  = " + grab("mdiis_nbases").value + "\n";

  grab("cfgout").value = s;
}



function init()
{
  change_ns();
  change_nprs();
  change_nbonds();
  change_unit_eps();
}



