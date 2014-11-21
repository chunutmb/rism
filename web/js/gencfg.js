function change_ns()
{
  var ns = $("#ns").val();
  if ( !is_int(ns) ) return;
  ns = parseInt(ns);
  var tab = grab("siteTable");
  var tbody = tab.lastChild;
  
  // determine the number of existing rows
  for ( var i = 0; ; i++ ) {
    var rowid = "row_" + (i+1);
    var row = document.getElementById(rowid);
    if ( row == null ) break;
  }
  var nrows = i; // number of existing rows

  // remove redundant elements, if any
  for ( var i = ns; i < nrows; i++ ) {
    var rowid = "row_" + (i+1);
    var row = document.getElementById(rowid);
    if ( row.parentNode == tbody )
      tbody.removeChild(row);
  }

  // create nonexisting elements
  for ( var i = nrows; i < ns; i++ ) {
    var rowid = "row_" + (i+1);
    var row = document.createElement("tr");
    var td, inp;

    row.setAttribute("id", rowid);
    tbody.appendChild(row);

    // site id
    td = document.createElement("td");
    td.innerHTML = "" + (i+1);
    row.appendChild(td);
    
    // sigma
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="sigma_' + (i+1) + '">';
    row.appendChild(td);
    
    // eps6
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="eps6_' + (i+1)
      + '" onChange="change_eps6(' + (i+1) + ')">';
    row.appendChild(td);
    
    // eps12
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="eps12_' + (i+1) + '">'
      + '<input type="checkbox" id="sameeps_' + (i+1) + '"> Same as '
      + '<span class="math"><i>&epsilon;</i><sub>&#8326;</sub></span>';
    row.appendChild(td);
    
    // density, rho
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="14" value="0" id="rho_' + (i+1) + '">';
    row.appendChild(td);
    
    // charge
    td = document.createElement("td");
    td.innerHTML = '<input type="text" size="10" value="0" id="charge_' + (i+1) + '">';
    row.appendChild(td);
  }
}

function change_unit_eps()
{
  $("#unit_eps12").val( $("#unit_eps6").val() );
}

function change_eps6(id)
{
  if ( $("#sameeps_" + id).is(':checked') ) {
    $(("#eps12_" + id)).val( $(("#eps6_" + id)).val() );
  }
}

$(document).ready(function() {
  change_ns();
});

function gencfg()
{
  var s = "# configuration file for rism0\n";

  var ns = $("#ns").val();
  if ( !is_int(ns) || ns <= 0 ) {
    alert("# of sites [" + $("#ns").val() + "] is invalid");
    return;
  }
  s += "ns            = " + ns + "\n";

  var ljtype = $("#ljtype").val();
  for ( var i = 0; i < ns; i++ ) {
    var i1 = i + 1;
    s += "\n";
    s += "sigma(" + i1 + ")      = " + $(("#sigma_"  + i1)).val() + "\n";
    var sameeps = $("#sameeps_" + i1).is(':checked');
    if ( ljtype != "Hard-sphere" ) {
      if ( sameeps ) {
        s += "eps(" + i1 + ")        = " + $(("#eps6_"   + i1)).val() + "\n";
      } else {
        s += "eps6(" + i1 + ")       = " + $(("#eps6_"   + i1)).val() + "\n";
        s += "eps12(" + i1 + ")      = " + $(("#eps12_"  + i1)).val() + "\n";
      }
    }
    s += "rho(" + i1 + ")        = " + $(("#rho_"    + i1)).val() + "\n";
    s += "charge(" + i1 + ")     = " + $(("#charge_" + i1)).val() + "\n";
  }
  s += "\n";

  /* parse chemical bonds */
  cb = $("#chembonds").val().split("\n");
  for ( var i = 0; i < cb.length; i++ ) {
    var ln = cb[i].trim();
    if ( ln.charAt(0) == "#" || ln.length == 0 ) continue;
    var arr = ln.split(" ");
    if ( arr.length < 3 ) {
      console.log("invalid line " + ln);
      continue;
    }
    s += "dis(" + arr[0] + ", " + arr[1] + ")     = " + arr[2] + "\n";
  }
  s += "\n";

  var unit_eps = $("#unit_eps6").val();
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

  s += "T             = " + $("#temp").val() + "\n";
  s += "kBT           = " + kBT + "\n";
  s += "kBU           = " + kBU + "\n";
  s += "ljtype        = " + $("#ljtype").val() + "\n";
  s += "ampch         = " + ampch + "\n";
  s += "rscreen       = " + $("#rscreen").val() + "\n";
  s += "closure       = " + $("#ietype").val() + "\n";
  s += "rmax          = " + $("#rmax").val() + "\n";
  s += "npt           = " + $("#npt").val() + "\n";
  s += "nlambdas      = " + $("#nlambdas").val() + "\n";
  s += "itmax         = " + $("#itmax").val() + "\n";
  s += "tol           = " + $("#tol").val() + "\n";
  s += "solver        = " + $("#solver").val() + "\n";
  s += "picard_damp   = " + $("#picard_damp").val() + "\n";
  s += "lmv_damp      = " + $("#lmv_damp").val() + "\n";
  s += "lmv_M         = " + $("#lmv_M").val() + "\n";
  s += "mdiis_damp    = " + $("#mdiis_damp").val() + "\n";
  s += "mdiis_nbases  = " + $("#mdiis_nbases").val() + "\n";

  $("#cfgout").val( s );
}
