var nsmax = 8;

/* check if `n' is a valid integer */
function isint(n)
{
  return !isNaN(parseInt(n)) && isFinite(n);
}

/* check if `n' is a valid floating number */
function isfloat(n)
{
  return !isNaN(parseFloat(n)) && isFinite(n);
}

function change_ns()
{
  var ns = $("#ns").val();
  if ( !isint(ns) ) return;
  if ( ns > nsmax ) ns = nsmax;
  for ( var i = 0; i < nsmax; i++ ) {
    var disabled = (i >= ns);
    $(("#sigma_"  + (i+1))).prop('disabled', disabled);
    $(("#eps6_"   + (i+1))).prop('disabled', disabled);
    $(("#eps12_"  + (i+1))).prop('disabled', disabled);
    $(("#rho_"    + (i+1))).prop('disabled', disabled);
    $(("#charge_" + (i+1))).prop('disabled', disabled);
  }
}

function change_unit_eps()
{
  $("#unit_eps12").val( $("#unit_eps6").val() );
}

$(document).ready(function() {
  $("#ns").val("0");
  change_ns();
});

function gencfg()
{
  var s = "# configuration file for rism0\n";

  var ns = $("#ns").val();
  if ( !isint(ns) || ns > nsmax || ns <= 0 ) {
    alert("# of sites [" + $("#ns").val() + "] is invalid");
    return;
  }
  s += "ns            = " + ns + "\n";

  for ( var i = 0; i < ns; i++ ) {
    var i1 = i + 1;
    s += "\n";
    s += "sigma(" + i1 + ")      = " + $(("#sigma_"  + i1)).val() + "\n";
    s += "eps6(" + i1 + ")       = " + $(("#eps6_"   + i1)).val() + "\n";
    s += "eps12(" + i1 + ")      = " + $(("#eps12_"  + i1)).val() + "\n";
    s += "rho(" + i1 + ")        = " + $(("#rho_"    + i1)).val() + "\n";
    s += "charge(" + i1 + ")     = " + $(("#charge_" + i1)).val() + "\n";
  }
  s += "\n";

  /* parse chemical bonds */
  cb = $("#chembonds").val().split("\n");
  for ( var i = 0; i < cb.length; i++ ) {
    var ln = cb[i].trim();
    if ( ln.charAt(0) == "#" ) continue;
    var arr = ln.split(" ");
    s += "dis(" + arr[0] + ", " + arr[1] + ")     = " + arr[2] + "\n";
  }
  s += "\n";

  var unit_eps = $("#unit_eps6").val();
  var kBT, kBU, ampch;
  if ( unit_eps == "K" ) {
    kBT = 1;
    kBU = "KBNA";
    ampch = "KE2PK";
  } else if ( unit_eps == "kJpermol" ) {
    kBT = "KBNA";
    kBU = 1;
    ampch = "KE2";
  } else if ( unit_eps == "kcalpermol" ) {
    kBT = "KBNAC";
    kBU = 1;
    ampch = "KE2C";
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
  s += "tol           = " + $("#tol").val() + "\n";
  s += "solver        = " + $("#solver").val() + "\n";
  s += "picard_damp   = " + $("#picard_damp").val() + "\n";
  s += "lmv_damp      = " + $("#lmv_damp").val() + "\n";
  s += "lmv_M         = " + $("#lmv_M").val() + "\n";
  s += "mdiis_damp    = " + $("#mdiis_damp").val() + "\n";
  s += "mdiis_nbases  = " + $("#mdiis_nbases").val() + "\n";

  $("#cfgout").text( s );
}
