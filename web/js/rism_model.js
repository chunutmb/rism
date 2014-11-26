var constants = {
  "pi":               PI,
  "cal_to_j":         CAL_TO_J,
  "j_to_cal":         J_TO_CAL,
  "cal_to_kj":        CAL_TO_KJ,
  "j_to_kcal":        J_TO_KCAL,
  "kcal_to_j":        KCAL_TO_J,
  "kj_to_cal":        KJ_TO_CAL,
  "na":               NA,
  "ec":               EC,
  "eps0_si":          EPS0_SI,
  "kb_si":            KB_SI,
  "kb_kj":            KB_KJ,
  "kb_j":             KB_J,
  "kb_erg":           KB_ERG,
  "kb_kcal":          KB_KCAL,
  "kb_cal":           KB_CAL,
  "kb":               KB,
  "kbna_si":          KBNA_SI,
  "kbna_kj":          KBNA_KJ,
  "kbna_j":           KBNA_J,
  "kbna_erg":         KBNA_ERG,
  "kbna_kcal":        KBNA_KCAL,
  "kbna_cal":         KBNA_CAL,
  "kbna":             KBNA,
  "kbnac":            KBNAC,
  "ke2_si":           KE2_SI,
  "ke2_akj":          KE2_AKJ,
  "ke2_aj":           KE2_AJ,
  "ke2_aerg":         KE2_AERG,
  "ke2_akcal":        KE2_AKCAL,
  "ke2_acal":         KE2_ACAL,
  "ke2":              KE2,
  "ke2na_si":         KE2NA_SI,
  "ke2na_akj":        KE2NA_AKJ,
  "ke2na_aj":         KE2NA_AJ,
  "ke2na_aerg":       KE2NA_AERG,
  "ke2na_akcal":      KE2NA_AKCAL,
  "ke2na_acal":       KE2NA_ACAL,
  "ke2na":            KE2NA,
  "ke2nac":           KE2NAC,
  "ke2pk_si":         KE2PK_SI,
  "ke2pk_a":          KE2PK_A,
  "ke2pk":            KE2PK,
  "erg_to_j":         ERG_TO_J,
  "erg_to_cal":       ERG_TO_CAL,
  "erg_to_kj":        ERG_TO_KJ,
  "erg_to_kcal":      ERG_TO_KCAL,
  "erg_to_jpmol":     ERG_TO_JPMOL,
  "erg_to_calpmol":   ERG_TO_CALPMOL,
  "erg_to_kjpmol":    ERG_TO_KJPMOL,
  "erg_to_kcalpmol":  ERG_TO_KCALPMOL,
  "j_to_erg":         J_TO_ERG,
  "cal_to_erg":       CAL_TO_ERG,
  "kj_to_erg":        KJ_TO_ERG,
  "kcal_to_erg":      KCAL_TO_ERG,
  "jpmol_to_erg":     JPMOL_TO_ERG,
  "calpmol_to_erg":   CALPMOL_TO_ERG,
  "kjpmol_to_erg":    KJPMOL_TO_ERG,
  "kcalpmol_to_erg":  KCALPMOL_TO_ERG,
};



/* map the string value to the numerical value */
function model_map(s, table)
{
  s = s.toLowerCase();
  return ( s in table ) ? table[s] : parseFloat(s);
}



/* choose the right option for a `select' */
function model_select(id, val)
{
  var node = grab(id);
  var opts = node.options;

  for ( var i = 0; i < opts.length; i++ )
    if ( striieq(opts[i].value, val) )
      return node.selectedIndex = i;
  console.log("cannot find match:", id, opts, val);
}



/* return the index of an array index */
function model_getidx(s, n)
{
  var p = s.indexOf('('), q = s.indexOf(')'), i;

  if ( n <= 0 )
    throw new Error("getidx has array size " + n + " of " + s + ", have you set `ns'?");
  if ( p < 0 ) p = s.indexOf('[');
  if ( p < 0 )
    throw new Error("getidx cannot find the index of [" + s + "]");
  if ( q < 0 ) q = s.indexOf('[');
  if ( q < 0 )
    throw new Error("getidx cannot find the index of [" + s + "]");
  var t = s.substr(p+1, q-p-1);
  var i = parseInt(t);
  if ( isNaN(i) || i <= 0 || i > n )
    throw new Error("getidx has bad index for " + s + ", i " + i + " > " + n);
  return i;
}




/* return the index of an array index */
function model_getidx2(s, n, hasii)
{
  var p = s.indexOf('('), q = s.indexOf(')'), i;

  if ( n <= 0 )
    throw new Error("getidx2 has array size " + n + " of " + s + ", have you set `ns'?");
  if ( p < 0 ) p = s.indexOf('[');
  if ( p < 0 )
    throw new Error("getidx2 cannot find the index of [" + s + "]");
  if ( q < 0 ) q = s.indexOf('[');
  if ( q < 0 )
    throw new Error("getidx2 cannot find the index of [" + s + "]");
  var t = s.substr(p+1, q-p-1);
  if ( t.indexOf(",") < 0 )
    throw new Error(s + " has no pair indices");
  var arr = t.split(",");
  var i = parseInt( arr[0].trim() );
  var j = parseInt( arr[1].trim() );
  if ( isNaN(i) || i <= 0 || i > n
    || isNaN(j) || j <= 0 || j > n
    || (i == j && !hasii) )
    throw new Error("getidx2 has bad index for " + s + ", i " + i + ", j " + j);
  return (i < j) ? [i, j] : [j, i];
}




/* find the index in the pair table */
function model_findpair(i, j)
{
  var npairs = get_int("npairs"), tmp;

  if ( i > j ) tmp = i, i = j, j = tmp;
  // search the current pair table
  for ( var ipr = 1; ipr <= npairs; ipr++ ) {
    var pairi = get_int("pairi_" + ipr);
    var pairj = get_int("pairj_" + ipr);
    if ( pairi > pairj )
      tmp = pairi, pairi = pairj, pairj = tmp;
    if ( pairi == i && pairj == j )
      return ipr;
  }

  // the pair (i, j) does not exist
  change_npairs(npairs += 1);
  grab("pairi_" + npairs).value = i;
  grab("pairj_" + npairs).value = j;
  return npairs;
}



/* load the configuration file from the textarea */
function loadcfg(s)
{
  if ( s == null || s == undefined ) {
    s = grab("cfgout").value.trim();
  } else {
    grab("cfgout").value = s;
  }

  var i, l, ln, ij, ipr, ns = 0, nps = 0, nbonds = 0;
  var kBT = 1, kBU = 1, ampch = 1;

  // split the file into lines
  s = s.split("\n");

  change_npairs(0);
  change_nbonds(0);
  for ( l = 0; l < s.length; l++ ) {
    ln = s[l].trim();
    if ( ln.charAt(0) == "#" || ln.length < 3
      || ln.indexOf("=") < 0 )
      continue;
    arr = ln.split("=");
    key = arr[0].trim();
    val = arr[1].trim();

    //console.log(key, val);

    if ( key == "ns" ) {
      grab("ns").value = val;
      ns = parseInt(val);
      change_ns();
    } else if ( strstartswith(key, "sigma(") ) {
      i = model_getidx(key, ns);
      grab("sigma_" + i).value = parseFloat(val);
    } else if ( strstartswith(key, "eps(") ) {
      i = model_getidx(key, ns);
      grab("eps6_" + i).value = grab("eps12_" + i).value = parseFloat(val);
      grab("sameeps_" + i).checked = true;
      change_eps6(i);
    } else if ( strstartswith(key, "eps6(") ) {
      i = model_getidx(key, ns);
      grab("eps6_" + i).value = parseFloat(val);
      grab("sameeps_" + i).checked = false;
      change_eps6(i);
    } else if ( strstartswith(key, "eps12(") ) {
      i = model_getidx(key, ns);
      grab("eps12_" + i).value = parseFloat(val);
      grab("sameeps_" + i).checked = false;
      change_eps6(i);
    } else if ( strstartswith(key, "rho(") ) {
      i = model_getidx(key, ns);
      grab("rho_" + i).value = parseFloat(val);
    } else if ( strstartswith(key, "charge") ) {
      i = model_getidx(key, ns);
      grab("charge_" + i).value = parseFloat(val);
    } else if ( strstartswith(key, "dis(") ) {
      ij = model_getidx2(key, ns, false);
      grab("nbonds").value = (nbonds += 1);
      change_nbonds();
      grab("bondi_" + nbonds).value = ij[0];
      grab("bondj_" + nbonds).value = ij[1];
      grab("bondlen_" + nbonds).value = parseFloat(val);
    } else if ( strstartswith(key, "c6ij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("pairC6_" + ipr).value = val;
    } else if ( strstartswith(key, "c12ij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("pairC12_" + ipr).value = val;
    } else if ( strstartswith(key, "sigmaij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("pairsigma_" + ipr).value = val;
    } else if ( strstartswith(key, "epsij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("paireps6_" + ipr).value = val;
      grab("paireps12_" + ipr).value = val;
      grab("pairsameeps_" + ipr).checked = true;
    } else if ( strstartswith(key, "eps6ij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("paireps6_" + ipr).value = val;
      grab("pairsameeps_" + ipr).checked = false;
    } else if ( strstartswith(key, "eps12ij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("paireps12_" + ipr).value = val;
      grab("pairsameeps_" + ipr).checked = false;
    } else if ( strstartswith(key, "bij(") ) {
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("pairB_" + ipr).value = val;
    } else if ( strstartswith(key, "rhoij(") ) { // radius, not the density
      ij = model_getidx2(key, ns, true);
      ipr = model_findpair(ij[0], ij[1]);
      val = parseFloat(val);
      grab("pairrho_" + ipr).value = val;
    } else if ( strieq(key, "t")
             || strieq(key, "temp") ) {
      grab("temp").value = parseFloat(val);
    } else if ( strieq(key, "kbt") ) {
      kBT = model_map(val, constants);
    } else if ( strieq(key, "kb") ||
                strieq(key, "kbu") ||
                strieq(key, "kbe") ) {
      kBU = model_map(val, constants);
    } else if ( strieq(key, "ljtype") ) {
      model_select("ljtype", val);
    } else if ( strieq(key, "ampch") ) {
      ampch = model_map(val, constants);
    } else if ( strieq(key, "rscreen") ) {
      grab("rscreen").value = parseFloat(val);
    } else if ( strieq(key, "closure")
             || strieq(key, "ietype") ) {
      model_select("ietype", val);
    } else if ( strieq(key, "rmax") ) {
      grab("rmax").value = parseFloat(val);
    } else if ( strieq(key, "npt")
             || strieq(key, "n-pts") ) {
      grab("npt").value = parseInt(val);
    } else if ( strstartswith(key, "nlambda") ) {
      grab("nlambdas").value = parseInt(val);
    } else if ( strieq(key, "itmax") ) {
      grab("itmax").value = parseInt(val);
    } else if ( strieq(key, "tol") ) {
      grab("tol").value = parseFloat(val);
    } else if ( strieq(key, "solver") ) {
      model_select("solver", val);
    } else if ( striieq(key, "picard_damp") ) {
      grab("picard_damp").value = parseFloat(val);
    } else if ( striieq(key, "lmv_m") ) {
      grab("lmv_M").value = parseInt(val);
    } else if ( striieq(key, "lmv_damp") ) {
      grab("lmv_damp").value = parseFloat(val);
    } else if ( striieq(key, "mdiis_nbases") ) {
      grab("mdiis_nbases").value = parseInt(val);
    } else if ( striieq(key, "mdiis_damp") ) {
      grab("mdiis_damp").value = parseFloat(val);
    } else {
      console.log("Warning: unknown option", key, " = ", val);
    }
  }

  // try to figure out the unit system
  var eps = 1e-4;
  if ( Math.abs(kBT - 1) < eps
    && Math.abs(ampch - 1) < eps) {
    model_select("unit_eps6", "reduced");
  } else if ( Math.abs(kBT - 1) < eps
    && Math.abs(ampch/KE2PK - 1) < eps) {
    model_select("unit_eps6", "K");
  } else if ( Math.abs(kBT/KB_ERG - 1) < eps
    && Math.abs(ampch/KE2_AERG - 1) < eps) {
    model_select("unit_eps6", "erg");
  } else if ( Math.abs(kBT/KBNA - 1) < eps
    && Math.abs(ampch/KE2NA - 1) < eps) {
    model_select("unit_eps6", "kJpermol");
  } else if ( Math.abs(kBT/KBNAC - 1) < eps
    && Math.abs(ampch/KE2NAC - 1) < eps) {
    model_select("unit_eps6", "kcalpermol");
  } else {
    console.log("unknown unit system:", "kBT", kBT, "ampch", ampch);
  }
  change_unit_eps();
}




// use the python script data/cfg2str.py to generate the strings from configuration file
var stockmodels = {
  "empty_reduced": "ns=0\nT=1\nkBT=1\nkBU=1\nljtype=Hard-sphere\nampch=1\nrscreen=1\nclosure=PY\nrmax=20.48\nnpt=1024\nnlambdas=1\nitmax=100000\ntol=1e-6\nsolver=MDIIS\npicard_damp=1.0\nlmv_damp=0.5\nlmv_M=25\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "empty_K":  "ns=0\nT=300\nkBT=1\nkBU=KBNA\nljtype=LJ-full\nampch=KE2PK\nrscreen=1.0\nclosure=HNC\nrmax=20.48\nnpt=1024\nnlambdas=10\nitmax=100000\ntol=1e-6\nsolver=MDIIS\npicard_damp=0.01\nlmv_damp=0.5\nlmv_M=25\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "empty_erg": "ns=0\nT=300\nkBT=KB_ERG\nkBU=1\nljtype=LJ-full\nampch=KE2_AERG\nrscreen=1.0\nclosure=HNC\nrmax=20.48\nnpt=1024\nnlambdas=10\nitmax=100000\ntol=1e-6\nsolver=MDIIS\npicard_damp=0.01\nlmv_damp=0.5\nlmv_M=25\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "empty_kJpermol": "ns=0\nT=300\nkBT=KBNA\nkBU=1\nljtype=LJ-full\nampch=KE2NA\nrscreen=1.0\nclosure=HNC\nrmax=20.48\nnpt=1024\nnlambdas=10\nitmax=100000\ntol=1e-6\nsolver=MDIIS\npicard_damp=0.01\nlmv_damp=0.5\nlmv_M=25\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "empty_kcalpermol": "ns=0\nT=300\nkBT=KBNAC\nkBU=1\nljtype=LJ-full\nampch=KE2NAC\nrscreen=1.0\nclosure=HNC\nrmax=20.48\nnpt=1024\nnlambdas=10\nitmax=100000\ntol=1e-6\nsolver=MDIIS\npicard_damp=0.01\nlmv_damp=0.5\nlmv_M=25\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "lc1973_hs": "ns=2\nsigma(1)=1\neps(1)=1\nrho(1)=0.5\ncharge(1)=0\nsigma(2)=1\neps(2)=1\nrho(2)=0.5\ncharge(2)=0\ndis(1,2)=0.6\nT=1\nkBT=1\nkBU=1\nljtype=Hard-sphere\nampch=1\nrscreen=1\nclosure=PY\nrmax=5.12\nnpt=256\nnlambdas=1\ntol=1e-7\nsolver=Picard\npicard_damp=1.0\nlmv_damp=0.5\nlmv_M=15\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "chs1977_modelII": "ns=2\nsigma(1)=0.79\nrho(1)=0.686\ncharge(1)=0\nsigma(2)=1\nrho(2)=0.686\ncharge(2)=0\ndis(1,2)=0.49\nT=1\nkBT=1\nkBU=1\nljtype=Hard-sphere\nampch=1\nrscreen=1\nclosure=PY\nrmax=5.12\nnpt=256\nnlambdas=1\nitmax=1000\ntol=1e-7\nsolver=MDIIS\npicard_damp=1.0\nlmv_damp=0.5\nlmv_M=15\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "chs1977_modelIII": "ns=2\nsigma(1)=0.675\nrho(1)=0.825\ncharge(1)=0\nsigma(2)=1\nrho(2)=0.825\ncharge(2)=0\ndis(1,2)=0.346\nT=1\nkBT=1\nkBU=1\nljtype=Hard-sphere\nampch=1\nrscreen=1\nclosure=PY\nrmax=5.12\nnpt=256\nnlambdas=1\nitmax=10000\ntol=1e-7\nsolver=MDIIS\npicard_damp=0.3\nlmv_damp=0.5\nlmv_M=15\nmdiis_damp=0.5\nmdiis_nbases=5\n",
  "lc1973_n2ljrep": "ns=2\nsigma(1)=3.341\neps(1)=1\nrho(1)=0.01866\ncharge(1)=0\nsigma(2)=3.341\neps(2)=1\nrho(2)=0.01866\ncharge(2)=0\ndis(1,2)=1.100\nT=1.46\nkBT=1\nkBU=1\nljtype=LJ-repulsive\nampch=1\nrscreen=1.0\nclosure=PY\nrmax=30.72\nnpt=1024\nnlambdas=1\ntol=1e-7\nsolver=LMV\npicard_damp=0.15\nlmv_damp=0.4\nlmv_M=20\nmdiis_damp=0.5\nmdiis_nbases=3\n",
};



function loadstockmodel()
{
  var mdl = grab("stockmodel").value;

  if ( strstartswith(mdl, "model_") )
    mdl = mdl.substr(6, mdl.length - 6)
  loadcfg( stockmodels[mdl] );
}

