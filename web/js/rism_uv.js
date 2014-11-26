var DOUU_NEVER = 0;
var DOUU_ALWAYS = 1;
var DOUU_ATOMIC = 2;

var SOLVENT_SOLVENT = 0;
var SOLUTE_SOLVENT = 1;
var SOLUTE_SOLUTE = 2;
var LAST_STAGE = 3;



/* compute the number solvents
 * here, solvent are defined as molecules with nonzero density
 * zero-density molecules must appear at the end of model */
function getnsv()
{
  for ( var i = 0; i < ns; i++ )
    if ( rho[i] <= 0 ) return i;
  return ns;
}



function uv_open(uutype)
{
  var i, j, ipr;

  var uv = {
    prmask: newarr(ns*ns),
    ns: ns,
    nsv: getnsv(), // number of solvent sites
    npr: 0, // number of active pairs
    stage: SOLVENT_SOLVENT,
    infdil: false, // infinite dilution
    atomicsolute: false, // solvent
    uutype: 0,
  };

  for ( i = 0; i < ns; i++ )
    for ( j = 0; j < ns; j++ )
      uv.prmask[i*ns + j] = (i < uv.nsv && j < uv.nsv);
  uv.npr = uv.nsv * (uv.nsv + 1) / 2;
  for ( i = uv.nsv; i < uv.ns; i++ )
    if ( rho[i] > 0 )
      break;
  uv.infdil = (i == uv.ns);
  // see if the solute is atomic
  uv.atomicsolute = true;
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns; j++, ipr++ )
      if ( i >= uv.nsv && j >= uv.nsv && dis[i*ns + j] > 0 ) {
        uv.atomicsolute = false;
        break;
      }
  uv.uutype = uutype;
  return uv;
}



/* switch between stages
 * solvent-solvent --> solvent-solute --> solute-solute
 * return 1 if the iteration should stop */
function uv_switch(uv)
{
  var i, j, nsv = uv.nsv;

  /* switch between stages */
  if ( uv.stage == SOLVENT_SOLVENT ) { // turn on solute-solvent interaction
    if ( nsv == ns ) return 1;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( uv.infdil ) {
          // turn on only solvent-solute interaction
          // shut down the solvent-solvent interaction
          uv.prmask[i*ns + j] = ((i < nsv && j >= nsv)
                               || (j < nsv && i >= nsv));
        } else {
          // keep the solvent-solvent interaction
          uv.prmask[i*ns + j] = (i < nsv || j < nsv);
        }
    uv.npr = nsv * (ns - nsv) + ( !uv.infdil ? nsv * (nsv + 1) / 2 : 0 );
    console.log("turning on solute-solvent interaction");
    uv.stage = SOLUTE_SOLVENT;
  } else if ( uv.stage == SOLUTE_SOLVENT ) { // turn on solute-solute interaction
    if ( uv.uutype == DOUU_NEVER ) return 1;
    if ( uv.uutype == DOUU_ATOMIC && !uv.atomicsolute ) return 1;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( uv.infdil ) { // only update u-u interactions
          uv.prmask[i*ns + j] = (i >= nsv && j >= nsv);
        } else { // update all interactions
          uv.prmask[i*ns + j] = 1;
        }
    uv.npr = ( !uv.infdil ? nsv * (nsv + 1) / 2 : ns * (ns + 1) / 2 );
    uv.stage = SOLUTE_SOLUTE;
    console.log("turning on solute-solute interaction");
  } else {
    uv.stage = LAST_STAGE;
    return 1;
  }
  return 0;
}

