var SOLVENT_SOLVENT = 0;
var SOLUTE_SOLVENT = 1;
var SOLUTE_SOLUTE = 2;
var LAST_STAGE = 3;



/* compute the number solvents
 * here, solvent are defined as molecules with nonzero density
 * zero-density molecules must appear at the end of model */
function getnsv(uutype)
{
  if ( uutype == "All-solvent" ) // explicitly set all atoms as the solvent
    return ns;
  for ( var i = 0; i < ns; i++ )
    if ( rho[i] <= 0 ) return i;
  return i;
}



function UV(ns)
{
  var i, j, ipr;

  this.ns = ns;
  this.prmask = newarr(ns*ns);
  this.uutype = grab("douu").value;
  this.nsv = getnsv(this.uutype); // number of solvent sites
  this.stage = SOLVENT_SOLVENT;
  for ( i = 0; i < ns; i++ )
    for ( j = 0; j < ns; j++ )
      this.prmask[i*ns + j] = (i < this.nsv && j < this.nsv);
  this.npr = this.nsv * (this.nsv + 1) / 2; // number of active pairs

  for ( i = this.nsv; i < this.ns; i++ )
    if ( rho[i] > 0 )
      break;
  this.infdil = (i == this.ns); // infinite dilution

  // see if the solute is atomic
  this.atomicsolute = true;
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns && this.atomicsolute; j++, ipr++ )
      if ( i >= this.nsv && j >= this.nsv && dis[i*ns + j] > 0 ) {
        this.atomicsolute = false;
        break;
      }
}



/* switch between stages
 * solvent-solvent --> solvent-solute --> solute-solute
 * return 1 if the iteration should stop */
UV.prototype.switchstage = function()
{
  var i, j, nsv = this.nsv;

  /* switch between stages */
  if ( this.stage == SOLVENT_SOLVENT ) { // turn on solute-solvent interaction
    if ( nsv == ns ) return 1;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( this.infdil ) {
          // turn on only solvent-solute interaction
          // shut down the solvent-solvent interaction
          this.prmask[i*ns + j] = ((i < nsv && j >= nsv)
                               || (j < nsv && i >= nsv));
        } else {
          // keep the solvent-solvent interaction
          this.prmask[i*ns + j] = (i < nsv || j < nsv);
        }
    this.npr = nsv * (ns - nsv) + ( !this.infdil ? nsv * (nsv + 1) / 2 : 0 );
    console.log("turning on solute-solvent interaction");
    this.stage = SOLUTE_SOLVENT;
  } else if ( this.stage == SOLUTE_SOLVENT ) { // turn on solute-solute interaction
    if ( this.uutype == "Never" ) return 1;
    if ( this.uutype == "Atomic" && !this.atomicsolute ) return 1;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( this.infdil ) { // only update u-u interactions
          this.prmask[i*ns + j] = (i >= nsv && j >= nsv);
        } else { // update all interactions
          this.prmask[i*ns + j] = 1;
        }
    this.npr = ( !this.infdil ? nsv * (nsv + 1) / 2 : ns * (ns + 1) / 2 );
    this.stage = SOLUTE_SOLUTE;
    console.log("turning on solute-solute interaction");
  } else {
    this.stage = LAST_STAGE;
    return 1;
  }
  return 0;
}

