#ifndef UV_H__
#define UV_H__



enum { SOLVENT_SOLVENT, SOLUTE_SOLVENT, SOLUTE_SOLUTE, LAST_STAGE };



/* structure that controls the solvent/solute interaction */
typedef struct {
  int *prmask;
  int ns, nsv; /* numbers of sites and solvent sites */
  int stage;
  int npr; /* number of active pairs */
  int infdil; /* infinite dilution */
  int atomicsolute; /* solute molecules are atomic, not molecular */
  int douu; /* how to do solute-solute interaction */
  int uu1step; /* stop after a single uu step */
} uv_t;



/* compute the number solvents
 * here, solvent are defined as molecules with nonzero density
 * zero-density molecules must appear at the end of model */
static int getnsv(model_t *m)
{
  int i;

  if ( m->douu == DOUU_ALLSOLVENT ) /* explicitly set all atoms as the solvent */
    return m->ns;
  for ( i = 0; i < m->ns; i++ )
    if ( m->rho[i] <= 0 )
      break;
  return i;
}



static uv_t *uv_open(model_t *m)
{
  uv_t *uv;
  int i, j, ipr, ns = m->ns;

  xnew(uv, 1);
  uv->ns = ns;
  uv->nsv = getnsv(m);
  xnew(uv->prmask, ns * ns);
  for ( i = 0; i < ns; i++ )
    for ( j = 0; j < ns; j++ )
      uv->prmask[i*ns + j] = (i < uv->nsv && j < uv->nsv);
  uv->npr = uv->nsv * (uv->nsv + 1) / 2;
  uv->stage = SOLVENT_SOLVENT;
  for ( i = uv->nsv; i < uv->ns; i++ )
    if ( m->rho[i] > 0 )
      break;
  uv->infdil = (i == uv->ns);
  /* check if the solute is atomic */
  uv->atomicsolute = 1;
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns && uv->atomicsolute; j++, ipr++ )
      if ( i >= uv->nsv && j >= uv->nsv && m->dis[ipr] > 0 ) {
        uv->atomicsolute = 0;
        break;
      }
  uv->douu = m->douu;
  return uv;
}



/* switch between stages
 * solvent-solvent --> solvent-solute --> solute-solute
 * return 1 if the iteration should stop */
static int uv_switch(uv_t *uv)
{
  int i, j, ns = uv->ns, nsv = uv->nsv;

  /* switch between stages */
  if ( uv->stage == SOLVENT_SOLVENT ) { /* turn on solute-solvent interaction */
    if ( nsv == ns ) return 1;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( uv->infdil ) {
          /* turn on only solvent-solute interaction
           * shut down the solvent-solvent interaction */
          uv->prmask[i*ns + j] = ((i < nsv && j >= nsv)
                               || (j < nsv && i >= nsv));
        } else {
          /* keep the solvent-solvent interaction */
          uv->prmask[i*ns + j] = (i < nsv || j < nsv);
        }
    uv->npr = nsv * (ns - nsv) + ( !uv->infdil ? nsv * (nsv + 1) / 2 : 0 );
    fprintf(stderr, "turning on solute-solvent interaction\n");
    uv->stage = SOLUTE_SOLVENT;
  } else if ( uv->stage == SOLUTE_SOLVENT ) { /* turn on solute-solute interaction */
    if ( uv->douu == DOUU_NEVER ) return 1;
    if ( uv->douu == DOUU_ATOMIC && !uv->atomicsolute ) return 1;
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( uv->infdil ) {
          /* only update u-u interactions */
          uv->prmask[i*ns + j] = (i >= nsv && j >= nsv);
        } else {
          /* update all interactions */
          uv->prmask[i*ns + j] = 1;
        }
    uv->npr = !uv->infdil ? nsv*(nsv + 1)/2 : ns*(ns + 1)/2;
    uv->stage = SOLUTE_SOLUTE;
    uv->uu1step = (uv->infdil && uv->atomicsolute && uv->douu != DOUU_ALWAYS);
    if ( uv->uu1step ) /* stop after a single uu step */
      return 1;
    fprintf(stderr, "turning on solute-solute interaction\n");
  } else {
    uv->stage = LAST_STAGE;
    return 1;
  }
  return 0;
}



static void uv_close(uv_t *uv)
{
  free(uv->prmask);
  free(uv);
}



#endif /* UV_H__ */
