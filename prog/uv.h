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
  int simplesolvent; /* solvent molecules are atomic, not molecular */
} uv_t;



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
  uv->stage = 0;
  for ( i = uv->nsv; i < uv->ns; i++ )
    if ( m->rho[i] > 0 )
      break;
  uv->infdil = (i == uv->ns);
  /* check if the solute is simple */
  uv->simplesolvent = 1;
  for ( ipr = 0, i = 0; i < ns; i++ )
    for ( j = i + 1; j < ns; j++, ipr++ )
      if ( i >= uv->nsv && j >= uv->nsv && m->dis[ipr] > 0 ) {
        uv->simplesolvent = 0;
        break;
      }
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
    for ( i = 0; i < ns; i++ )
      for ( j = 0; j < ns; j++ )
        if ( uv->infdil ) {
          uv->prmask[i*ns + j] = (i >= nsv && j >= nsv);
        } else {
          uv->prmask[i*ns + j] = 1;
        }
    uv->npr = ( !uv->infdil ? nsv * (nsv + 1) / 2 : ns * (ns + 1) / 2 );
    uv->stage = SOLUTE_SOLUTE;
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
