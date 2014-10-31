#ifndef MODEL_H__
#define MODEL_H__



/* input parameters */



#define MAXATOM 6 /* maximal number of atoms/sites */

enum { HARD_SPHERE, LJ_FULL, LJ_REPULSIVE };
enum { IE_PY, IE_HNC };
enum { SOLVER_PICARD, SOLVER_LMV, SOLVER_MDIIS };

typedef struct {
  int ns;
  double sigma[MAXATOM];
  double eps6_12[MAXATOM][2];
  double C6_12[MAXATOM*(MAXATOM+1)/2][2]; /* alternative to sigma/epsilon */
  double rho[MAXATOM];
  double Lpm[MAXATOM*(MAXATOM-1)/2];
  double beta;
  int ljtype;

  double charge[MAXATOM];
  double ampch;
  double radch; /* screening distance */

  int ietype;
  double rmax; /* radius cutoff */
  int npt; /* number of sampling points along r */

  int nlambda; /* number of intermediate stages */
  int itmax; /* maximial number of iterations in each stage */
  double tol; /* tolerance of error */
  int solver; /* solver */

  double picdamp; /* damping factor for the Picard solver */

  int Mpt; /* number of points for the Newton-Raphson method (LMV) */
  double lmvdamp; /* damping factor for the LMV solver */

  int nbases; /* number of bases for the MDIIS solver */
  double mdiisdamp; /* damping factor for the MDIIS solver */
} model_t;



#endif /* MODEL_H__ */
