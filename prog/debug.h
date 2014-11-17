#ifndef DEBUG_H__
#define DEBUG_H__




/* Debugging routines */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>



/* Print routines */


#define printmat(m, n, name) fprintmat(stderr, m, n, name)

__inline static void fprintmat(FILE *fp, double *m, int n, const char *name)
{
  int i, j;

  if ( name != NULL ) fprintf(fp, "%s:\n", name);
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      fprintf(fp, "%+10.3e ", m[i*n+j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
}



#endif

