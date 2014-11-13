#ifndef DEBUG_H__
#define DEBUG_H__




/* Debugging routines */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>



/* Print routines */



__inline static void printmat(double *m, int n, const char *name)
{
  int i, j;

  if ( name != NULL ) printf("%s:\n", name);
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      printf("%9g ", m[i*n+j]);
    printf("\n");
  }
  printf("\n");
}



#endif

