/* 
   Function to calculate Chinchilli correlation coeficients (Chinchilli et al., 
   Biometrical Journal 47(5), 644-653, 2005).
   
   With some modifications (in 18/09/2012) from Roberto Hirata Jr / IME-USP.
   
   Gustavo H. Esteves & Raydonal Ospina
   sáb 12 out 2024 19:42:16 -03
*/


#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
//#include <gsl/gsl_math.h>
//#include "gsl_math.h"
//#include "stats.h"
/*
int sign(double x) {
  int res;
  if (x < 0) {
    res = -1;
  else {
    res = 1;
  }
  return res;
}
*/

double cor(double *vx, double *vy, int n, double gamma) {
  /* In: vx vector of doubles
     In: vy vector of doubles
     In: n number of elements of the vectors
     In: gamma 
  */

  int i, j;
  double constante, gx, gy, deltax, deltay, tmp1=0, tmp2=0, tmp3=0;
  
  //constante = (double) 2/(n*(n-1));
  
  for (i=0; i<(n-1); i++) {
    for(j=(i+1); j<n; j++) {
      deltax = vx[i]-vx[j] ;
      deltay = vy[i]-vy[j] ;
      gx = sign(deltax)*R_pow(fabs(deltax), gamma);
      gy = sign(deltay)*R_pow(fabs(deltay), gamma);
      tmp1 += gx*gy;
      tmp2 += gx*gx;
      tmp3 += gy*gy;
    }
  }

  /*
  tmp1 = constante*tmp1;
  tmp2 = constante*tmp2;
  tmp3 = constante*tmp3;
  */

  return (double) tmp1/sqrt(tmp2*tmp3);

}


void corCh(double *mx, int *nr, int *nc, double *gamma, double *res) {
  /* In: mx matrix
     In: nr number of rows --- porque eh um ponteiro?
     In: nc number of columns --- porque eh um ponteiro?
     In: gamma
     Out: resulting matrix
  */

    int i, j, k;
    double tmp, *vtmp1, *vtmp2;
    vtmp1 = (double *) R_alloc((*nr), sizeof(double));
    vtmp2 = (double *) R_alloc((*nr), sizeof(double));
    
    for (i=0; i<(*nc-1); i++) {
	  for(j=(i+1); j<*nc; j++) {
	    for(k=0; k<*nr; k++) {
		    vtmp1[k] = mx[((*nr)*i)+k];
		    vtmp2[k] = mx[((*nr)*j)+k];
	    }
	    tmp = cor(vtmp1, vtmp2, *nr, *gamma);
	    res[((*nc)*i)+j] = tmp;
	    res[((*nc)*j)+i] = tmp;
	  }  
    }
}
