/*Funcao para calcular a soma da correlacao de Chinchilli*/

/*Comandos:
 R CMD SHLIB somac.c
 dyn.load("somac.so")
 .C("chinc",...)*/

#include <R.h>
#include <gsl/gsl_math.h>

void chinc(double *x,double *y, int *n, double *gamma, double *coeficiente_chinchilli)
{
   double Uxx=0.0,Uyy=0.0,Uxy=0.0;
   double gx,gy;
   double resultado1=0.0;
   double resultado2=0.0;

   //constante=(double) 2/(*n *(*n-1));	

	for (int i=0; i<(*n-1); i++){
		for (int j=i+1; j<*n; j++){
                  resultado1 = x[i]-x[j];
                  resultado2 = y[i]-y[j];
                  gx = GSL_SIGN(resultado1)*pow(fabs(resultado1),*gamma);
                  gy = GSL_SIGN(resultado2)*pow(fabs(resultado2),*gamma);
                  Uxx += gx*gx;
		              Uyy += gy*gy;
		              Uxy += gx*gy;
		}			
	}


	*coeficiente_chinchilli=Uxy/sqrt(Uxx*Uyy);
}
