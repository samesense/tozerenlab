#include<math.h>


void lpr(double *E,
          double *s2,
          double *t2, 
          double *a, double *b, double *theta,
          int *n, int *m,
          double *th, double *sm)
{
   int i, j;
   double lp0,lp1,e,e1,e2;
   *sm=0.0;
   for(j = 0 ; j < *m ; j++ ) {
       e1=0.0; e2=0.0;
       for(i=0 ; i< *n ; i++ ) { e=E[ i*(*m) + j ] ; e1=e1+e;  e2=e2+e*e ; }
       e1=e1/(*n*1.0) ;  e2=e2/(*n*1.0) ;
       lp0= -.5*(*n)*( log(s2[j]) +  e2/s2[j] ) ;
       lp1= -1.0*(  .5*( (*n)*log(s2[j]) + log(1+ t2[j]*(*n)/(s2[j]))) )  +
          lgamma(*a+ *n/(2.0)) - lgamma(*a) +
          (*a)*log(*b) - (*a+*n/(2.0))*
          log(*b+.5*( (*n)*(e2-e1*e1*t2[j]/(s2[j]/((*n)*1.0) +t2[j]))/s2[j]));
       th[j]= lp1-lp0; 
       *sm=*sm+th[j]+log(exp(-th[j])+exp(theta[j]))-log(1.0+exp(theta[j])) ;
                                }
}

                                                                                

