// Usage : gcc -o glikemain  glikemain.c glike.c ; ./glikemain


#include <stdio.h>
#include "glike.h" 


int main(void)
{
 	double L = 100;
 	int i, j ; 
 	double cov1, cov2;
 	int ell1, ell2 ;
 	double logPdiff;

	for (j = 0 ; j<5 ; j++){
		for (i= 0; i< 5; i++){
			cov1 = 50.0  ;
			cov2 = 50.0 + 2.0*i ; 
			ell1 = 30*pow(10,j); 
			ell2 = 30*pow(10,j);

			logPdiff = compute_log_ratio(ell1, ell2, cov1, cov2, L);
			printf("cov1, cov2, l1, l2, logPdiff: %7f \t %7f \t %7d \t %7d \t %7f \n", cov1, cov2, ell1, ell2, logPdiff);
		} 
			
		printf("\n"); 
	}
		

    return 0;
}