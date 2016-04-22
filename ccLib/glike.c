#include "glike.h"  /* Include the header (not strictly necessary here) */


double find_y(double chat, int ell, int L ){
	return chat*ell/L ;
}

double compute_log_likelihood(double c1, int ell1,  double c2, int ell2, 
	int L, double c1hat, double c2hat){
	
	double y1hat = find_y(c1hat, ell1, L), y2hat = find_y(c2hat, ell2, L) ; 
	double y1 = c1*ell1*1.0/L ,  y2 = c2*ell2*1.0/L ;

	return -y1 - y2 + y1hat*log(y1) + y2hat*log(y2);
}

double find_MLE_fit(double c1hat, double c2hat, int ell1, int ell2, int L, 
	double eps, double interval_begin, double interval_end, int is_inside_interval){


	int check = 0 ;

	if (is_inside_interval && (interval_begin < c1hat/c2hat) && (c1hat/c2hat < interval_end) ) {
		check = 1 ; 
	}else if ((! is_inside_interval ) && ( interval_begin > c1hat/c2hat || interval_end < c1hat/c2hat)){
		check = 1;
	}



	if (check){
		return compute_log_likelihood(c1hat,  ell1,  c2hat, ell2, L, c1hat, c2hat);
	}
	else{
		double c1, c2 , ml1, ml2; 

		c2 = (c1hat*ell1 + c2hat*ell2)*1.0/((1-eps)*ell1 + ell2);
		c1 = c2*(1- eps);
		ml1 = compute_log_likelihood(c1, ell1, c2, ell2, L, c1hat, c2hat);

		c2 = (c1hat*ell1 + c2hat*ell2)*1.0/((1+eps)*ell1 + ell2);
		c1 = c2*(1+ eps);
		ml2 = compute_log_likelihood(c1, ell1, c2, ell2, L, c1hat, c2hat);
		
		return  ml1 > ml2 ? ml1 : ml2 ;	
	}

}

double compute_log_ratio(int ell1,int ell2,double c1hat, double c2hat, int L){
	double eps = 0.05 ;
	double interval_begin = 1- eps, interval_end = 1+ eps; 

	return find_MLE_fit(c1hat, c2hat, ell1, ell2, L,  eps, interval_begin,  interval_end, 1) 
			- find_MLE_fit(c1hat, c2hat, ell1, ell2, L,  eps, interval_begin,  interval_end, 0)  ;

}



























