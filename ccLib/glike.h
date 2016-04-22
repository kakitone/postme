#include <math.h>

double find_y(double chat, int ell, int L );

double compute_log_likelihood(double c1, int ell1,  double c2, int ell2, 
	int L, double c1hat, double c2hat);

double find_MLE_fit(double c1hat, double c2hat, int ell1, int ell2, int L, 
	double eps, double interval_begin, double interval_end, int is_inside_interval);

double compute_log_ratio(int ell1,int ell2,double c1hat, double c2hat, int L);


