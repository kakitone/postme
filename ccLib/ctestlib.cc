#include "ctestlib.h"

CTestRobot::CTestRobot(double cov1_in, double cov2_in, int l1_in, int l2_in){
	CTestRobot::setParameters(cov1_in,cov2_in, l1_in, l2_in)	;
}

void CTestRobot::setParameters(double cov1_in, double cov2_in, int l1_in, int l2_in){
	cov1 = cov1_in;  
	cov2 = cov2_in;
	l1 = l1_in; 
	l2 = l2_in; 
}
	
double CTestRobot::computePValue(){

	double read_length = 150;
    double read_count_1, read_count_2 ; 
    
    // Preprocessing to get read counts 
    read_count_1 = cov1*l1/read_length;
    read_count_2 = cov2*l2/read_length; 

    // CTest steps
    int xTmp, nTmp; 
    double  pTmp, pvalue;
    
    xTmp = (int) read_count_1;
    nTmp = (int) (read_count_1 + read_count_2);
    pTmp = l1*1.0/(l1 + l2);

    double tmp = cdf(binomial(nTmp, pTmp), xTmp); 
	pvalue = 2*std::min(1 - tmp, tmp);

    return pvalue;
}

