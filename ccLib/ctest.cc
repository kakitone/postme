// You will need Boost library : sudo port install boost 
// g++ ctest.cc ctestlib.cc -o ctest 
// ./ctest 
#include <iostream>
using namespace std;
#include "ctestlib.h"

int main()
{
  CTestRobot crobot(50.0 , 50.0, 300, 300);
  double pvalue = crobot.computePValue(); 
  double cov1_in, cov2_in;
  int l1_in, l2_in ; 
  	
  std::cout << "Sample run p-value : " << pvalue << std::endl;

  for (int j =0; j<= 4; j++){
	  for (int i =0 ; i <= 4 ; i++){
  		cov1_in = 50.0 - 10*i ;
		cov2_in = 52.0 - 10*i;
		l1_in = 300*std::pow(10, j);
		l2_in = 300*std::pow(10, j) ;
		crobot.setParameters(cov1_in, cov2_in, l1_in, l2_in);
		pvalue = crobot.computePValue();
		std::cout << "cov1_in, cov2_in, l1_in, l2_in, pvalue ";
		std::cout << cov1_in <<  "\t" << cov2_in <<  "\t" <<  l1_in <<  "\t" <<  l2_in <<  "\t" <<  pvalue << endl;  

  	}
	std::cout << endl;
  }
  return 0;
}

