#include <boost/math/distributions/binomial.hpp>
using ::boost::math::cdf;
using ::boost::math::binomial;

class CTestRobot {

public:  
	CTestRobot(double cov1_in, double cov2_in, int l1_in, int l2_in);	
	double computePValue();
	void setParameters(double cov1_in, double cov2_in, int l1_in, int l2_in);
	
private:
	double cov1, cov2;
	int l1, l2; 
};
