// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/students_t.hpp>

using namespace Rcpp;
using namespace boost::math;

//Based upon code from http://www.boost.org/doc/libs/1_65_0/libs/math/doc/html/math_toolkit/stat_tut/weg/st_eg/two_sample_students_t.html
// and https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes.2C_unequal_variances

// [[Rcpp::export]]
long double ttest_cpp(NumericVector beta1, NumericVector beta2) {
	unsigned Sn1 = beta1.size(), Sn2 = beta2.size();
	if((Sn1 < 3) ||(Sn2 < 3)){
		return(1);
	}
		
	double Sm1 = mean(beta1);
	double Sd1 = sd(beta1);
	double Sm2 = mean(beta2);
	double Sd2 = sd(beta2);
		
	//assume unequal variance
	double t_stat = (Sm1 - Sm2) / sqrt(Sd1/Sn1 + Sd2/Sn2);
	double df = ((Sd1/Sn1 + Sd2/Sn2)*(Sd1/Sn1 + Sd2/Sn2)) /
				((Sd1/Sn1)*(Sd1/Sn1)/(Sn1-1) + (Sd2/Sn2)*(Sd2/Sn2)/(Sn2-1));
	students_t dist(df);
	double p1 = 1- cdf(dist, fabs(t_stat));
	//this code was faster than custom t-test function, but using fastLm() + dt() was slight faster than using fastLm() with Boost p-value calculation
	//so, keep this code if using t-test, but just use fastlm() and pt() for linear regression
	//however, without assuming equal variance, c++ p-value tends to be lower (even without doubling for two-sided test)
	p1 = 2*p1;
	if(p1 <= 1){
		return(p1);
	}else{
		return(1);
	}
}//end def ttest_cpp

// [[Rcpp::export]]
long double ttest_cpp_paired(NumericVector paired_diff) {
	//unlike fastLm code, this is faster than calculating within R using pt(), or with t.test
	unsigned Sn = paired_diff.size();
	if(Sn < 2){
		return(1);
	}
		
	double t_stat = mean(paired_diff) / sqrt(sd(paired_diff)/Sn);
	double df = Sn-1;
	students_t dist(df);
	double p1 = 1- cdf(dist, fabs(t_stat));

	p1 = 2*p1;
	if(p1 <= 1){
		return(p1);
	}else{
		return(1);
	}
}//end def ttest_cpp_paired