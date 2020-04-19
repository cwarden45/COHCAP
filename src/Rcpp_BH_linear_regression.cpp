// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/students_t.hpp>

using namespace Rcpp;
using namespace boost::math;

//Based upon calculation for pBeta1 from https://github.com/bcaffo/courses/blob/master/07_RegressionModels/01_07_inference/index.Rmd
// [[Rcpp::export]]
long double lmResidual_cpp_1var(NumericVector full_beta, NumericVector full_continuous) {
		int N = full_beta.size();

		if(N < 3){
			return(1);
		}

		if(sd(full_continuous) == 0){
			return(1);
		}

		if(sd(full_beta) == 0){
			return(1);
		}
				
		//beta or methylation = y, continuous variable = x (which would matter if this was multi-variate)
		
		double meanY = mean(full_beta);
		double meanX = mean(full_continuous);

		double sdY = sd(full_beta);
		double sdX = sd(full_continuous);

		//calculate correlation coefficient
		double cov = 0;
		//covariance for loop
		for (int i = 0; i < N; i++){
			cov += (full_beta[i] - meanY)*(full_continuous[i] - meanX);
		}
		double cor = cov / (sdY * sdX);
		
		//calculate slope (beta1)
		double slope =  cor * sdY / sdX;
		
		//calculate intercept (beta0)
		double intercept = meanY - slope * meanX;
		
		//calculuate error / residual (e)
		NumericVector error (N);
		//error for loop
		for (int i = 0; i < N; i++){
			error[i] = full_beta[i] - intercept - slope * full_continuous[i];
		}

		//follow p-value caluation from example
		double sumErrorSquared = 0;
		//sigma for loop
		for (int i = 0; i < N; i++){
			sumErrorSquared += error[i] * error[i];
		}
		double sigma = sqrt(sumErrorSquared / (N-2));

		double SSX = 0;
		//SSX for loop
		for (int i = 0; i < N; i++){
			SSX += (full_continuous[i] - meanX)*(full_continuous[i] - meanX);
		}

		double seSlope = sigma / sqrt(SSX);
		double tSlope = slope / seSlope;
		
		int df = N - 2;
		students_t dist(df);
		double pSlope = 1- cdf(dist, fabs(tSlope));

		pSlope = 2*pSlope;
		if(pSlope <= 1){
			return(pSlope);
		}else{
			return(1);
		}		
		return(pSlope);
}//end def ANOVA_cpp_2group
