// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/fisher_f.hpp>

using namespace Rcpp;
using namespace boost::math;

//Based upon code from http://www.boost.org/doc/libs/1_65_0/libs/math/doc/html/math_toolkit/stat_tut/weg/f_eg.html
// and https://en.wikipedia.org/wiki/Analysis_of_variance and http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_HypothesisTesting-ANOVA/BS704_HypothesisTesting-Anova3.html

// [[Rcpp::export]]
long double ANOVA_cpp_2group(NumericVector full_beta, NumericVector betaT, NumericVector betaR) {
		int N = full_beta.size(), n1 = betaT.size(), n2 = betaR.size();
		int k = 2;
		if((n1 ==0) ||(n2 == 0)){
			return(1);
		}

		if(sd(full_beta) == 0){
			return(1);
		}
		
		double overall_mean = mean(full_beta);
		double m1 = mean(betaT);
		double m2 = mean(betaR);

		double SST = n1*(m1 - overall_mean)*(m1 - overall_mean) + n2*(m2 - overall_mean)*(m2 - overall_mean);
		double MST = SST / (k-1);
		
		double SS1 = 0;
		//SS1 for loop
		for (int i = 0; i < n1; i++){
			SS1 += (betaT[i] - m1)*(betaT[i] - m1);
		}
		double SS2 = 0;
		//SS2 for loop
		for (int i = 0; i < n2; i++){
			SS2 += (betaR[i] - m2)*(betaR[i] - m2);
		}
		double MSE = (SS1 + SS2) / (N-k);

		if(MSE == 0){
			MSE = 0.00001;
		}//F-statistic can't be infinity, so provide a small value
		
		double f_stat = MST / MSE;
		
		boost::math::fisher_f dist(k-1, N-k);
		double p1 = 1 - cdf(dist, f_stat);
		return(p1);
}//end def ANOVA_cpp_2group

// [[Rcpp::export]]
long double ANOVA_cpp_2group_2way(NumericVector full_beta,
									NumericVector betaT, NumericVector betaR,
									NumericVector interact_var) {
		int N = full_beta.size(), n1 = betaT.size(), n2 = betaR.size();
		int k = 2;
		if((n1 ==0) ||(n2 == 0)){
			return(1);
		}

		if(sd(full_beta) == 0){
			return(1);
		}
		
		double overall_mean = mean(full_beta);
		double m1 = mean(betaT);
		double m2 = mean(betaR);

		double SST = n1*(m1 - overall_mean)*(m1 - overall_mean) + n2*(m2 - overall_mean)*(m2 - overall_mean);
		double MST = SST / (k-1);
		
		double SSE = 0;
		for (int i = 1; i <= max(interact_var); i++){
			NumericVector temp_var;
			double temp_sum = 0;
			//collect values for interaction variable and calculate mean
			for(int j = 0; j < interact_var.size(); j++){
				if(interact_var[j] == i){
					temp_var.push_back(full_beta[j]);
					temp_sum += full_beta[j];
				}// end if(interact_var[j] == i)
			}//end for(int j = 0; j < interact_var.size(); j++)
			double interact_mean = temp_sum / temp_var.size();
		
			//calculate sum of errors
			for(int j = 0; j < temp_var.size(); j++){
				SSE += (temp_var[j] - interact_mean)*(temp_var[j] - interact_mean);
			}//end for(int j = 0; j < interact_var.size(); j++)
		}//end for (int i = 1; i <= max(interact_var); i++)
		double MSE = SSE / (N-max(interact_var));
		
		if(MSE == 0){
			MSE = 0.00001;
		}//F-statistic can't be infinity, so provide a small value
		
		double f_stat = MST / MSE;
		
		//use maximum numeric interaction value as the product of number of levels for each variable
		boost::math::fisher_f dist(k-1, N-max(interact_var));
		double p1 = 1 - cdf(dist, f_stat);
		return(p1);
}//end def ANOVA_cpp_2group_2way