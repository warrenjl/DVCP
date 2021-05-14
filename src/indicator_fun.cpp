#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec indicator_fun(int h_model,
                        arma::vec distance_to_ps,
                        arma::vec angle_key,
                        double lambda,
                        arma::vec eta){
   
int n = distance_to_ps.size();
arma::vec indicator_temp(n); indicator_temp.fill(0.00);
double eta_full = 0.00;
for(int j = 0; j < n; ++j){
   
   eta_full = eta(angle_key(j) - 1);
   indicator_temp(j) = (distance_to_ps(j) <= lambda*exp(eta_full));
   if(distance_to_ps(j) == 0.00){
     indicator_temp(j) = 1.00;
     }
   
   }

arma::vec indicator(n); indicator.fill(0.00);
if(h_model == 0){ //Indicator
  indicator = indicator_temp;
  }
   
if(h_model == 1){ //Linear
  indicator = (1.00 - distance_to_ps)%indicator_temp;
  }
   
if(h_model == 2){ //Exponential
  indicator = exp(-distance_to_ps)%indicator_temp;
  }
   
if(h_model == 3){ //Gaussian
  indicator = exp(-(distance_to_ps%distance_to_ps))%indicator_temp;
  }

if(h_model == 4){ //Spherical
  indicator = (1.00 +
               -1.50*distance_to_ps/exp(eta_full) +
               0.50*pow((distance_to_ps/exp(eta_full)), 3))*(distance_to_ps < exp(eta_full));
  }
 
return indicator;
    
}
   
   