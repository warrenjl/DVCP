#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double theta_update(arma::mat x, 
                    arma::vec indicator,
                    arma::vec omega,
                    arma::vec kappa,
                    arma::vec beta,
                    double sigma2_theta){

double var_theta = 1.00/(dot(indicator, (omega%indicator)) + 
                         1.00/sigma2_theta);

double mean_theta = var_theta*dot(indicator, (omega%(kappa - x*beta)));

double theta = R::rnorm(mean_theta,
                        sqrt(var_theta));

return theta;

}



