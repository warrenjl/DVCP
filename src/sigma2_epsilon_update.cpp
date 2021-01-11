#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x,
                             arma::mat indicator,
                             arma::vec beta,
                             double theta,
                             double alpha_sigma2_epsilon,
                             double beta_sigma2_epsilon){

int n = y.size();
double alpha_sigma2_epsilon_update = 0.50*n + 
                                     alpha_sigma2_epsilon;

double beta_sigma2_epsilon_update = 0.50*dot((y - x*beta - theta*indicator), (y - x*beta - theta*indicator)) + 
                                    beta_sigma2_epsilon;

double sigma2_epsilon = 1/R::rgamma(alpha_sigma2_epsilon_update,
                                    (1.00/beta_sigma2_epsilon_update));

return(sigma2_epsilon);

}





