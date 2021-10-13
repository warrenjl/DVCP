#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List omega_update(arma::vec y,
                        arma::mat x,
                        arma::vec indicator,
                        arma::vec beta,
                        double theta){
  
arma::vec mean_omega = x*beta +
                       theta*indicator;

arma::vec input(1); input.fill(1.00);
arma::vec omega = rcpp_pgdraw(input,
                              mean_omega);

arma::vec kappa = (y - 0.50)/omega;

return Rcpp::List::create(Rcpp::Named("omega") = omega,
                          Rcpp::Named("kappa") = kappa);

}
































































