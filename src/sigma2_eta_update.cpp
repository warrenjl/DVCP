#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_eta_update(arma::mat eta_approx_corr_inv,
                         arma::vec eta_approx,
                         double alpha_sigma2_eta,
                         double beta_sigma2_eta){

int p_eta = eta_approx.size();
  
double alpha_sigma2_eta_update = 0.50*p_eta + 
                                 alpha_sigma2_eta;

double beta_sigma2_eta_update = 0.50*dot(eta_approx, (eta_approx_corr_inv*eta_approx)) + 
                                beta_sigma2_eta;

double sigma2_eta = 1/R::rgamma(alpha_sigma2_eta_update,
                                (1.00/beta_sigma2_eta_update));

return(sigma2_eta);

}





