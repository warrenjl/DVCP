#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List eta_update(int k_approx,
                      arma::mat cross_a_approx,
                      arma::mat eta_approx_corr_inv,
                      int h_model,
                      arma::vec distance_to_ps,
                      arma::vec angle_key,
                      arma::mat x,
                      arma::vec indicator,
                      arma::vec omega,
                      arma::vec kappa,
                      arma::vec beta,
                      double theta,
                      double lambda,
                      arma::vec eta,
                      arma::vec eta_approx,
                      double sigma2_eta,
                      double phi_eta,
                      arma::vec metrop_var_eta,
                      arma::vec acctot_eta){

arma::vec eta_approx_old = eta_approx;

for(int j = 0; j < k_approx; ++j){
  
   /*Second*/
   arma::vec eta_old = eta;
   arma::vec indicator_old = indicator;
   
   double second = -0.50*dot((kappa - x*beta - theta*indicator_old), (omega%(kappa - x*beta - theta*indicator_old))) -
                   0.50*dot(eta_approx, (eta_approx_corr_inv*eta_approx))/sigma2_eta;
  
   /*First*/
   eta_approx(j) = R::rnorm(eta_approx_old(j), 
                            sqrt(metrop_var_eta(j)));
   eta = exp(-phi_eta*cross_a_approx)*(eta_approx_corr_inv*eta_approx);
   indicator = indicator_fun(h_model,
                             distance_to_ps,
                             angle_key,
                             lambda,
                             eta);
   
   double first = -0.50*dot((kappa - x*beta - theta*indicator), (omega%(kappa - x*beta - theta*indicator))) -
                  0.50*dot(eta_approx, (eta_approx_corr_inv*eta_approx))/sigma2_eta;
  
   /*Decision*/
   double ratio = exp(first - second);   
   int acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
    
     eta_approx(j) = eta_approx_old(j);
     eta = eta_old;
     indicator = indicator_old;
     acc = 0;
    
     }
   acctot_eta(j) = acctot_eta(j) + 
                   acc;
  
   }
  
return Rcpp::List::create(Rcpp::Named("eta") = eta,
                          Rcpp::Named("eta_approx") = eta_approx,
                          Rcpp::Named("indicator") = indicator,
                          Rcpp::Named("acctot_eta") = acctot_eta);
  
}

