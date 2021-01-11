#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_eta_update(arma::mat d_a_approx,
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
                          double alpha_phi_eta,
                          double beta_phi_eta,
                          double metrop_var_phi_eta,
                          int acctot_phi_eta){

double sign = 0.00;
  
/*Second*/
double phi_eta_old = phi_eta;
double phi_eta_trans_old = log(phi_eta_old);
arma::mat eta_approx_corr_inv_old = eta_approx_corr_inv;
arma::vec eta_old = eta;
arma::vec indicator_old = indicator;
double log_deter_old = 0.00;
log_det(log_deter_old, sign, eta_approx_corr_inv_old);

double second = -0.50*dot((kappa - x*beta - theta*indicator_old), (omega%(kappa - x*beta - theta*indicator_old))) + 
                0.50*log_deter_old - 
                0.50*dot(eta_approx, (eta_approx_corr_inv_old*eta_approx))/sigma2_eta + 
                R::dgamma(phi_eta_old,
                          alpha_phi_eta,
                          (1.00/beta_phi_eta),
                          TRUE) +
                phi_eta_trans_old;

/*First*/
double phi_eta_trans = R::rnorm(phi_eta_trans_old, 
                                sqrt(metrop_var_phi_eta));
phi_eta = exp(phi_eta_trans);
eta_approx_corr_inv = inv_sympd(exp(-phi_eta*d_a_approx));
eta = exp(-phi_eta*cross_a_approx)*(eta_approx_corr_inv*eta_approx);
indicator = indicator_fun(h_model,
                          distance_to_ps,
                          angle_key,
                          lambda,
                          eta);
double log_deter = 0.00;
log_det(log_deter, sign, eta_approx_corr_inv);

double first = -0.50*dot((kappa - x*beta - theta*indicator), (omega%(kappa - x*beta - theta*indicator))) + 
               0.50*log_deter - 
               0.50*dot(eta_approx, (eta_approx_corr_inv*eta_approx))/sigma2_eta + 
               R::dgamma(phi_eta,
                         alpha_phi_eta,
                         (1.00/beta_phi_eta),
                         TRUE) +
               phi_eta_trans;

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){

  phi_eta = phi_eta_old;
  eta_approx_corr_inv = eta_approx_corr_inv_old;
  eta = eta_old;
  indicator = indicator_old;
  acc = 0;
  
  }
acctot_phi_eta = acctot_phi_eta + 
                 acc;

return Rcpp::List::create(Rcpp::Named("phi_eta") = phi_eta,
                          Rcpp::Named("eta_approx_corr_inv") = eta_approx_corr_inv,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("indicator") = indicator,
                          Rcpp::Named("acctot_phi_eta") = acctot_phi_eta);

}
                 
  
