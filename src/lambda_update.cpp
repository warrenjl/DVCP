#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List lambda_update(int h_model,
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
                         double a_lambda,
                         double b_lambda,
                         double metrop_var_lambda,
                         int acctot_lambda){

/*Second*/
double lambda_old = lambda;
double lambda_trans_old = log(lambda_old/(1.00 - lambda_old));
arma::vec indicator_old = indicator;

double second = -0.50*dot((kappa - x*beta - theta*indicator), (omega%(kappa - x*beta - theta*indicator))) +
                R::dbeta(lambda_old,
                         a_lambda,
                         b_lambda,
                         TRUE) +
                lambda_trans_old -
                2.00*log(1.00 + exp(lambda_trans_old));

/*First*/
double lambda_trans = R::rnorm(lambda_trans_old, 
                               sqrt(metrop_var_lambda));
lambda = exp(lambda_trans)/(1.00 + exp(lambda_trans));
indicator = indicator_fun(h_model,
                          distance_to_ps,
                          angle_key,
                          lambda,
                          eta);

double first = -0.50*dot((kappa - x*beta - theta*indicator), (omega%(kappa - x*beta - theta*indicator))) +
               R::dbeta(lambda,
                        a_lambda,
                        b_lambda,
                        TRUE) +
               lambda_trans -
               2.00*log(1.00 + exp(lambda_trans));

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  
  lambda = lambda_old;
  indicator = indicator_old;
  acc = 0;
  
  }
acctot_lambda = acctot_lambda + 
                acc;

return Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                          Rcpp::Named("indicator") = indicator,
                          Rcpp::Named("acctot_lambda") = acctot_lambda);

}
                 
  
