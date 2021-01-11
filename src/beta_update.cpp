#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_update(arma::mat x, 
                      arma::vec indicator,
                      arma::vec omega,
                      arma::vec kappa,
                      double theta,
                      double sigma2_beta){

int p_x = x.n_cols;
int n = omega.size();

arma::mat omega_mat(n, p_x);
for(int j = 0; j < p_x; ++j){
   omega_mat.col(j) = omega;
   }

arma::mat x_trans = trans(x);

arma::mat cov_beta = inv_sympd(x_trans*(omega_mat%x) + 
                               eye(p_x, p_x)/sigma2_beta);

arma::vec mean_beta = cov_beta*(x_trans*(omega%(kappa - theta*indicator)));

arma::mat ind_norms = arma::randn(1, p_x);
arma::vec beta = mean_beta + 
                 trans(ind_norms*arma::chol(cov_beta));

return beta;

}



