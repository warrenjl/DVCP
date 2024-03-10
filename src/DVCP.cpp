#include "RcppArmadillo.h"
#include "DVCP.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List DVCP(int mcmc_samples,
                int burnin,
                int thin,
                int adapt,
                int likelihood_indicator, //0: Bernoulli, 1: Gaussian
                int h_model, //0: Indicator; 1: Linear, 2: Exponential; 3: Gaussian; 4: Spherical
                arma::vec approx_angles, //Degrees, Not Radians
                arma::vec y,
                arma::mat x,
                arma::vec distance_to_ps,
                arma::vec unique_angles, //Degrees, Not Radians
                arma::vec angle_key,
                double metrop_var_lambda,
                arma::vec metrop_var_eta,
                double metrop_var_phi_eta,
                double adapt_lambda,
                double adapt_eta,
                double adapt_phi_eta,
                Rcpp::Nullable<Rcpp::NumericVector> trials = R_NilValue,
                Rcpp::Nullable<double> alpha_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<double> beta_sigma2_epsilon_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_theta_prior = R_NilValue,
                Rcpp::Nullable<double> a_lambda_prior = R_NilValue,
                Rcpp::Nullable<double> b_lambda_prior = R_NilValue,
                Rcpp::Nullable<double> alpha_sigma2_eta_prior = R_NilValue,
                Rcpp::Nullable<double> beta_sigma2_eta_prior = R_NilValue,
                Rcpp::Nullable<double> alpha_phi_eta_prior = R_NilValue,
                Rcpp::Nullable<double> beta_phi_eta_prior = R_NilValue,
                Rcpp::Nullable<double> sigma2_epsilon_init = R_NilValue,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                Rcpp::Nullable<double> theta_init = R_NilValue,
                Rcpp::Nullable<double> lambda_init = R_NilValue,
                Rcpp::Nullable<double> sigma2_eta_init = R_NilValue,
                Rcpp::Nullable<double> phi_eta_init = R_NilValue){
  
//Defining Parameters and Quantities of Interest
int n = y.size();
int m = unique_angles.size();
int k_approx = approx_angles.size();
int p_x = x.n_cols;
arma::vec sigma2_epsilon_keep((mcmc_samples - burnin)/thin); sigma2_epsilon_keep.fill(0.00);
arma::mat beta_keep(p_x, (mcmc_samples - burnin)/thin); beta_keep.fill(0.00);
arma::vec theta_keep((mcmc_samples - burnin)/thin); theta_keep.fill(0.00);
arma::vec lambda_keep((mcmc_samples - burnin)/thin); lambda_keep.fill(0.00);
arma::mat eta_keep(m, (mcmc_samples - burnin)/thin); eta_keep.fill(0.00);
arma::vec sigma2_eta_keep((mcmc_samples - burnin)/thin); sigma2_eta_keep.fill(0.00);
arma::vec phi_eta_keep((mcmc_samples - burnin)/thin); phi_eta_keep.fill(0.00);
arma::vec neg_two_loglike_keep((mcmc_samples - burnin)/thin); neg_two_loglike_keep.fill(0.00);

arma::vec tri_als(n); tri_als.fill(1);
if(trials.isNotNull()){
  tri_als = Rcpp::as<arma::vec>(trials);
  }

//Approximation Information
arma::mat d_a_approx(k_approx, k_approx); d_a_approx.fill(0.00);
for(int j = 0; j < k_approx; ++j){
   for(int k = 0; k < k_approx; ++k){
      d_a_approx(j,k) = std::min(abs(approx_angles(j) - approx_angles(k)),
                                 360.00 - std::max(approx_angles(j), approx_angles(k)) + std::min(approx_angles(j), approx_angles(k)));
      }
   }  
d_a_approx = d_a_approx/100.00; //Scaling

arma::mat cross_a_approx(m, k_approx); cross_a_approx.fill(0.00);
for(int j = 0; j < m; ++j){
  for(int k = 0; k < k_approx; ++k){
     cross_a_approx(j,k) = std::min(abs(unique_angles(j) - approx_angles(k)),
                                    360.00 - std::max(unique_angles(j), approx_angles(k)) + std::min(unique_angles(j), approx_angles(k)));
     }
   }  
cross_a_approx = cross_a_approx/100.00; //Scaling
  
//Prior Information
double alpha_sigma2_epsilon = 0.01;
if(alpha_sigma2_epsilon_prior.isNotNull()){
  alpha_sigma2_epsilon = Rcpp::as<double>(alpha_sigma2_epsilon_prior);
  }

double beta_sigma2_epsilon = 0.01;
if(beta_sigma2_epsilon_prior.isNotNull()){
  beta_sigma2_epsilon = Rcpp::as<double>(beta_sigma2_epsilon_prior);
  }

double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double sigma2_theta = 10000.00;
if(sigma2_theta_prior.isNotNull()){
  sigma2_theta = Rcpp::as<double>(sigma2_theta_prior);
  }

double a_lambda = 1.00;
if(a_lambda_prior.isNotNull()){
  a_lambda = Rcpp::as<double>(a_lambda_prior);
  }

double b_lambda = 1.00;
if(b_lambda_prior.isNotNull()){
  b_lambda = Rcpp::as<double>(b_lambda_prior);
  }

double alpha_sigma2_eta = 0.01;
if(alpha_sigma2_eta_prior.isNotNull()){
  alpha_sigma2_eta = Rcpp::as<double>(alpha_sigma2_eta_prior);
  }

double beta_sigma2_eta = 0.01;
if(beta_sigma2_eta_prior.isNotNull()){
  beta_sigma2_eta = Rcpp::as<double>(beta_sigma2_eta_prior);
  }
  
double alpha_phi_eta = 1.00;
if(alpha_phi_eta_prior.isNotNull()){
  alpha_phi_eta = Rcpp::as<double>(alpha_phi_eta_prior);
  }
  
double beta_phi_eta = 1.00;
if(beta_phi_eta_prior.isNotNull()){
  beta_phi_eta = Rcpp::as<double>(beta_phi_eta_prior);
  }

//Initial Values
double sigma2_epsilon = 1.00;
if(sigma2_epsilon_init.isNotNull()){
  sigma2_epsilon = Rcpp::as<double>(sigma2_epsilon_init);
  }

arma::vec beta(p_x); beta.fill(0.00);
if(beta_init.isNotNull()){
  beta = Rcpp::as<arma::vec>(beta_init);
  }

double theta = 1.00;
if(theta_init.isNotNull()){
  theta = Rcpp::as<double>(theta_init);
  }

double lambda = 0.01;
if(lambda_init.isNotNull()){
  lambda = Rcpp::as<double>(lambda_init);
  }

double sigma2_eta = 1.00;
if(sigma2_eta_init.isNotNull()){
  sigma2_eta = Rcpp::as<double>(sigma2_eta_init);
  }

double phi_eta = 1.00;
if(phi_eta_init.isNotNull()){
  phi_eta = Rcpp::as<double>(phi_eta_init);
  }

arma::mat eta_approx_corr_inv = inv_sympd(exp(-phi_eta*d_a_approx));
arma::vec eta_approx(k_approx); eta_approx.fill(0.00);
arma::vec eta = exp(-phi_eta*cross_a_approx)*(eta_approx_corr_inv*eta_approx);

arma::vec indicator = indicator_fun(h_model,
                                    distance_to_ps,
                                    angle_key,
                                    lambda,
                                    eta);

double neg_two_loglike = neg_two_loglike_update(likelihood_indicator,
                                                y,
                                                x,
                                                tri_als,
                                                indicator,
                                                sigma2_epsilon,
                                                beta,
                                                theta);

//Metropolis Settings
int acctot_lambda = 1;
arma::vec acctot_eta(k_approx); acctot_eta.fill(1);
int acctot_phi_eta = 1;

//Main Sampling Loop
arma::vec omega(n); omega.fill(0.00);
arma::vec kappa = y;
int counter = 0;

for(int j = 1; j < mcmc_samples; ++j){
   
   if(likelihood_indicator == 0){ //Bernoulli
  
     //omega Update
     Rcpp::List omega_output = omega_update(y,
                                            x,
                                            tri_als,
                                            indicator,
                                            beta,
                                            theta);
  
     omega = Rcpp::as<arma::vec>(omega_output[0]);
     kappa = Rcpp::as<arma::vec>(omega_output[1]);
  
     }
   
   if(likelihood_indicator == 1){ //Gaussian
      
     //sigma2_epsilon Update
     sigma2_epsilon = sigma2_epsilon_update(y,
                                            x,
                                            indicator,
                                            beta,
                                            theta,
                                            alpha_sigma2_epsilon,
                                            beta_sigma2_epsilon);
     omega.fill(1.00/sigma2_epsilon);
      
     }
   
   //beta Update
   beta = beta_update(x, 
                      indicator,
                      omega,
                      kappa,
                      theta,
                      sigma2_beta);
   
   //theta Update
   theta = theta_update(x, 
                        indicator,
                        omega,
                        kappa,
                        beta,
                        sigma2_theta);
   
   //lambda Update
   Rcpp::List lambda_output = lambda_update(h_model,
                                            distance_to_ps,
                                            angle_key,
                                            x,
                                            indicator,
                                            omega,
                                            kappa,
                                            beta,
                                            theta,
                                            lambda,
                                            eta,
                                            a_lambda,
                                            b_lambda,
                                            metrop_var_lambda,
                                            acctot_lambda);
   
   lambda = Rcpp::as<double>(lambda_output[0]);
   indicator = Rcpp::as<arma::vec>(lambda_output[1]);
   acctot_lambda = Rcpp::as<int>(lambda_output[2]);
   
   //eta Update
   Rcpp::List eta_output = eta_update(k_approx,
                                      cross_a_approx,
                                      eta_approx_corr_inv,
                                      h_model,
                                      distance_to_ps,
                                      angle_key,
                                      x,
                                      indicator,
                                      omega,
                                      kappa,
                                      beta,
                                      theta,
                                      lambda,
                                      eta,
                                      eta_approx,
                                      sigma2_eta,
                                      phi_eta,
                                      metrop_var_eta,
                                      acctot_eta);
   
   eta = Rcpp::as<arma::vec>(eta_output[0]);
   eta_approx = Rcpp::as<arma::vec>(eta_output[1]);
   indicator = Rcpp::as<arma::vec>(eta_output[2]);
   acctot_eta = Rcpp::as<arma::vec>(eta_output[3]);
   
   //sigma2_eta update
   sigma2_eta = sigma2_eta_update(eta_approx_corr_inv,
                                  eta_approx,
                                  alpha_sigma2_eta,
                                  beta_sigma2_eta);
   
   //phi_eta Update
   Rcpp::List phi_eta_output = phi_eta_update(d_a_approx,
                                              cross_a_approx,
                                              eta_approx_corr_inv,
                                              h_model,
                                              distance_to_ps,
                                              angle_key,
                                              x,
                                              indicator,
                                              omega,
                                              kappa,
                                              beta,
                                              theta,
                                              lambda,
                                              eta,
                                              eta_approx,
                                              sigma2_eta,
                                              phi_eta,
                                              alpha_phi_eta,
                                              beta_phi_eta,
                                              metrop_var_phi_eta,
                                              acctot_phi_eta);

   phi_eta = Rcpp::as<double>(phi_eta_output[0]);
   eta_approx_corr_inv = Rcpp::as<arma::mat>(phi_eta_output[1]);
   eta = Rcpp::as<arma::vec>(phi_eta_output[2]);
   indicator= Rcpp::as<arma::vec>(phi_eta_output[3]);
   acctot_phi_eta = Rcpp::as<int>(phi_eta_output[4]);
    
   //neg_two_loglike Update
   neg_two_loglike = neg_two_loglike_update(likelihood_indicator, 
                                            y,
                                            x,
                                            tri_als,
                                            indicator,
                                            sigma2_epsilon,
                                            beta,
                                            theta);
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
     
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
     
     if(h_model == 0){
       Rcpp::Rcout << "Indicator" << std::endl;
       }
     
     if(h_model == 1){
       Rcpp::Rcout << "Linear" << std::endl;
       }
     
     if(h_model == 2){
       Rcpp::Rcout << "Exponential" << std::endl;
       }
     
     if(h_model == 3){
       Rcpp::Rcout << "Gaussian" << std::endl;
       }
     
     if(h_model == 4){
       Rcpp::Rcout << "Spherical" << std::endl;
       }
     
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     
     double accrate_lambda = round(100*(acctot_lambda/(double)j));
     Rcpp::Rcout << "lambda Acceptance: " << accrate_lambda << "%" << std::endl;
     
     double accrate_eta_min = round(100*(min(acctot_eta)/(double)j));
     Rcpp::Rcout << "eta Acceptance (min): " << accrate_eta_min << "%" << std::endl;
     
     double accrate_eta_max = round(100*(max(acctot_eta)/(double)j));
     Rcpp::Rcout << "eta Acceptance (max): " << accrate_eta_max << "%" << std::endl;
     
     double accrate_phi_eta = round(100*(acctot_phi_eta/(double)j));
     Rcpp::Rcout << "phi_eta Acceptance: " << accrate_phi_eta << "%" << std::endl;
     
     Rcpp::Rcout << "*************************" << std::endl;
     
     }
   
   if((j + 1) == adapt){
      
     metrop_var_lambda = metrop_var_lambda*adapt_lambda;
     metrop_var_eta = metrop_var_eta*adapt_eta;
     metrop_var_phi_eta = metrop_var_phi_eta*adapt_phi_eta;
      
     }
   
   if((j >= burnin) & (((j - burnin)) % thin == 0)){
     
     sigma2_epsilon_keep(counter) = sigma2_epsilon;
     beta_keep.col(counter) = beta;
     theta_keep(counter) = theta;
     lambda_keep(counter) = lambda;
     eta_keep.col(counter) = eta;
     sigma2_eta_keep(counter) = sigma2_eta;
     phi_eta_keep(counter) = phi_eta;
     neg_two_loglike_keep(counter) = neg_two_loglike;
     
     counter = counter + 
               1;
     
     }
   
   }
  
return Rcpp::List::create(Rcpp::Named("sigma2_epsilon_keep") = sigma2_epsilon_keep,
                          Rcpp::Named("beta_keep") = beta_keep,
                          Rcpp::Named("theta_keep") = theta_keep,
                          Rcpp::Named("lambda_keep") = lambda_keep,
                          Rcpp::Named("eta_keep") = eta_keep,
                          Rcpp::Named("sigma2_eta_keep") = sigma2_eta_keep,
                          Rcpp::Named("phi_eta_keep") = phi_eta_keep,
                          Rcpp::Named("neg_two_loglike_keep") = neg_two_loglike_keep,
                          Rcpp::Named("acctot_lambda") = acctot_lambda,
                          Rcpp::Named("acctot_phi_eta") = acctot_phi_eta,
                          Rcpp::Named("acctot_eta") = acctot_eta);

}
