#ifndef __DVCP__
#define __DVCP__

arma::vec indicator_fun(int h_model,
                        arma::vec distance_to_ps,
                        arma::vec angle_key,
                        double lambda,
                        arma::vec eta);

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x,
                             arma::mat indicator,
                             arma::vec beta,
                             double theta,
                             double alpha_sigma2_epsilon,
                             double beta_sigma2_epsilon);

Rcpp::List omega_update(arma::vec y,
                        arma::mat x,
                        arma::vec indicator,
                        arma::vec beta,
                        double theta);

arma::vec beta_update(arma::mat x, 
                      arma::vec indicator,
                      arma::vec omega,
                      arma::vec kappa,
                      double theta,
                      double sigma2_beta);

double theta_update(arma::mat x, 
                    arma::vec indicator,
                    arma::vec omega,
                    arma::vec kappa,
                    arma::vec beta,
                    double sigma2_theta);

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
                         int acctot_lambda);

Rcpp::List theta_update(arma::mat x,
                        arma::mat z,
                        arma::vec distance_to_ps,
                        double theta,
                        arma::vec w_aux,
                        arma::vec gamma,
                        arma::vec beta,
                        double lambda,
                        arma::vec w,
                        arma::vec spillover_covar_temp,
                        int spillover_covar_def,
                        double a_theta,
                        double b_theta,
                        double metrop_var_theta_trans,
                        int acctot_theta_trans);

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
                      arma::vec acctot_eta);

double sigma2_eta_update(arma::mat eta_approx_corr_inv,
                         arma::vec eta_approx,
                         double alpha_sigma2_eta,
                         double beta_sigma2_eta);

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
                          int acctot_phi_eta);

double neg_two_loglike_update(int likelihood_indicator,
                              arma::vec y,
                              arma::mat x,
                              arma::vec indicator,
                              double sigma2_epsilon, 
                              arma::vec beta,
                              double theta);

Rcpp::List DVCP(int mcmc_samples,
                int burnin,
                int thin,
                int adapt,
                int likelihood_indicator, //0: Bernoulli, 1: Gaussian
                int h_model, //0: Indicator; 1: Linear, 2: Exponential; 3: Gaussian
                int k_approx,
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
                Rcpp::Nullable<double> alpha_sigma2_epsilon_prior,
                Rcpp::Nullable<double> beta_sigma2_epsilon_prior,
                Rcpp::Nullable<double> sigma2_beta_prior,
                Rcpp::Nullable<double> sigma2_theta_prior,
                Rcpp::Nullable<double> a_lambda_prior,
                Rcpp::Nullable<double> b_lambda_prior,
                Rcpp::Nullable<double> alpha_sigma2_eta_prior,
                Rcpp::Nullable<double> beta_sigma2_eta_prior,
                Rcpp::Nullable<double> alpha_phi_eta_prior,
                Rcpp::Nullable<double> beta_phi_eta_prior,
                Rcpp::Nullable<double> sigma2_epsilon_init,
                Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                Rcpp::Nullable<double> theta_init,
                Rcpp::Nullable<double> lambda_init,
                Rcpp::Nullable<double> sigma2_eta_init,
                Rcpp::Nullable<double> phi_eta_init); 

#endif // __DVCP__
