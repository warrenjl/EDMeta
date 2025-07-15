#ifndef __EDMeta__
#define __EDMeta__

double neg_two_loglike_update(arma::vec beta_hat,
                              arma::vec se,
                              int n,
                              double mu,
                              double tau2,
                              arma::vec eta_star,
                              arma::mat Sigma_eta_star_inv,
                              arma::mat x);

arma::vec beta_true_update(arma::vec beta_hat,
                           arma::mat inv_var_mat,
                           int n,
                           double mu_old,
                           arma::vec eta_star_old,
                           double tau2_old,
                           arma::mat Sigma_eta_star_inv,
                           arma::mat x);

double mu_update(int n,
                 double sigma2_mu_inv,
                 arma::vec beta_true,
                 arma::vec eta_star_old,
                 double tau2_old,
                 arma::mat Sigma_eta_star_inv,
                 arma::mat x);
  
arma::vec eta_star_update(arma::vec beta_hat,
                          int m_star,
                          arma::vec beta_true,
                          double mu,
                          double tau2_old,
                          double sigma2_eta_old,
                          arma::mat Sigma_eta_star_inv,
                          arma::mat x);
  
double tau2_update(int n,
                   double a_tau2,
                   double b_tau2,
                   arma::vec beta_true,
                   double mu,
                   arma::vec eta_star,
                   arma::mat Sigma_eta_star_inv,
                   arma::mat x);

double sigma2_eta_update(int m_star,
                         double a_sigma2_eta,
                         double b_sigma2_eta,
                         arma::vec eta_star,
                         arma::mat Sigma_eta_star_inv);

Rcpp::List phi_eta_update(int n,
                          int m,
                          arma::vec min_exposure,
                          arma::vec max_exposure,
                          std::vector<arma::mat> abs_d_diffs,
                          arma::mat abs_d_star_diffs,
                          std::vector<arma::mat> exposure_dens,
                          double a_phi_eta,
                          double b_phi_eta,
                          arma::vec beta_true,
                          double mu,
                          double tau2,
                          arma::vec eta_star,
                          double sigma2_eta,
                          double phi_eta_old,
                          arma::mat Sigma_eta_star_inv,
                          arma::mat x,
                          double Sigma_eta_star_inv_log_deter,
                          double metrop_sd_phi_eta_trans,
                          int acctot_phi_eta_trans);

Rcpp::List EDMeta(int mcmc_samples,
                  arma::vec beta_hat,
                  arma::vec se,
                  arma::vec min_exposure,
                  arma::vec max_exposure,
                  arma::vec mid_exposure,
                  arma::vec sd_exposure,
                  int exposure_distribution,         //0: Uniform; 1: Gaussian
                  double metrop_sd_phi_eta_trans,
                  Rcpp::Nullable<int> m_approx,      //Integral Approximation
                  Rcpp::Nullable<int> m_star_approx, //Predictive Process Approximation
                  Rcpp::Nullable<double> sigma2_mu_prior,
                  Rcpp::Nullable<double> a_tau2_prior,
                  Rcpp::Nullable<double> b_tau2_prior,
                  Rcpp::Nullable<double> a_sigma2_eta_prior,
                  Rcpp::Nullable<double> b_sigma2_eta_prior,
                  Rcpp::Nullable<double> a_phi_eta_prior,
                  Rcpp::Nullable<double> b_phi_eta_prior,
                  Rcpp::Nullable<double> mu_init,
                  Rcpp::Nullable<arma::vec> eta_star_init,
                  Rcpp::Nullable<double> tau2_init,
                  Rcpp::Nullable<double> sigma2_eta_init,
                  Rcpp::Nullable<double> phi_eta_init);

#endif // __EDMeta__
