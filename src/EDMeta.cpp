#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List EDMeta(int mcmc_samples,
                  arma::vec beta_hat,
                  arma::vec se,
                  arma::vec min_exposure,
                  arma::vec max_exposure,
                  arma::vec mid_exposure,
                  arma::vec sd_exposure,
                  int exposure_distribution,                      //0: Uniform; 1: Gaussian
                  double metrop_sd_phi_eta_trans,
                  Rcpp::Nullable<int> m_approx = R_NilValue,      //Integral Approximation
                  Rcpp::Nullable<int> m_star_approx = R_NilValue, //Predictive Process Approximation
                  Rcpp::Nullable<double> sigma2_mu_prior = R_NilValue,
                  Rcpp::Nullable<double> a_tau2_prior = R_NilValue,
                  Rcpp::Nullable<double> b_tau2_prior = R_NilValue,
                  Rcpp::Nullable<double> a_sigma2_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> b_sigma2_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_phi_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> b_phi_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> mu_init = R_NilValue,
                  Rcpp::Nullable<arma::vec> eta_star_init = R_NilValue,
                  Rcpp::Nullable<double> tau2_init = R_NilValue,
                  Rcpp::Nullable<double> sigma2_eta_init = R_NilValue,
                  Rcpp::Nullable<double> phi_eta_init = R_NilValue){
 
//Defining Parameters and Quantities of Interest
int n = beta_hat.size();
  
arma::vec a = min_exposure;
arma::vec b = max_exposure;
  
double a_star = min(min_exposure);
double b_star = max(max_exposure);

//Approximation Information
int m = 100;
if(m_approx.isNotNull()){
  m = Rcpp::as<int>(m_approx);
  }

int m_star = 100;
if(m_star_approx.isNotNull()){
  m_star = Rcpp::as<int>(m_star_approx);
  }

arma::mat d(n, m); d.fill(0.00);
std::vector<arma::mat> exposure_dens(n);
for(int j = 0; j < n; ++j){
  
   d.row(j) = linspace<rowvec>(min_exposure(j), 
                               max_exposure(j), 
                               m);
  
   if(exposure_distribution == 0){  //Uniform distribution
    
     arma::mat dens(m, m_star); dens.fill(0.00);
     for(int k = 0; k < m; ++k){
        for(int l = 0; l < m_star; ++l){
           dens(k, l) = 1.00/(max_exposure(j) - min_exposure(j));
           }
        }
     exposure_dens[j] = dens;
     
     }
    
   if(exposure_distribution == 1){  //Gaussian distribution  
     
     arma::mat dens(m, m_star); dens.fill(0.00);
     for(int k = 0; k < m; ++k){
       
        double val = d(j, k);
        for(int l = 0; l < m_star; ++l){
           dens(k, l) = R::dnorm(val, 
                                 mid_exposure(j), 
                                 sd_exposure(j), 
                                 false);
           }
        
        }
     exposure_dens[j] = dens;
    
     }
   
   }
  
arma::vec d_star = linspace<vec>(a_star, 
                                 b_star, 
                                 m_star);
std::vector<arma::mat> abs_d_diffs(n);
for(int j = 0; j < n; ++j){
  
   arma::mat diffs(m, m_star); diffs.fill(0.00);
   for(int k = 0; k < m; ++k){
      for(int l = 0; l < m_star; ++l){
         diffs(k, l) = std::abs(d(j, k) - d_star(l))/b_star;
         }
      }
   abs_d_diffs[j] = diffs;
   
   }
arma::mat abs_d_star_diffs = arma::abs(arma::repmat(d_star, 1, m_star) - arma::repmat(d_star.t(), m_star, 1))/b_star;

arma::mat inv_var_mat(n, n); inv_var_mat.fill(0.00);
for(int j = 0; j < n; ++ j){
   inv_var_mat(j,j) = (1.00/pow(se(j), 2)); 
   }
  
arma::mat beta_true(n, mcmc_samples); beta_true.fill(0.00);
arma::vec mu(mcmc_samples); mu.fill(0.00);
arma::mat eta_star(m_star, mcmc_samples); eta_star.fill(0.00);
arma::vec tau2(mcmc_samples); tau2.fill(0.00);
arma::vec sigma2_eta(mcmc_samples); sigma2_eta.fill(0.00);
arma::vec phi_eta(mcmc_samples); phi_eta.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_mu = 10000.00;
if(sigma2_mu_prior.isNotNull()){
  sigma2_mu = Rcpp::as<double>(sigma2_mu_prior);
  }
double sigma2_mu_inv = 1.00/sigma2_mu;

double a_tau2 = 0.01;
if(a_tau2_prior.isNotNull()){
  a_tau2 = Rcpp::as<double>(a_tau2_prior);
  }

double b_tau2 = 0.01;
if(b_tau2_prior.isNotNull()){
  b_tau2 = Rcpp::as<double>(b_tau2_prior);
  }

double a_sigma2_eta = 0.01;
if(a_sigma2_eta_prior.isNotNull()){
  a_sigma2_eta = Rcpp::as<double>(a_sigma2_eta_prior);
  }

double b_sigma2_eta = 0.01;
if(b_sigma2_eta_prior.isNotNull()){
  b_sigma2_eta = Rcpp::as<double>(b_sigma2_eta_prior);
  }

double a_phi_eta = 1.00;
if(a_phi_eta_prior.isNotNull()){
  a_phi_eta = Rcpp::as<double>(a_phi_eta_prior);
  }

double b_phi_eta = 1.00;
if(b_phi_eta_prior.isNotNull()){
  b_phi_eta = Rcpp::as<double>(b_phi_eta_prior);
  }

//Initial Values
beta_true.col(0) = beta_hat;

mu(0) = 0.00;
if(mu_init.isNotNull()){
  mu(0) = Rcpp::as<double>(mu_init);
  }

eta_star.col(0).fill(0.00);
if(eta_star_init.isNotNull()){
  eta_star.col(0) = Rcpp::as<arma::vec>(eta_star_init);
  }

tau2(0) = 1.00;
if(tau2_init.isNotNull()){
  tau2(0) = Rcpp::as<double>(tau2_init);
  }

sigma2_eta(0) = 1.00;
if(sigma2_eta_init.isNotNull()){
  sigma2_eta(0) = Rcpp::as<double>(sigma2_eta_init);
  }

phi_eta(0) = 1.00;
if(phi_eta_init.isNotNull()){
  phi_eta(0) = Rcpp::as<double>(phi_eta_init);
  }

arma::mat Sigma_eta_star_inv = inv_sympd(exp(-phi_eta(0)*abs_d_star_diffs));

double Sigma_eta_star_inv_log_deter = 0.00; 
double sign = 0.00;     
log_det(Sigma_eta_star_inv_log_deter, sign, Sigma_eta_star_inv);

arma::mat x(n, m_star); x.fill(0.00);
for(int j = 0; j < n; ++j){
  
   arma::rowvec weighted_sum = sum(exposure_dens[j]%arma::exp(-phi_eta(0)*abs_d_diffs[j]), 0)*(max_exposure(j) - min_exposure(j))/m;
   x.row(j) = weighted_sum;
   
   }

neg_two_loglike(0) = neg_two_loglike_update(beta_hat,
                                            se,
                                            n,
                                            mu(0),
                                            tau2(0),
                                            eta_star.col(0),
                                            Sigma_eta_star_inv,
                                            x);
//Metropolis Settings
int acctot_phi_eta_trans = 0;

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){
   
   //beta_true
   beta_true.col(j) = beta_true_update(beta_hat,
                                       inv_var_mat,
                                       n,
                                       mu(j-1),
                                       eta_star.col(j-1),
                                       tau2(j-1),
                                       Sigma_eta_star_inv,
                                       x);
   
   //mu
   mu(j) = mu_update(n,
                     sigma2_mu_inv,
                     beta_true.col(j),
                     eta_star.col(j-1),
                     tau2(j-1),
                     Sigma_eta_star_inv,
                     x);
   
   //eta_star
   eta_star.col(j) = eta_star_update(beta_hat,
                                     m_star,
                                     beta_true.col(j),
                                     mu(j),
                                     tau2(j-1),
                                     sigma2_eta(j-1),
                                     Sigma_eta_star_inv,
                                     x);
   
   //tau2
   tau2(j) = tau2_update(n,
                         a_tau2,
                         b_tau2,
                         beta_true.col(j),
                         mu(j),
                         eta_star.col(j),
                         Sigma_eta_star_inv,
                         x);
     
   //sigma2_eta Update
   sigma2_eta(j) = sigma2_eta_update(m_star,
                                     a_sigma2_eta,
                                     b_sigma2_eta,
                                     eta_star.col(j),
                                     Sigma_eta_star_inv);
      
   //phi_eta Update
   Rcpp::List phi_eta_output = phi_eta_update(n,
                                              m,
                                              min_exposure,
                                              max_exposure,
                                              abs_d_diffs,
                                              abs_d_star_diffs,
                                              exposure_dens,
                                              a_phi_eta,
                                              b_phi_eta,
                                              beta_true.col(j),
                                              mu(j),
                                              tau2(j),
                                              eta_star.col(j),
                                              sigma2_eta(j),
                                              phi_eta(j-1),
                                              Sigma_eta_star_inv,
                                              x,
                                              Sigma_eta_star_inv_log_deter,
                                              metrop_sd_phi_eta_trans,
                                              acctot_phi_eta_trans);
   phi_eta(j) = Rcpp::as<double>(phi_eta_output[0]);
   acctot_phi_eta_trans = Rcpp::as<int>(phi_eta_output[1]);
   Sigma_eta_star_inv = Rcpp::as<arma::mat>(phi_eta_output[2]);
   Sigma_eta_star_inv_log_deter = Rcpp::as<double>(phi_eta_output[3]);
   x = Rcpp::as<arma::mat>(phi_eta_output[4]);
     
   //neg_two_loglike
   neg_two_loglike(j) = neg_two_loglike_update(beta_hat,
                                               se,
                                               n,
                                               mu(j),
                                               tau2(j),
                                               eta_star.col(j),
                                               Sigma_eta_star_inv,
                                               x);
      
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     
     double accrate_phi_eta_trans = round(100*(acctot_phi_eta_trans/(double)j));
     Rcpp::Rcout << "phi_eta Acceptance: " << accrate_phi_eta_trans << "%" << std::endl;
     
     Rcpp::Rcout << "************************" << std::endl;
    
     }
  
   }
                                  
  return Rcpp::List::create(Rcpp::Named("beta_true") = beta_true,
                            Rcpp::Named("mu") = mu,
                            Rcpp::Named("eta_star") = eta_star,
                            Rcpp::Named("tau2") = tau2,
                            Rcpp::Named("sigma2_eta") = sigma2_eta,
                            Rcpp::Named("phi_eta") = phi_eta,
                            Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                            Rcpp::Named("acctot_phi_eta_trans") = acctot_phi_eta_trans);
   
  }


