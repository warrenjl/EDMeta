#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

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
                          int acctot_phi_eta_trans){

arma::mat Sigma_eta_star_inv_old = Sigma_eta_star_inv;
double Sigma_eta_star_inv_log_deter_old = Sigma_eta_star_inv_log_deter;
arma::mat x_old = x;
double phi_eta_trans_old = log(phi_eta_old);

double denom = -0.50*dot((beta_true - mu - x_old*(Sigma_eta_star_inv_old*eta_star)), (beta_true - mu - x_old*(Sigma_eta_star_inv_old*eta_star)))/tau2 +
               0.50*Sigma_eta_star_inv_log_deter_old +
               -0.50*dot(eta_star, (Sigma_eta_star_inv_old*eta_star))/sigma2_eta +
               a_phi_eta*phi_eta_trans_old +
               -b_phi_eta*exp(phi_eta_trans_old);
      
double phi_eta_trans = R::rnorm(phi_eta_trans_old, 
                                metrop_sd_phi_eta_trans);
double phi_eta = exp(phi_eta_trans);
Sigma_eta_star_inv = inv_sympd(exp(-phi_eta*abs_d_star_diffs));
double sign = 0.00;     
log_det(Sigma_eta_star_inv_log_deter, sign, Sigma_eta_star_inv);

for(int j = 0; j < n; ++j){
  
   arma::rowvec weighted_sum = sum(exposure_dens[j]%arma::exp(-phi_eta*abs_d_diffs[j]), 0)*(max_exposure(j) - min_exposure(j))/m;
   x.row(j) = weighted_sum;
  
   } 
      
double numer = -0.50*dot((beta_true - mu - x*(Sigma_eta_star_inv*eta_star)), (beta_true - mu - x*(Sigma_eta_star_inv*eta_star)))/tau2 +
               0.50*Sigma_eta_star_inv_log_deter +
               -0.50*dot(eta_star, (Sigma_eta_star_inv*eta_star))/sigma2_eta +
               a_phi_eta*phi_eta_trans +
               -b_phi_eta*exp(phi_eta_trans);
        
int accept = 1;
double ratio = exp(numer - denom);
if(ratio < R::runif(0.00, 1.00)){
          
  phi_eta = phi_eta_old;
  Sigma_eta_star_inv = Sigma_eta_star_inv_old;
  Sigma_eta_star_inv_log_deter = Sigma_eta_star_inv_log_deter_old;
  x = x_old;
  accept = 0;
          
  }
acctot_phi_eta_trans = acctot_phi_eta_trans + 
                       accept;   

return Rcpp::List::create(Rcpp::Named("phi_eta") = phi_eta,
                          Rcpp::Named("acctot_phi_eta_trans") = acctot_phi_eta_trans,
                          Rcpp::Named("Sigma_eta_star_inv") = Sigma_eta_star_inv,
                          Rcpp::Named("Sigma_eta_star_inv_log_deter") = Sigma_eta_star_inv_log_deter,
                          Rcpp::Named("x") = x);

}

