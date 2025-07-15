#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_true_update(arma::vec beta_hat,
                           arma::mat inv_var_mat,
                           int n,
                           double mu_old,
                           arma::vec eta_star_old,
                           double tau2_old,
                           arma::mat Sigma_eta_star_inv,
                           arma::mat x){
  
arma::mat cov_beta_true = inv_sympd(inv_var_mat + 
                                    (1.00/tau2_old)*eye(n, n));
arma::mat mean_beta_true = cov_beta_true*(inv_var_mat*beta_hat + 
                                          (mu_old + x*(Sigma_eta_star_inv*eta_star_old))/tau2_old);
arma::mat ind_norms = arma::randn(1, n);
arma::vec beta_true = mean_beta_true + 
                      trans(ind_norms*arma::chol(cov_beta_true));

return(beta_true);

}






