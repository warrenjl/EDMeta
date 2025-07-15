#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_star_update(arma::vec beta_hat,
                          int m_star,
                          arma::vec beta_true,
                          double mu,
                          double tau2_old,
                          double sigma2_eta_old,
                          arma::mat Sigma_eta_star_inv,
                          arma::mat x){
  
arma::mat temp_mat = x*Sigma_eta_star_inv;
arma::mat cov_eta_star = inv_sympd((temp_mat.t()*temp_mat/tau2_old) +
                                   (Sigma_eta_star_inv/sigma2_eta_old));
arma::vec mean_eta_star = cov_eta_star*(temp_mat.t()*(beta_true - mu))/tau2_old;
  
arma::mat ind_norms = arma::randn(1, m_star);
arma::vec eta_star = mean_eta_star + 
                     trans(ind_norms*arma::chol(cov_eta_star));

return(eta_star);

}






