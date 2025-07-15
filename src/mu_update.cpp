#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double mu_update(int n,
                 double sigma2_mu_inv,
                 arma::vec beta_true,
                 arma::vec eta_star_old,
                 double tau2_old,
                 arma::mat Sigma_eta_star_inv,
                 arma::mat x){
  
double var_mu = 1.00/((n/tau2_old) + sigma2_mu_inv);
double mean_mu = var_mu*sum(beta_true - x*(Sigma_eta_star_inv*eta_star_old))/tau2_old;
double mu = mean_mu +
            arma::randn()*sqrt(var_mu);
  
return(mu);
  
}


