#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double tau2_update(int n,
                   double a_tau2,
                   double b_tau2,
                   arma::vec beta_true,
                   double mu,
                   arma::vec eta_star,
                   arma::mat Sigma_eta_star_inv,
                   arma::mat x){

double a_tau2_update = n/2.00 +
                       a_tau2;
arma::vec temp_vec = beta_true +
                     -mu +
                     -x*(Sigma_eta_star_inv*eta_star);
double b_tau2_update = 0.50*dot(temp_vec, temp_vec) +
                       b_tau2;
double tau2 = 1.00/R::rgamma(a_tau2_update,
                             (1.00/b_tau2_update));

return(tau2);

}





