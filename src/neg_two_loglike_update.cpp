#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update(arma::vec beta_hat,
                              arma::vec se,
                              int n,
                              double mu,
                              double tau2,
                              arma::vec eta_star,
                              arma::mat Sigma_eta_star_inv,
                              arma::mat x){
  
arma::vec dens(n); dens.fill(0.00);
  
arma::vec mean = mu +
                 x*(Sigma_eta_star_inv*eta_star);
  
for(int j = 0; j < n; ++j){
   dens(j) = R::dnorm(beta_hat(j),
                      mean(j),
                      sqrt(tau2 + pow(se(j), 2.00)),
                      TRUE);
   }
  
double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;
  
}


















