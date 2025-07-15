#include "RcppArmadillo.h"
#include "EDMeta.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_eta_update(int m_star,
                         double a_sigma2_eta,
                         double b_sigma2_eta,
                         arma::vec eta_star,
                         arma::mat Sigma_eta_star_inv){

double a_sigma2_eta_update = m_star/2.00 +
                             a_sigma2_eta;
double b_sigma2_eta_update = 0.50*dot(eta_star, (Sigma_eta_star_inv*eta_star)) +
                             b_sigma2_eta;
double sigma2_eta = 1.00/R::rgamma(a_sigma2_eta_update,
                                   (1.00/b_sigma2_eta_update));

return(sigma2_eta);

}





